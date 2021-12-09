## packages
library(sf)
library(maps)
library(sp)
library(INLA)
require(ggplot2)
require(patchwork) ## just for aligning plots
## source required functions
source("functions.r")

## read in data
ufo <- readr::read_csv("data/ufo.csv")
## usa map
states <- map_data("state")
usa <- map_data('usa')
usa_region <- data.frame(Longitude = usa$long, Latitude = usa$lat)

## define region and points
region <- as(sf::st_as_sf(maps::map("usa", fill = TRUE, plot = FALSE)), "Spatial")
locs <- data.frame(x = ufo$city_longitude, y = ufo$city_latitude)

## convert the existing lat-long coordinates to UTM (easting-northing system)
ufo_locs <- ufo
# Conversion into SpatialPoints
coordinates(ufo_locs) <- c("city_longitude", "city_latitude")
# Setting default projection
proj4string(ufo_locs) <- CRS('+init=epsg:4326')
ufo_sp <- spTransform(ufo_locs, CRS("+init=epsg:5070"))


# US region / domain
# Conversion into SpatialPoints
usa_bdy <- usa_region
coordinates(usa_bdy) <- c("Longitude", "Latitude")
# Setting default projection
proj4string(usa_bdy) <- CRS('+init=epsg:4326')
usa_utm <- spTransform(usa_bdy, CRS('+init=epsg:5070'))


mesh05 <- inla.mesh.2d(loc = coordinates(ufo_sp), max.edge = c(150000, 500000), cutoff = 5000)

plot(mesh05, asp=1)

## number of vertices in the mesh05
(nv <- mesh05$n)
## number of observations
(n <- nrow(ufo))

## Create SPDE
ufo_spde <- inla.spde2.pcmatern(mesh = mesh05, 
                                alpha = 2,
                                prior.range = c(0.1, 0.01), # P(range < 0.1) = 0.01
                                prior.sigma = c(0.01, 0.01)) # P(sigma > 0.01) = 0.01


library("deldir")
library("SDraw")
library("rgeos")
mytiles <- voronoi.polygons(SpatialPoints(mesh05$loc[, 1:2]))
usabdy.sp <- SpatialPolygons(list(Polygons(list(Polygon(usa_utm)), ID = "1")))
usabdy.sp <- gBuffer(usabdy.sp, byid=TRUE, width=0)
## Compute weights
w <- sapply(1:length(mytiles), function(p) {
  aux <- mytiles[p, ]  
  if(gIntersects(aux, usabdy.sp) ) {
    return(gArea(gIntersection(aux, usabdy.sp)))
  } else {
    return(0)
  }
})

sum(w)
table(w>0)

plot(mytiles, border ="gray")
points(mesh05$loc[, 1:2], pch = 19, cex = 0.25)
lines(coordinates(usa_utm), col = "green", lwd = 2)



## projection matrix
y.pp <- rep(0:1, c(nv, n))
e.pp <- c(w, rep(0, n))
imat <- Diagonal(nv, rep(1, nv))
lmat <- inla.spde.make.A(mesh05, coordinates(ufo_sp))
A.pp <- rbind(imat, lmat)


## the data stack
stk.pp <- inla.stack(
  data = list(y = y.pp, e = e.pp),
  A = list(1, A.pp),
  effects = list(list(b0 = rep(1, nv + n)), list(i = 1:nv)),
  tag = 'pp')


## model with no covariate
pp.res <- inla(y ~ 0 + b0 + f(i, model = ufo_spde),
               family = 'poisson', 
               data = inla.stack.data(stk.pp),
               control.predictor = list(A = inla.stack.A(stk.pp)),
               control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
               E = inla.stack.data(stk.pp)$e)

summary(pp.res)
pp.res$summary.fixed
pp.res$summary.hyperpar





## Adding a covariate, population

library(raster)
library(rgdal)
GDALinfo("pden2010_block/pden2010_60m.tif")
pop_ras <- raster("pden2010_block/pden2010_60m.tif")
pop_ras
library(inlmisc)
spdf_1 <- Grid2Polygons(pop_ras)
spdf_2 <- as(pop_ras,'SpatialPolygonsDataFrame')


library(stars)
r <- stack("pden2010_block/pden2010_60m.tif")
plot(r)



## data stack
stk.cov <- inla.stack(
  data = list(y = y.pp, e = e.pp), 
  A = list(1, A.pp),
  effects = list(list(b0 = rep(1, nv + n), pop = covs$pop), list(i = 1:nv)),
  tag = 'pp')

## fit model with a covariate
pp.cov <- inla(y ~ 0 + b0 + pop + f(i, model = spde), 
               family = 'poisson', data = inla.stack.data(stk.cov), 
               control.predictor = list(A = inla.stack.A(stk.cov)), 
               E = inla.stack.data(stk.pp)$e)

## coefficients of the fixed effects
pp.cov$summary.fixed

idx <- inla.stack.index(stk.pp, 'pred')$data
## Model with no covariates
spdf$SPDE0 <- pp.res$summary.fitted.values[idx, "mean"]
## Model with covariates
spdf$SPDE <- pp.cov$summary.fitted.values[idx, "mean"]

spplot(spdf, c("SPDE0", "SPDE"), 
       names.attr = c("No covariates", "With covariates"),
       col.regions = viridis::plasma(16))







