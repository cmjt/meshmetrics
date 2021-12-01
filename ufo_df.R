library(spatstat)
library(INLA)
library(maps) ## Display of maps
library(ggplot2)
library(sp)
library(patchwork)
## read in data
ufo <- readr::read_csv("data/ufo.csv")
## source required functions
source("functions.r")
## define region and points
region <-  as(sf::st_as_sf(maps::map("usa", fill = TRUE, plot = FALSE)), "Spatial")
locs <- data.frame(x = ufo$city_longitude, y = ufo$city_latitude)

## convert the existing lat-long coordinates to UTM (easting-northing system)
ufo_locs <- ufo
coordinates(ufo_locs) <- c("city_longitude", "city_latitude")
proj4string(ufo_locs) <- CRS("+proj=longlat +datum=WGS84")
ufo_sp <- spTransform(ufo_locs, CRS("+init=esri:102003"))

## create an interactive map using use a basemap given by the OpenStreetmap provider
library(tmap)
tmap_mode("view")
tm_basemap(leaflet::providers$OpenStreetMap) +
  tm_shape(cord.dec) + tm_dots()


## Display of maps
library(maps) 
states <- map_data("state")
usa <- map_data('usa')

ggplot() +
  # geom_polygon(data=states, aes(x=long, y=lat, group=group),fill="white", colour="black") +
  geom_path(data=usa, aes(x=long, y=lat, group=group)) +
  geom_point(data = ufo, aes(x=city_longitude, y=city_latitude), alpha=0.5) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()) +
  coord_map("mercator")


## mesh examples (first three obvioulsy ridiculous mesh)

mesh01 <- inla.mesh.2d(boundary = inla.sp2segment(region),
                       max.edge = c(7.5,15))

## add locs
mesh02 <- inla.mesh.2d(loc = locs, boundary = inla.sp2segment(region),
                       max.edge = c(7.5,15))

## add cutoff
mesh03 <- inla.mesh.2d(boundary = inla.sp2segment(region),
                       max.edge = c(7.5,15), cutoff = 10)

## make edges smaller
mesh04 <- inla.mesh.2d(loc = locs, boundary = inla.sp2segment(region),
                       max.edge = c(1.5,2), cutoff = 3)

## bigger difference internally
mesh05 <- inla.mesh.2d(loc = locs, boundary = inla.sp2segment(region),
                       max.edge = c(0.75,2), cutoff = 0.5)

## fine scale
mesh06 <- inla.mesh.2d(loc = locs, boundary = inla.sp2segment(region),
                       max.edge = c(0.75,1), cutoff = 0.5)


## Use mesh05 to fit the model
plot(mesh05, asp=1)


## Number of vertices in the mesh05
(nv <- mesh05$n)
## number of observations
(n <- nrow(ufo_sp))

## Create SPDE
ufo_spde <- inla.spde2.pcmatern(mesh = mesh05, alpha = 2,
                                prior.range = c(0.1, 0.01), # Pro(range < 0.1) = 0.01
                                prior.sigma = c(0.01, 0.01)) # Pro(sigma > 0.01) = 0.01


## mesh weights
source("book.mesh.dual.R")
dmesh <- book.mesh.dual(mesh05)
usa_bdy <- SpatialPoints(region, proj4string = CRS("+proj=longlat +datum=WGS84"))
usa_UTM <- spTransform(usa_bdy, CRS("+init=esri:102003"))
domainSP <- SpatialPolygons(list(Polygons(list(Polygon(usa_bdy)), '0')))
library(rgeos)
w <- sapply(1:length(dmesh), function(i) {
  if (gIntersects(dmesh[i, ], domainSP))
    return(gArea(gIntersection(dmesh[i, ], domainSP)))
  else return(0)
})

sum(w)
table(w > 0)


## projection matrix

y.pp <- rep(0:1, c(nv, n))
e.pp <- c(w, rep(0, n))
imat <- Diagonal(nv, rep(1, nv))
lmat <- inla.spde.make.A(mesh05, coordinates(ufo_sp))
A.pp <- rbind(imat, lmat)

stk.pp <- inla.stack(data = list(y = y.pp, e = e.pp),
                     A = list(1, A.pp),
                     effects = list(list(b0 = rep(1, nv + n)), list(i = 1:nv)),
                     tag = 'pp')


## model without covariates
pp.res <- inla(y ~ 0 + b0 + f(i, model = ufo_spde),
               family = 'poisson', data = inla.stack.data(stk.pp),
               control.predictor = list(A = inla.stack.A(stk.pp)),
               control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
               E = inla.stack.data(stk.pp)$e)
summary(pp.res)
pp.res$summary.fixed
pp.res$summary.hyperpar



## inlabru for model fitting
library(inlabru)
sp <- SpatialPolygons(list(Polygons(list(Polygon(usa_bdy)), '0')))
ufo_pts <- data.frame(x=ufo[, "city_longitude"], y=ufo[, "city_latitude"])
coordinates(ufo_pts) <- c("x", "y")


## Construct latent model components

matern <- inla.spde2.pcmatern(mesh05,
                              prior.sigma = c(0.1, 0.01),
                              prior.range = c(0.01, 0.01))

cmp <- coordinates ~ mySmooth(main = coordinates, model = matern) + Intercept

fits <- lgcp(cmp, ufo_pts,
             samplers = sp,
             domain = list(coordinates = mesh05),
             options = list(control.inla = list(int.strategy = "eb")))

fits$summary.fixed
fits$summary.hyperpar


### Model with Population as covariates, 2019 population





