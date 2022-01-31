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
ufo_sp <- spTransform(ufo_locs, CRS("+init=epsg:2163"))


# US region / domain
# Conversion into SpatialPoints
usa_bdy <- usa_region
coordinates(usa_bdy) <- c("Longitude", "Latitude")
# Setting default projection
proj4string(usa_bdy) <- CRS('+init=epsg:4326')
usa_utm <- spTransform(usa_bdy, CRS("+init=epsg:2163"))


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


source("book.mesh.dual.R")
dmesh <- book.mesh.dual(mesh05)
usabdy.sp <- SpatialPolygons(list(Polygons(list(Polygon(usa_utm)), ID = "1")))
library(rgeos)
usabdy.sp <- gBuffer(usabdy.sp, byid=TRUE, width=0)
w <- sapply(1:length(dmesh), function(i) {
  if (gIntersects(dmesh[i, ], usabdy.sp)) 
    return(gArea(gIntersection(dmesh[i, ], usabdy.sp)))
  else 
    return(0) 
})

sum(w)
table(w > 0)


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
library(sf)
library(dplyr)
### US 2020 Cartographic Boundary Files
### downloaded from https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_us_county_20m.zip
US20 <- st_read("data/cb_2020_us_county_20m/cb_2020_us_county_20m.shp")
# leave out AK, HI, and PR (state FIPS: 02, 15, and 72)
conti_US20 <- US20[!(US20$STATEFP %in% c("02","15","72")), ]


### 2020 Census Demographic Data By County, downloaded from https://www.dataemporium.com/dataset/254/
census20 <- read.csv("data/US_2020_census.csv")
contiguous <- subset(census20, census20$STATE_ABBR != "AK" & census20$STATE_ABBR != "HI" & census20$STATE_ABBR != "PR")
pop20 <- contiguous[, 1:3]

### append total population to US20
conti_US20 <- arrange(conti_US20, NAMELSAD)
pop20 <- arrange(pop20, COUNTY)
conti_US20$TOT_POP <- pop20$TOTAL_POPULATION



### create shape object with state polygons
us_pop20 <- conti_US20
# us_pop20 <- unionSpatialPolygons(as(us_pop20, "Spatial"), IDs = us_pop20$STATEFP)
us_pop20 <- as(us_pop20, "Spatial")
us_pop20 <- spTransform(us_pop20, CRS("+init=epsg:2163"))
pop_mesh <- sp::over(SpatialPoints(mesh05$loc[,1:2], 
                                   proj4string = CRS(proj4string(ufo_sp))), us_pop20)$TOT_POP
pop_obs <- sp::over(ufo_sp, us_pop20)$TOT_POP
covs <- data.frame(pop = c(pop_mesh, pop_obs))

## data stack
stk.cov <- inla.stack(
  data = list(y = y.pp, e = e.pp), 
  A = list(1, A.pp),
  effects = list(list(b0 = rep(1, nv + n), pop = covs$pop), list(i = 1:nv)),
  tag = 'pp')

## fit model with a covariate
pp.cov <- inla(y ~ 0 + b0 + pop + f(i, model = ufo_spde), 
               family = 'poisson', data = inla.stack.data(stk.cov), 
               control.predictor = list(A = inla.stack.A(stk.cov)), 
               E = inla.stack.data(stk.pp)$e)

## coefficients of the fixed effects
pp.cov$summary.fixed
pp.cov$summary.hyperpar


