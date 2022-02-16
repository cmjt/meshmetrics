## packages
library(sf)
library(maps)
library(sp)
library(INLA)
require(ggplot2)
require(patchwork) ## just for aligning plots
## source required functionsR
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


### --------------------------------------------------------------------------->> Mesh Construction 
### We compare six results for the ufo dataset based on the six different meshes
### mesh_val <- c(loc, boundary, max.edge, cutoff)
source("mesh_sim.R")
mesh_val <- c(0, "region", c(7.5,15), 0,
              "locs", "region", c(7.5,15), 0,
              0, "region", c(7.5,15), 10,
              "locs", "region", c(1.5,2), 3,
              "locs", "region", c(0.75,2), 0.5,
              "locs", "region", c(0.75,1), 0.5)
mesh_mat <- matrix(mesh_val, ncol=5, byrow = TRUE)
mesh_n <- mesh_sim(mesh_mat)



### --------------------------------------------------------------------------->> Dual Mesh and Weights
temp <- vector("list", length = nrow(mesh_mat))
n <- c() # number of observations
nv <- c() # number of vertices in the mesh
dual_mesh <- temp
mesh_weight <- temp
ufo_spde <- temp

for (i in 1:nrow(mesh_mat)) {
  # number of observations
  n <- nrow(ufo)
  # number of vertices in the mesh
  nv[i] <- mesh_n[[i]]$n
  # create the dual mesh polygons
  source("functions.r")
  dual_mesh[[i]] <- book.mesh.dual(mesh_n[[i]])
  # convert domain polygon into a Spatial Polygons
  usabdy_sp <- SpatialPolygons(list(Polygons(list(Polygon(usa_utm)), ID = "1")))
  library(rgeos)
  usabdy_sp <- gBuffer(usabdy_sp, byid = TRUE, width = 0)
  # compute intersection between each polygon
  mesh_weight[[i]] <- sapply(1:length(dual_mesh[[i]]), 
                             function(x) {
                               if (gIntersects(dual_mesh[[i]][x, ], usabdy_sp)) 
                                 return(gArea(gIntersection(dual_mesh[[i]][x, ], usabdy_sp)))
                               else 
                                 return(0) 
                             })
  # set up the SPDE model
  ufo_spde[[i]] <- inla.spde2.pcmatern(mesh = mesh_n[[i]], 
                                       alpha = 2,
                                       prior.range = c(0.1, 0.01), # P(range < 0.1) = 0.01
                                       prior.sigma = c(0.01, 0.01))
}


### --------------------------------------------------------------------------->> Projection Matrices
y.pp <- temp
e.pp <- temp
imat <- temp
lmat <- temp
A.pp <- temp
stk.pp <- temp


for (j in 1:nrow(mesh_mat)){
  # define a vector of ones of the observations and zeros for the mesh nodes
  y.pp[[j]] <- rep(0:1, c(nv[j], n))
  # define the exposure vector 
  e.pp[[j]] <- c(mesh_weight[[j]], rep(0, n))
  # projection matrix
  imat[[j]] <- Diagonal(nv[j], rep(1, nv[j]))
  lmat[[j]]<- inla.spde.make.A(mesh = mesh_n[[j]], loc = coordinates(ufo_sp))
  A.pp[[j]] <- rbind(imat[[j]], lmat[[j]])
  # set up the data stack
  stk.pp[[j]] <- inla.stack(
    data = list(y = y.pp[[j]], e = e.pp[[j]]),
    A = list(1, A.pp[[j]]),
    effects = list(list(b0 = rep(1, nv[j] + n)), list(i = 1:nv[j])),
    tag = "pp")
}



### --------------------------------------------------------------------------->> Model without covariate
pp.res <- temp
pp.res.est <- temp

for (k in 1:nrow(mesh_mat)){
  pp.res[[k]] <- inla(y ~ 0 + b0 + f(i, model = ufo_spde[[k]]),
                      family = 'poisson', 
                      data = inla.stack.data(stk.pp[[k]]),
                      control.predictor = list(A = inla.stack.A(stk.pp[[k]])),
                      control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
                      E = inla.stack.data(stk.pp[[k]])$e)
  # transform
  pp.res.est[[k]] <- inla.spde.result(inla = pp.res[[k]], name = "i", spde = ufo_spde[[k]], do.transf = TRUE)
}



### --------------------------------------------------------------------------->> Adding a covariate

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




pop_mesh <- temp
pp.cov <- temp
pp.cov.est <- temp

for (l in 1: nrow(mesh_mat)){
  pop_mesh[[l]] <- sp::over(SpatialPoints(mesh_n[[l]]$loc[,1:2], 
                                          proj4string = CRS(proj4string(ufo_sp))), us_pop20)$TOT_POP
  pop_obs <- sp::over(ufo_sp, us_pop20)$TOT_POP
  covs[[l]] <- data.frame(pop = c(pop_mesh[[l]], pop_obs))
  # data stack includes the covariate
  stk.cov <- inla.stack(
    data = list(y = y.pp[[l]], e = e.pp[[l]]), 
    A = list(1, A.pp[[l]]),
    effects = list(list(b0 = rep(1, nv[l] + n), pop = covs[[l]]$pop), 
                   list(i = 1:nv[l])),
    tag = 'pp_cov')
  # fit model with a covariate
  pp.cov[[l]] <- inla(y ~ 0 + b0 + pop + f(i, model = ufo_spde), 
                      family = 'poisson', data = inla.stack.data(stk.cov), 
                      control.predictor = list(A = inla.stack.A(stk.cov)), 
                      E = inla.stack.data(stk.cov)$e)
  pp.cov.est[[l]] <- inla.spde.result(inla = pp.cov, name = "i", spde = ufo_spde, do.transf = TRUE)
  
}




### --------------------------------------------------------------------------->> Final output CSV file

list_fixed <- tmp
marg_fixed <- tmp
marg_kappa <- tmp
marg_variance <- tmp
marg_range <- tmp
list_dic <- tmp
list_waic <- tmp

for (f in 1:nrow(mesh_mat)){
  list_fixed[[f]] <- pp.res[[f]]$summary.fixed
  marg_fixed[[f]] <- pp.res[[f]]$marginals.fixed[[1]]
  marg_kappa[[f]] <- pp.res.est[[f]]$marginals.kappa[[1]]
  marg_variance[[f]] <- pp.res.est[[f]]$marginals.variance.nominal[[1]]
  marg_range[[f]] <- pp.res.est[[f]]$marginals.range.nominal[[1]]
  list_dic[[f]] <- pp.res[[f]]$dic$dic
  list_waic[[f]] <- pp.res[[f]]$waic$waic
}








