
## packages
library(sf)
library(qpdf)
library(maps)
library(sp)
library(INLA)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork)
library(maptools)
library(rgeos)
library(lattice)
# INLA:::inla.binary.install()

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
coordinates(ufo_locs) <- c("city_longitude", "city_latitude")
# Setting default projection
proj4string(ufo_locs) <- CRS('+init=epsg:4326')
## using North America Lambert Conformal Conic
ufo_sp <- spTransform(ufo_locs,
                      CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs "))
# US region / domain
usa_utm <- spTransform(region, 
                       CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))


### --------------------------------------------------------------------------->> Mesh Construction 
### We compare six results for the ufo dataset based on the six different meshes
### mesh_val <- c(loc, boundary, max.edge, cutoff)
source("ufo_mesh.R")


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
  dual_mesh[[i]] <- book.mesh.dual(mesh_n[[i]])
  # convert domain polygon into a Spatial Polygons
  usabdy_sp <- SpatialPolygons(list(Polygons(list(Polygon(usa_utm)), ID = "1")))
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
times.INLA <- temp
times.stelfi <- temp
stelfi.fit <- temp
valid_meshes_INLA <- c()
valid_meshes_stelfi <- c()
locs_stelfi <- coordinates(ufo_sp)
locs_stelfi <- data.frame(x = locs_stelfi[,1], y = locs_stelfi[,2])
names(locs_stelfi) <- c("x", "y")

for (k in 1:nrow(mesh_mat)){
  # Fit using INLA
  start <- as.numeric(Sys.time())
  out <- tryCatch(
    {
    pp.res[[k]] <- inla(y ~ 0 + b0 + f(i, model = ufo_spde[[k]]),
                      family = 'poisson', 
                      data = inla.stack.data(stk.pp[[k]]),
                      control.predictor = list(A = inla.stack.A(stk.pp[[k]])),
                      control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
                      E = inla.stack.data(stk.pp[[k]])$e)
    # transform
    pp.res.est[[k]] <- inla.spde.result(inla = pp.res[[k]], name = "i", spde = ufo_spde[[k]], do.transf = TRUE)
    out <- k
    },
    error = function(cond) {
      return(0)
    }
  )
  times.INLA[[k]] <- as.numeric(Sys.time()) - start
  if (out == k) {
    valid_meshes_INLA = append(valid_meshes_INLA, k)
  }
  # Fit using stelfi
  start <- as.numeric(Sys.time())
  out <- tryCatch(
    {
    res <- stelfi::fit_lgcp(locs = locs_stelfi, sp = usa_utm,
                          smesh = mesh_n[[k]])
    stelfi.fit[[k]] <- TMB::sdreport(res)
    out <- k
    },
    error = function(cond) {
      return(0)
    }
  )
  times.stelfi <- as.numeric(Sys.time()) - start
  if (out == k) {
    valid_meshes_stelfi = append(valid_meshes_stelfi, k)
  }
  
}

tempn <- vector("list", length = length(valid_meshes_INLA))

## Projection on a grid

gproj <- tempn
g_mean <- tempn
g_sd <- tempn

for (k in 1:length(valid_meshes_INLA)) {
  l = valid_meshes_INLA[k]
  gproj[[k]] <- inla.mesh.projector(mesh_n[[l]], dims = c(1000, 1000))
  g_mean[[k]] <- inla.mesh.project(gproj[[k]], pp.res[[l]]$summary.random$i$mean)
  g_sd[[k]] <- inla.mesh.project(gproj[[k]], pp.res[[l]]$summary.random$i$sd)
}

plot_1 <- do.call('grid.arrange',
                  lapply(g_mean,
                         levelplot, col.regions=terrain.colors(16), scales=list(draw=FALSE), 
                         main='latent field mean', xlab='', ylab=''))
ggsave(path = "ufo_sim_files", filename = "latent_field_mean_INLA.png", plot = plot_1,
       units = "in", width = 10, height = 10)


plot_2 <- do.call('grid.arrange',
                  lapply(g_sd,
                         levelplot, col.regions=terrain.colors(16), scales=list(draw=FALSE), 
                         main='latent field SD', xlab='', ylab=''))
ggsave(path = "ufo_sim_files", filename = "latent_field_sd_INLA.png", plot = plot_2,
       units = "in", width = 10, height = 10)

for(l in valid_meshes_stelfi) {
  plot_3 = stelfi::show_field(x = stelfi.fit[[l]]$par.random, smesh = mesh_n[[l]],
                              border = usa_utm)
  filename = paste("latent_field_stelfi_",l, ".png",sep="")
  ggsave(path = "ufo_sim_files", filename = filename, plot = plot_3,
         units = "in", width = 10, height = 10)
}

dev.off()

### --------------------------------------------------------------------------->> Final output CSV file

dir.create("CSV_files")

list_fixed <- tempn
marg_fixed <- tempn
marg_kappa <- tempn
marg_variance <- tempn
marg_range <- tempn
list_dic <- tempn
list_waic <- tempn


for (k in 1:length(valid_meshes_INLA)){
  f = valid_meshes_INLA[k]
  list_fixed[[k]] <- pp.res[[f]]$summary.fixed
  marg_fixed[[k]] <- pp.res[[f]]$marginals.fixed[[1]]
  marg_kappa[[k]] <- pp.res.est[[f]]$marginals.kappa[[1]]
  marg_variance[[k]] <- pp.res.est[[f]]$marginals.variance.nominal[[1]]
  marg_range[[k]] <- pp.res.est[[f]]$marginals.range.nominal[[1]]
  list_dic[[k]] <- pp.res[[f]]$dic$dic
  list_waic[[k]] <- pp.res[[f]]$waic$waic
}

### --------------------------------------------------------------------------->> CSV summary.fixed
fixed <- data.frame(matrix(unlist(list_fixed), ncol = 7, byrow = TRUE))
colnames(fixed) <- names(pp.res[[1]]$summary.fixed)
fixed$Mesh <- paste("Mesh", valid_meshes_INLA, sep = "_")
write.csv(fixed, file = "CSV_files/summary_fixed.csv")


### --------------------------------------------------------------------------->> CSV marginals.fixed
margfixed <- data.frame(do.call('rbind', marg_fixed))
margfixed$Mesh <- rep(paste("Mesh", valid_meshes_INLA, sep = "_"), each = nrow(pp.res[[1]]$marginals.fixed[[1]]))
write.csv(margfixed, file = "CSV_files/marginals_fixed.csv")


### --------------------------------------------------------------------------->> CSV marginals.log.kappa
margkappa <- data.frame(do.call('rbind', marg_kappa))
margkappa$Mesh <- rep(paste("Mesh", valid_meshes_INLA, sep = "_"), each = nrow(pp.res.est[[1]]$marginals.kappa[[1]]))
write.csv(margkappa, file = "CSV_files/marginals_kappa.csv")


### --------------------------------------------------------------------------->> CSV marginals.log.variance.nominal
margvar <- data.frame(do.call('rbind', marg_variance))
margvar$Mesh <- rep(paste("Mesh", valid_meshes_INLA, sep = "_"), each = nrow(pp.res.est[[1]]$marginals.variance.nominal[[1]]))
write.csv(margvar, file = "CSV_files/marginals_variance.csv")


### --------------------------------------------------------------------------->> CSV marginals.log.range.nominal
margrange <- data.frame(do.call('rbind', marg_range))
margrange$Mesh <- rep(paste("Mesh", valid_meshes_INLA, sep = "_"), each = nrow(pp.res.est[[1]]$marginals.range.nominal[[1]]))
write.csv(margrange, file = "CSV_files/marginals_range.csv")


### --------------------------------------------------------------------------->> CSV DIC
dic <- data.frame(matrix(unlist(list_dic), ncol = 1, byrow = TRUE))
colnames(dic) <- "DIC"
dic$Mesh <- paste("Mesh", valid_meshes_INLA, sep = "_")
write.csv(dic, file = "CSV_files/DIC.csv")


### --------------------------------------------------------------------------->> CSV WAIC
waic <- data.frame(matrix(unlist(list_waic), ncol = 1, byrow = TRUE))
colnames(waic) <- "WAIC"
waic$Mesh <- paste("Mesh", valid_meshes_INLA, sep = "_")
write.csv(waic, file = "CSV_files/WAIC.csv")

### --------------------------------------------------------------------------->> CSV INLA Computation time
times <- data.frame(matrix(unlist(times.INLA), ncol = 1, byrow = TRUE))
colnames(times) <- "CompTime"
times$Mesh <- paste("Mesh", 1:nrow(mesh_mat), sep = "_")
write.csv(times, file = "CSV_files/Times_INLA.csv")

### --------------------------------------------------------------------------->> CSV stelfi Computation time
times <- data.frame(matrix(unlist(times.stelfi), ncol = 1, byrow = TRUE))
colnames(times) <- "CompTime"
times$Mesh <- paste("Mesh", 1:nrow(mesh_mat), sep = "_")
write.csv(times, file = "CSV_files/Times_stelfi.csv")


### --------------------------------------------------------------------------->> Adding a covariate

### US 2020 Cartographic Boundary Files
### downloaded from https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_us_county_20m.zip
US20 <- st_read("data/cb_2020_us_county_20m/cb_2020_us_county_20m.shp")
# leave out AK, HI, and PR (state FIPS: 02, 15, and 72)
conti_US20 <- US20[!(US20$STATEFP %in% c("02","15","72")), ]
### 2020 Census Demographic Data By County, downloaded from https://www.dataemporium.com/dataset/254/
census20 <- read.csv("data/US_2020_census.csv")
colnames(census20)[1] = "STATE_ABBR"
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
us_pop20 <- spTransform(us_pop20, 
                        CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs "))


### model with population as covariate
pop_mesh <- temp
covs <- temp
stk.cov <- temp
pp.cov <- temp
pp.cov.est <- temp
pp.cov.times <- temp
valid_meshes_INLA <- c()

for (l in 1:nrow(mesh_mat)){
  pop_mesh[[l]] <- sp::over(SpatialPoints(mesh_n[[l]]$loc[,1:2], 
                                          proj4string = CRS(proj4string(ufo_sp))), us_pop20)$TOT_POP
  pop_obs <- sp::over(ufo_sp, us_pop20)$TOT_POP
  covs[[l]] <- data.frame(pop = c(pop_mesh[[l]], pop_obs))
  # data stack includes the covariate
  stk.cov[[l]] <- inla.stack(data = list(y = y.pp[[l]], e = e.pp[[l]]), 
                             A = list(1, A.pp[[l]]),
                             effects = list(list(b0 = rep(1, nv[l] + n), pop = covs[[l]]$pop), 
                                            list(i = 1:nv[l])),
                             tag = 'pp_cov')
  # fit model with a covariate
  start <- Sys.time()
  out <- tryCatch(
    {
    pp.cov[[l]] <- inla(y ~ 0 + b0 + pop + f(i, model = ufo_spde[[l]]), 
                      family = 'poisson', 
                      data = inla.stack.data(stk.cov[[l]]), 
                      control.predictor = list(A = inla.stack.A(stk.cov[[l]])),
                      control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
                      E = inla.stack.data(stk.cov[[l]])$e)
    pp.cov.est[[l]] <- inla.spde.result(inla = pp.cov[[l]], name = "i", spde = ufo_spde[[l]], do.transf = TRUE)
    out <- l
    },
    error = function(cond) {
      return(0)
    }
  )
  pp.cov.times[[l]] <- Sys.time() - start
  if (out == l) {
    valid_meshes_INLA <- append(valid_meshes_INLA,l)
  }
}

tempn <- vector("list", length = length(valid_meshes_INLA))



gproj_cov <- tempn
g_mean_cov <- tempn
g_sd_cov <- tempn

for (k in 1:length(valid_meshes_INLA)) {
  g = valid_meshes_INLA[k]
  gproj_cov[[k]] <- inla.mesh.projector(mesh_n[[g]], dims = c(1000, 1000))
  g_mean_cov[[k]] <- inla.mesh.project(gproj_cov[[k]], pp.cov[[g]]$summary.random$i$mean)
  g_sd_cov[[k]] <- inla.mesh.project(gproj_cov[[k]], pp.cov[[g]]$summary.random$i$sd)
}

plot_m_cov <- do.call('grid.arrange',
                      lapply(g_mean_cov,
                             levelplot, col.regions=terrain.colors(16), scales=list(draw=FALSE), 
                             main='latent field mean', xlab='', ylab=''))
ggsave(path = "ufo_sim_files", filename = "latent_field_mean_cov.png", plot_m_cov,
       units = "in", width = 10, height = 10)


plot_sd_cov <- do.call('grid.arrange',
                       lapply(g_sd_cov,
                              levelplot, col.regions=terrain.colors(16), scales=list(draw=FALSE), 
                              main='latent field SD', xlab='', ylab=''))
ggsave(path = "ufo_sim_files", filename = "latent_field_sd_cov.png", plot_sd_cov,
       units = "in", width = 10, height = 10)



### --------------------------------------------------------------------------->> Final output w/ covariate CSV file

dir.create("cov_CSV_files")

list_fixed <- tempn
marg_fixed <- tempn
marg_kappa <- tempn
marg_variance <- tempn
marg_range <- tempn
list_dic <- tempn
list_waic <- tempn


for (k in 1:length(valid_meshes_INLA)){
  f = valid_meshes_INLA[k]
  list_fixed[[k]] <- pp.cov[[f]]$summary.fixed
  marg_fixed[[k]] <- pp.cov[[f]]$marginals.fixed[[1]]
  marg_kappa[[k]] <- pp.cov.est[[f]]$marginals.kappa[[1]]
  marg_variance[[k]] <- pp.cov.est[[f]]$marginals.variance.nominal[[1]]
  marg_range[[k]] <- pp.cov.est[[f]]$marginals.range.nominal[[1]]
  list_dic[[k]] <- pp.cov[[f]]$dic$dic
  list_waic[[k]] <- pp.cov[[f]]$waic$waic
}

### --------------------------------------------------------------------------->> CSV summary.fixed
fixed <- data.frame(matrix(unlist(list_fixed), ncol = 7, byrow = TRUE))
colnames(fixed) <- names(pp.res[[1]]$summary.fixed)
fixed$Mesh <- paste("Mesh", valid_meshes_INLA, sep = "_")
write.csv(fixed, file = "cov_CSV_files/summary_fixed_cov.csv")


### --------------------------------------------------------------------------->> CSV marginals.fixed
margfixed <- data.frame(do.call('rbind', marg_fixed))
margfixed$Mesh <- rep(paste("Mesh", valid_meshes_INLA, sep = "_"), each = nrow(pp.res[[1]]$marginals.fixed[[1]]))
write.csv(margfixed, file = "cov_CSV_files/marginals_fixed_cov.csv")


### --------------------------------------------------------------------------->> CSV marginals.log.kappa
margkappa <- data.frame(do.call('rbind', marg_kappa))
margkappa$Mesh <- rep(paste("Mesh", valid_meshes_INLA, sep = "_"), each = nrow(pp.res.est[[1]]$marginals.kappa[[1]]))
write.csv(margkappa, file = "cov_CSV_files/marginals_kappa_cov.csv")


### --------------------------------------------------------------------------->> CSV marginals.log.variance.nominal
margvar <- data.frame(do.call('rbind', marg_variance))
margvar$Mesh <- rep(paste("Mesh", valid_meshes_INLA, sep = "_"), each = nrow(pp.res.est[[1]]$marginals.variance.nominal[[1]]))
write.csv(margvar, file = "cov_CSV_files/marginals_variance_cov.csv")


### --------------------------------------------------------------------------->> CSV marginals.log.range.nominal
margrange <- data.frame(do.call('rbind', marg_range))
margrange$Mesh <- rep(paste("Mesh", valid_meshes_INLA, sep = "_"), each = nrow(pp.res.est[[1]]$marginals.range.nominal[[1]]))
write.csv(margrange, file = "cov_CSV_files/marginals_range_cov.csv")


### --------------------------------------------------------------------------->> CSV DIC
dic <- data.frame(matrix(unlist(list_dic), ncol = 1, byrow = TRUE))
colnames(dic) <- "DIC"
dic$Mesh <- paste("Mesh", valid_meshes_INLA, sep = "_")
write.csv(dic, file = "cov_CSV_files/DIC_cov.csv")


### --------------------------------------------------------------------------->> CSV WAIC
waic <- data.frame(matrix(unlist(list_waic), ncol = 1, byrow = TRUE))
colnames(waic) <- "WAIC"
waic$Mesh <- paste("Mesh", valid_meshes_INLA, sep = "_")
write.csv(waic, file = "cov_CSV_files/WAIC_cov.csv")

### --------------------------------------------------------------------------->> CSV computation time
times <- data.frame(matrix(unlist(pp.cov.times), ncol = 1, byrow = TRUE))
colnames(times) <- "CompTime"
times$Mesh <- paste("Mesh", 1:nrow(mesh_mat), sep = "_")
write.csv(times, file = "cov_CSV_files/Times.csv")






