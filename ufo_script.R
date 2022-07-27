
## packages
devtools::load_all("stelfi")
library(sf)
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
library(inlabru)
# INLA:::inla.binary.install()
#setwd("GitHub/meshmetrics")

dir.create("UFO_stelfi_maps")
dir.create("UFO_bru_maps")
dir.create("UFO_meshes")


## read in data
ufo <- readr::read_csv("data/ufo.csv")
## usa map
usa <- map_data('usa')
#usa_region <- data.frame(Longitude = usa$long, Latitude = usa$lat)

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
                      #CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs "))
                      CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs "))
# US region / domain
usa_utm <- spTransform(region, 
                       #CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))
                       CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"))
bound <- inla.sp2segment(usa_utm)
locs <- coordinates(ufo_sp)
locs <- data.frame(x = locs[,1], y = locs[,2])

### --------------------------------------------------------------------------->> Mesh Construction 

tmp <- lapply(seq(150, 1050, length.out = 10), function(x) INLA::inla.mesh.2d(loc.domain = bound$loc,
                                                                                   max.edge = c(x, 2 * x)))
# create random mesh nodes
win <- as.owin.SpatialPolygons(usa_utm)
## choose resolution
set.seed(4321)
points <- list()
for (i in 1:10) {
  X <- rLGCP("matern", mu = -10,
             var = 0.5, scale = 50, nu = 1,
             win = win)
  Lamda <- attr(X, 'Lambda')
  points[[i]] <- data.frame(x = X$x, y = X$y)
}

tmp2 = lapply(seq(150, 1050, length.out = 10), function(x) INLA::inla.mesh.2d(loc.domain = bound$loc,
                                                                              loc = points[[(as.integer(0.01*x - 0.5))]],
                                                                              cutoff = 0.075*x,
                                                                              max.edge = c(x, 2 * x)))

tmp = append(tmp, tmp2)

attrs <- lapply(tmp, stelfi:::meshmetrics)
results <- as.data.frame(matrix(0, nrow = length(tmp), ncol = 8))
names(results) <- c(paste("bru_", 1:3, sep = ""),"bru_time", paste("stelfi_", 1:3, sep = ""), "stelfi_time")

## append mesh attributes summaries
results$mean_radius_edge <- sapply(attrs, function(x) mean(x$triangles$radius_edge))
results$sd_radius_edge <- sapply(attrs, function(x) sd(x$triangles$radius_edge))
results$uq_radius_edge <- sapply(attrs, function(x) quantile(x$triangles$radius_edge)[4])
results$mean_radius_ratio <- sapply(attrs, function(x) mean(x$triangles$radius_ratio))
results$sd_radius_ratio <- sapply(attrs, function(x) sd(x$triangles$radius_ratio))
results$uq_radius_ratio <- sapply(attrs, function(x) quantile(x$triangles$radius_ratio)[4])
results$n_triangles <- sapply(tmp, function(x) x$n)


points <- locs
coordinates(locs) <- c("x", "y")
for (i in 1:length(tmp)) {
  # plot mesh
  filepath <- paste("UFO_meshes/UFO_mesh_",i, ".jpg", sep="")
  jpeg(filepath, width = 1440, height = 1440)
  plot(tmp[[i]])
  lines(coordinates(usa_utm@polygons[[1]]@Polygons[[1]])[,1], coordinates(usa_utm@polygons[[1]]@Polygons[[1]])[,2], col = "red")
  dev.off()
  
  # plot mesh attributes
  filename <- paste("UFO_mesh_attrs_",i, ".jpg", sep="")
  mesh_plot <- stelfi::plot_mesh(tmp[[i]])
  ggsave(path = "UFO_meshes", filename = filename, plot = mesh_plot,
         units = "in", width = 10, height = 10)
  
  
  ## Define SPDE prior
  matern <- INLA::inla.spde2.pcmatern(tmp[[i]],
                                      prior.sigma = c(0.1, 0.01),
                                      prior.range = c(100, 0.01)
  )
  
  ## Define domain of the LGCP as well as the model components (spatial SPDE
  ## effect and Intercept)
  cmp <- coordinates ~ mySmooth(main = coordinates, model = matern) + Intercept(1)
  ## Fit the model (with int.strategy = "eb" to make the example take less time)
  out <- tryCatch(
    {
    t_bru <- system.time(fit_bru <- lgcp(cmp, locs,
                                        samplers = usa_utm,
                                        domain = list(coordinates = tmp[[i]]),
                                        options = list(control.inla = list(int.strategy = "eb")))
    )
    results[i, 1:4] <- c(fit_bru$summary.fixed[, 1], fit_bru$summary.hyperpar[, 1], as.numeric(t_bru[3]))
    field_plot <- stelfi::show_field(x=fit_bru$summary.random$mySmooth$mean, smesh = tmp[[i]], border = usa_utm)
    filename <- paste("UFO_field_", i, ".png", sep = "")
    ggsave(path = "UFO_bru_maps", filename = filename, plot = field_plot,
           units = "in", width = 10, height = 10)
    },
    error = function(cond) {
      results[i, 1:4] <- rep(NA,4)
    }
  )
  ## Stelfi fit
  out <- tryCatch(
    {
    t_stelfi <- system.time(fit_stelfi <- stelfi::fit_lgcp(locs = points, sp = usa_utm, smesh = tmp[[i]]))
    results[i, 5:8] <- c(get_coefs(fit_stelfi)[c(1, 4:5), 1], as.numeric(t_stelfi[3]))
    res <- TMB::sdreport(fit_stelfi)
    field_plot <- stelfi::show_field(x = res$par.random, smesh = tmp[[i]], border = usa_utm)
    filename <- paste("UFO_field_", i, ".png", sep = "")
    ggsave(path = "UFO_stelfi_maps", filename = filename, plot = field_plot,
           units = "in", width = 10, height = 10)
    },
    error = function(cond) {
      results[i, 5:8] <- rep(NA,4)
    }
  )
}

write.csv(results, file = paste("res_ufo", ".csv", sep = ""))

# ### --------------------------------------------------------------------------->> Dual Mesh and Weights
# temp <- vector("list", length = nrow(mesh_mat))
# nv <- c() # number of vertices in the mesh
# dual_mesh <- temp
# mesh_weight <- temp
# ufo_spde <- temp
# n <- nrow(ufo) # number of observations
# 
# for (i in 1:nrow(mesh_mat)) {
#   # number of vertices in the mesh
#   nv[i] <- mesh_n[[i]]$n
#   # create the dual mesh polygons
#   dual_mesh[[i]] <- book.mesh.dual(mesh_n[[i]])
#   # convert domain polygon into a Spatial Polygons
#   usabdy_sp <- SpatialPolygons(list(Polygons(list(Polygon(usa_utm)), ID = "1")))
#   usabdy_sp <- gBuffer(usabdy_sp, byid = TRUE, width = 0)
#   # compute intersection between each polygon
#   mesh_weight[[i]] <- stelfi::get_weights(mesh_n[[i]], usa_utm)$weights
#   # set up the SPDE model
#   ufo_spde[[i]] <- inla.spde2.pcmatern(mesh = mesh_n[[i]], 
#                                        alpha = 2,
#                                        prior.range = c(0.1, 0.01), # P(range < 0.1) = 0.01
#                                        prior.sigma = c(0.01, 0.01))
# }
# 
# 
# ### --------------------------------------------------------------------------->> Projection Matrices
# y.pp <- temp
# e.pp <- temp
# imat <- temp
# lmat <- temp
# A.pp <- temp
# stk.pp <- temp
# 
# 
# for (j in 1:nrow(mesh_mat)){
#   # define a vector of ones of the observations and zeros for the mesh nodes
#   y.pp[[j]] <- rep(0:1, c(nv[j], n))
#   # define the exposure vector 
#   e.pp[[j]] <- c(mesh_weight[[j]], rep(0, n))
#   # projection matrix
#   imat[[j]] <- Diagonal(nv[j], rep(1, nv[j]))
#   lmat[[j]]<- inla.spde.make.A(mesh = mesh_n[[j]], loc = coordinates(ufo_sp))
#   A.pp[[j]] <- rbind(imat[[j]], lmat[[j]])
#   # set up the data stack
#   stk.pp[[j]] <- inla.stack(
#     data = list(y = y.pp[[j]], e = e.pp[[j]]),
#     A = list(1, A.pp[[j]]),
#     effects = list(list(b0 = rep(1, nv[j] + n)), list(i = 1:nv[j])),
#     tag = "pp")
# }
# 
# 
# 
# ### --------------------------------------------------------------------------->> Model without covariate
# pp.res <- temp
# pp.res.est <- temp
# times.INLA <- temp
# times.stelfi <- temp
# stelfi.fit <- temp
# valid_meshes_INLA <- c()
# valid_meshes_stelfi <- c()
# locs_stelfi <- coordinates(ufo_sp)
# locs_stelfi <- data.frame(x = locs_stelfi[,1], y = locs_stelfi[,2])
# names(locs_stelfi) <- c("x", "y")
# 
# for (k in 1:nrow(mesh_mat)){
#   # Fit using INLA
#   start <- as.numeric(Sys.time())
#   out <- tryCatch(
#     {
#     pp.res[[k]] <- inla(y ~ 0 + b0 + f(i, model = ufo_spde[[k]]),
#                       family = 'poisson', 
#                       data = inla.stack.data(stk.pp[[k]]),
#                       control.predictor = list(A = inla.stack.A(stk.pp[[k]])),
#                       control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
#                       E = inla.stack.data(stk.pp[[k]])$e)
#     # transform
#     pp.res.est[[k]] <- inla.spde.result(inla = pp.res[[k]], name = "i", spde = ufo_spde[[k]], do.transf = TRUE)
#     out <- k
#     },
#     error = function(cond) {
#       return(0)
#     }
#   )
#   times.INLA[[k]] <- as.numeric(Sys.time()) - start
#   if (out == k) {
#     valid_meshes_INLA = append(valid_meshes_INLA, k)
#   }
#   # Fit using stelfi
#   start <- as.numeric(Sys.time())
#   out <- tryCatch(
#     {
#     res <- stelfi::fit_lgcp(locs = locs_stelfi, sp = usa_utm,
#                           smesh = mesh_n[[k]])
#     stelfi.fit[[k]] <- TMB::sdreport(res)
#     out <- k
#     },
#     error = function(cond) {
#       return(0)
#     }
#   )
#   times.stelfi[[k]] <- as.numeric(Sys.time()) - start
#   if (out == k) {
#     valid_meshes_stelfi = append(valid_meshes_stelfi, k)
#   }
#   
# }
# 
# tempn <- vector("list", length = length(valid_meshes_INLA))
# 
# ## Projection on a grid
# 
# gproj <- tempn
# g_mean <- tempn
# g_sd <- tempn
# 
# for (k in 1:length(valid_meshes_INLA)) {
#   l = valid_meshes_INLA[k]
#   gproj[[k]] <- inla.mesh.projector(mesh_n[[l]], dims = c(1000, 1000))
#   g_mean[[k]] <- inla.mesh.project(gproj[[k]], pp.res[[l]]$summary.random$i$mean)
#   g_sd[[k]] <- inla.mesh.project(gproj[[k]], pp.res[[l]]$summary.random$i$sd)
# }
# 
# plot_1 <- do.call('grid.arrange',
#                   lapply(g_mean,
#                          levelplot, col.regions=terrain.colors(16), scales=list(draw=FALSE), 
#                          main='latent field mean', xlab='', ylab=''))
# ggsave(path = "ufo_sim_files", filename = "latent_field_mean_INLA.png", plot = plot_1,
#        units = "in", width = 10, height = 10)
# 
# 
# plot_2 <- do.call('grid.arrange',
#                   lapply(g_sd,
#                          levelplot, col.regions=terrain.colors(16), scales=list(draw=FALSE), 
#                          main='latent field SD', xlab='', ylab=''))
# ggsave(path = "ufo_sim_files", filename = "latent_field_sd_INLA.png", plot = plot_2,
#        units = "in", width = 10, height = 10)
# 
# for(l in valid_meshes_stelfi) {
#   plot_3 = stelfi::show_field(x = stelfi.fit[[l]]$par.random, smesh = mesh_n[[l]],
#                               border = usa_utm)
#   filename = paste("latent_field_stelfi_",l, ".png",sep="")
#   ggsave(path = "ufo_sim_files", filename = filename, plot = plot_3,
#          units = "in", width = 10, height = 10)
# }
# 
# dev.off()
# 
# ### --------------------------------------------------------------------------->> Final output CSV file
# 
# dir.create("UFO_CSV_files")
# 
# # Record fitted values by INLA
# 
# list_fixed_INLA <- tempn
# list_hyperpar_INLA <- tempn
# #marg_fixed <- tempn
# #marg_kappa <- tempn
# #marg_variance <- tempn
# #marg_range <- tempn
# list_dic <- tempn
# list_waic <- tempn
# 
# 
# for (k in 1:length(valid_meshes_INLA)){
#   f = valid_meshes_INLA[k]
#   list_fixed_INLA[[k]] <- pp.res[[f]]$summary.fixed
#   list_hyperpar_INLA[[k]] <- pp.res[[f]]$summary.hyperpar
#   #marg_fixed[[k]] <- pp.res[[f]]$marginals.fixed[[1]]
#   #marg_kappa[[k]] <- pp.res.est[[f]]$marginals.kappa[[1]]
#   #marg_variance[[k]] <- pp.res.est[[f]]$marginals.variance.nominal[[1]]
#   #marg_range[[k]] <- pp.res.est[[f]]$marginals.range.nominal[[1]]
#   list_dic[[k]] <- pp.res[[f]]$dic$dic
#   list_waic[[k]] <- pp.res[[f]]$waic$waic
# }
# 
# # Record fitted values by stelfi
# 
# tempn <- vector("list", length = length(valid_meshes_stelfi))
# 
# list_fixed_stelfi <- tempn
# list_sd_stelfi <- tempn
# for (k in 1:length(valid_meshes_stelfi)){
#   f = valid_meshes_INLA[k]
#   list_fixed_stelfi[[k]] <- stelfi.fit[[f]]$value
#   list_sd_stelfi[[k]] <- stelfi.fit[[f]]$sd
# }
# 
# ### --------------------------------------------------------------------------->> CSV summary.fixed (INLA)
# fixed <- data.frame(matrix(unlist(list_fixed_INLA), ncol = 7, byrow = TRUE))
# colnames(fixed) <- names(pp.res[[1]]$summary.fixed)
# fixed$Mesh <- paste("Mesh", valid_meshes_INLA, sep = "_")
# write.csv(fixed, file = "UFO_CSV_files/summary_fixed_INLA.csv")
# 
# ### --------------------------------------------------------------------------->> CSV hyperparameters (INLA)
# fixed <- data.frame(matrix(unlist(list_hyperpar_INLA), ncol = 12, byrow = TRUE))
# colnames(fixed) <- rep(names(pp.res[[1]]$summary.hyperpar),2)
# fixed$Mesh <- paste("Mesh", valid_meshes_INLA, sep = "_")
# write.csv(fixed, file = "UFO_CSV_files/summary_hyperpar_INLA.csv")
# 
# ### --------------------------------------------------------------------------->> CSV fitted values (stelfi)
# fixed <- data.frame(matrix(unlist(list_fixed_stelfi), ncol = 5, byrow = TRUE))
# colnames(fixed) <- names(stelfi.fit[[1]]$value)
# fixed$Mesh <- paste("Mesh", valid_meshes_stelfi, sep = "_")
# write.csv(fixed, file = "UFO_CSV_files/summary_fixed_stelfi.csv")
# 
# ### --------------------------------------------------------------------------->> CSV value sd (stelfi)
# fixed <- data.frame(matrix(unlist(list_sd_stelfi), ncol = 5, byrow = TRUE))
# colnames(fixed) <- names(stelfi.fit[[1]]$sd)
# fixed$Mesh <- paste("Mesh", valid_meshes_stelfi, sep = "_")
# write.csv(fixed, file = "UFO_CSV_files/summary_sd_stelfi.csv")
# 
# 
# ### --------------------------------------------------------------------------->> CSV marginals.fixed
# #margfixed <- data.frame(do.call('rbind', marg_fixed))
# #margfixed$Mesh <- rep(paste("Mesh", valid_meshes_INLA, sep = "_"), each = nrow(pp.res[[1]]$marginals.fixed[[1]]))
# #write.csv(margfixed, file = "UFO_CSV_files/marginals_fixed.csv")
# 
# 
# ### --------------------------------------------------------------------------->> CSV marginals.log.kappa
# #margkappa <- data.frame(do.call('rbind', marg_kappa))
# #margkappa$Mesh <- rep(paste("Mesh", valid_meshes_INLA, sep = "_"), each = nrow(pp.res.est[[1]]$marginals.kappa[[1]]))
# #write.csv(margkappa, file = "UFO_CSV_files/marginals_kappa.csv")
# 
# 
# ### --------------------------------------------------------------------------->> CSV marginals.log.variance.nominal
# #margvar <- data.frame(do.call('rbind', marg_variance))
# #margvar$Mesh <- rep(paste("Mesh", valid_meshes_INLA, sep = "_"), each = nrow(pp.res.est[[1]]$marginals.variance.nominal[[1]]))
# #write.csv(margvar, file = "UFO_CSV_files/marginals_variance.csv")
# 
# 
# ### --------------------------------------------------------------------------->> CSV marginals.log.range.nominal
# #margrange <- data.frame(do.call('rbind', marg_range))
# #margrange$Mesh <- rep(paste("Mesh", valid_meshes_INLA, sep = "_"), each = nrow(pp.res.est[[1]]$marginals.range.nominal[[1]]))
# #write.csv(margrange, file = "UFO_CSV_files/marginals_range.csv")
# 
# 
# ### --------------------------------------------------------------------------->> CSV DIC
# dic <- data.frame(matrix(unlist(list_dic), ncol = 1, byrow = TRUE))
# colnames(dic) <- "DIC"
# dic$Mesh <- paste("Mesh", valid_meshes_INLA, sep = "_")
# write.csv(dic, file = "UFO_CSV_files/DIC.csv")
# 
# 
# ### --------------------------------------------------------------------------->> CSV WAIC
# waic <- data.frame(matrix(unlist(list_waic), ncol = 1, byrow = TRUE))
# colnames(waic) <- "WAIC"
# waic$Mesh <- paste("Mesh", valid_meshes_INLA, sep = "_")
# write.csv(waic, file = "UFO_CSV_files/WAIC.csv")
# 
# ### --------------------------------------------------------------------------->> CSV INLA Computation time
# times <- data.frame(matrix(unlist(times.INLA), ncol = 1, byrow = TRUE))
# colnames(times) <- "CompTime"
# times$Mesh <- paste("Mesh", 1:nrow(mesh_mat), sep = "_")
# write.csv(times, file = "UFO_CSV_files/Times_INLA.csv")
# 
# ### --------------------------------------------------------------------------->> CSV stelfi Computation time
# times <- data.frame(matrix(unlist(times.stelfi), ncol = 1, byrow = TRUE))
# colnames(times) <- "CompTime"
# times$Mesh <- paste("Mesh", 1:nrow(mesh_mat), sep = "_")
# write.csv(times, file = "UFO_CSV_files/Times_stelfi.csv")


### --------------------------------------------------------------------------->> Adding a covariate






