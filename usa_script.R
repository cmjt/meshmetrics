# Simulate Point Process over a map of the USA
# Compare different meshes and parameter values

## read in "index" argument on NeSi
args <- commandArgs(trailingOnly = TRUE)
## index is the first element of this form slurm
## second for bash
k <- as.numeric(args[1])

devtools::load_all("stelfi")
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
library(inlabru)
library(spatstat)

## source required functions 
#source("functions.r")
## usa map and region
states <- map_data("state")
usa <- map_data('usa')
usa_region <- data.frame(Longitude = usa$long, Latitude = usa$lat)
region <- as(sf::st_as_sf(maps::map("usa", fill = TRUE, plot = FALSE)), "Spatial")

## convert the existing lat-long coordinates to UTM (easting-northing system)
usa_utm <- spTransform(region, 
                       CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))

sets <- expand.grid(data.frame(beta = c(-24,-21), ## per unit area
                               var = c(5, 0.1),
                               rho = c(50000, 10000)))


### Construct Simulations
#sp <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(usa_utm)), '0')))
bound <- inla.sp2segment(usa_utm)
proj4string(usa_utm) <- ""
win <- as.owin.SpatialPolygons(usa_utm)
## choose resolution
set.seed(4321 + k)
points <- list()
for (i in 1:nrow(sets)) {
  X <- rLGCP("matern", mu = sets[i, 1],
           var = sets[i, 2], scale = sets[i, 3], nu = 1,
           win = win)
  Lamda <- attr(X, 'Lambda')
  points[[i]] <- data.frame(x = X$x, y = X$y)
}
## mesh construction (9 for now)
tmp <- lapply(seq(200000, 2000000, length.out = 9), function(x) INLA::inla.mesh.2d(boundary = bound,
                                                                        max.edge = c(x, 2 * x)))
attrs <- lapply(tmp, stelfi:::meshmetrics)

##  model fitting
results <- as.data.frame(matrix(0, nrow = length(tmp)*length(points), ncol = 8))
names(results) <- c(paste("bru_", 1:3, sep = ""),"bru_time", paste("stelfi_", 1:3, sep = ""), "stelfi_time")
## append mesh attributes summaries
results$mean_radius_edge <- rep(sapply(attrs, function(x) mean(x$triangles$radius_edge)), each = length(points))
results$sd_radius_edge <- rep(sapply(attrs, function(x) sd(x$triangles$radius_edge)), each = length(points))
results$mean_radius_ratio <- rep(sapply(attrs, function(x) mean(x$triangles$radius_ratio)), each = length(points))
results$sd_radius_ratio <- rep(sapply(attrs, function(x) sd(x$triangles$radius_ratio)), each = length(points))

for (i in 1:length(tmp)) {
  for(j in 1:length(points)) {
    locs <- points[[j]]
    coordinates(locs) <- c("x", "y")
    ## Define SPDE prior
    matern <- INLA::inla.spde2.pcmatern(tmp[[i]],
                                        prior.sigma = c(0.1, 0.01),
                                        prior.range = c(0.01, 0.01)
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
        results[((i-1)*length(points) + j), 1:4] <- c(fit_bru$summary.fixed[, 1], fit_bru$summary.hyperpar[, 1], as.numeric(t_bru[3]))
      },
      error = function(cond) {
        results[((i-1)*length(points) + j), 1:4] <- rep(NA,4)
      }
    )
    ## Stelfi fit
    out <- tryCatch(
      {
        t_stelfi <- system.time(fit_stelfi <- stelfi::fit_lgcp(locs = points[[j]], sp = usa_utm, smesh = tmp[[i]]))
        results[((i-1)*length(points) + j), 5:8] <- c(get_coefs(fit_stelfi)[c(1, 4:5), 1], as.numeric(t_stelfi[3]))
      },
      error = function(cond) {
        results[((i-1)*length(points) + j), 5:8] <- rep(NA,4)
      }
    )
  }
}

write.csv(results, file = paste("res_", k, ".csv", sep = ""))


#### OLD INLA code

# rf <- log(Lamda$v)
# rf <- apply(rf, 2, t) ## field flips
# ex <- expand.grid(seq(win$xrange[1], win$xrange[2],length.out = npix),
#                   seq(win$yrange[1], win$yrange[2],length.out = npix))
# dat <- data.frame(x = ex$Var1, y = ex$Var2, rf = as.vector(t(rf)))
# ins <- inside.owin(dat$x, dat$y, win)
# dat$rf[!ins] <- NA
# # true <- ggplot(dat , aes(x = x, y = y, fill = rf)) +
# #   geom_tile() +
# #   labs(fill = "Expected (log) \nnumber \nof points") +
# #   geom_point(data = points, aes(x = x, y = y), inherit.aes = FALSE,
# #              shape = 18, size = 1, alpha = 0.5) +
# #   scale_fill_viridis_b(na.value = "transparent") +
# #   theme_void()  +
# #   geom_sf(data = sf::st_as_sf(sp), size = 1,
# #           inherit.aes = FALSE, fill = "transparent") +
# #   ggtitle("Simulated field and points")
# # 
# # #true
# 
# ## Fit LGCP, using both INLA and stelfi
# m1 <- INLA::inla.mesh.2d(loc.domain = bound$loc,
#                          max.edge = 1, max.n.strict = 1000)
# m2 <- INLA::inla.mesh.2d(loc.domain = bound$loc,
#                          max.edge = 6, max.n.strict = 200)
# m3 <- INLA::inla.mesh.2d(loc.domain = bound$loc,
#                          max.edge = 3, max.n.strict = 500)
# 
# mesh_n = list(m1,m2,m3)
# mesh_mat = matrix(0, nrow=length(mesh_n), ncol = 5)
# temp <- vector("list", length = nrow(mesh_mat))
# nv <- c() # number of vertices in the mesh
# dual_mesh <- temp
# mesh_weight <- temp
# usa_spde <- temp
# n <- nrow(points) # number of observations
# 
# for (i in 1:nrow(mesh_mat)) {
#   # number of vertices in the mesh
#   nv[i] <- mesh_n[[i]]$n
#   # create the dual mesh polygons
#   dual_mesh[[i]] <- book.mesh.dual(mesh_n[[i]])
#   # convert domain polygon into a Spatial Polygons
#   mesh_weight[[i]] <- stelfi::get_weights(mesh_n[[i]], usa_utm)$weights
#   # set up the SPDE model
#   usa_spde[[i]] <- inla.spde2.pcmatern(mesh = mesh_n[[i]], 
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
# locs = as.matrix(points)
# 
# 
# for (j in 1:nrow(mesh_mat)){
#   # define a vector of ones of the observations and zeros for the mesh nodes
#   y.pp[[j]] <- rep(0:1, c(nv[j], n))
#   # define the exposure vector 
#   e.pp[[j]] <- c(mesh_weight[[j]], rep(0, n))
#   # projection matrix
#   imat[[j]] <- Diagonal(nv[j], rep(1, nv[j]))
#   lmat[[j]]<- inla.spde.make.A(mesh = mesh_n[[j]], loc = locs)
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
# 
# for (k in 1:nrow(mesh_mat)){
#   # Fit using INLA
#   start <- as.numeric(Sys.time())
#   out <- tryCatch(
#     {
#       pp.res[[k]] <- inla(y ~ 0 + b0 + f(i, model = usa_spde[[k]]),
#                           family = 'poisson', 
#                           data = inla.stack.data(stk.pp[[k]]),
#                           control.predictor = list(A = inla.stack.A(stk.pp[[k]])),
#                           control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
#                           E = inla.stack.data(stk.pp[[k]])$e)
#       # transform
#       pp.res.est[[k]] <- inla.spde.result(inla = pp.res[[k]], name = "i", spde = usa_spde[[k]], do.transf = TRUE)
#       out <- k
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
#       res <- stelfi::fit_lgcp(locs = points, sp = usa_utm,
#                               smesh = mesh_n[[k]])
#       stelfi.fit[[k]] <- TMB::sdreport(res)
#       out <- k
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
# dir.create("USA_files")
# plot_1 <- do.call('grid.arrange',
#                   lapply(g_mean,
#                          levelplot, col.regions=terrain.colors(16), scales=list(draw=FALSE), 
#                          main='latent field mean', xlab='', ylab=''))
# ggsave(path = "USA_files", filename = "latent_field_mean_INLA.png", plot = plot_1,
#        units = "in", width = 10, height = 10)
# 
# 
# plot_2 <- do.call('grid.arrange',
#                   lapply(g_sd,
#                          levelplot, col.regions=terrain.colors(16), scales=list(draw=FALSE), 
#                          main='latent field SD', xlab='', ylab=''))
# ggsave(path = "USA_files", filename = "latent_field_sd_INLA.png", plot = plot_2,
#        units = "in", width = 10, height = 10)
# 
# for(l in valid_meshes_stelfi) {
#   plot_3 = stelfi::show_field(x = stelfi.fit[[l]]$par.random, smesh = mesh_n[[l]],
#                               border = usa_utm)
#   filename = paste("latent_field_stelfi_",l, ".png",sep="")
#   ggsave(path = "USA_files", filename = filename, plot = plot_3,
#          units = "in", width = 10, height = 10)
# }
# 
# dev.off()
# 
# ### --------------------------------------------------------------------------->> Final output CSV file
# 
# 
# # Record fitted values by INLA
# 
# list_fixed_INLA <- tempn
# list_hyperpar_INLA <- tempn
# list_dic <- tempn
# list_waic <- tempn
# 
# 
# for (k in 1:length(valid_meshes_INLA)){
#   f = valid_meshes_INLA[k]
#   list_fixed_INLA[[k]] <- pp.res[[f]]$summary.fixed
#   list_hyperpar_INLA[[k]] <- pp.res[[f]]$summary.hyperpar
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
# write.csv(fixed, file = "USA_files/summary_fixed_INLA.csv")
# 
# ### --------------------------------------------------------------------------->> CSV hyperparameters (INLA)
# fixed <- data.frame(matrix(unlist(list_hyperpar_INLA), ncol = 12, byrow = TRUE))
# colnames(fixed) <- rep(names(pp.res[[1]]$summary.hyperpar),2)
# fixed$Mesh <- paste("Mesh", valid_meshes_INLA, sep = "_")
# write.csv(fixed, file = "USA_files/summary_hyperpar_INLA.csv")
# 
# ### --------------------------------------------------------------------------->> CSV fitted values (stelfi)
# fixed <- data.frame(matrix(unlist(list_fixed_stelfi), ncol = 5, byrow = TRUE))
# colnames(fixed) <- names(stelfi.fit[[1]]$value)
# fixed$Mesh <- paste("Mesh", valid_meshes_stelfi, sep = "_")
# write.csv(fixed, file = "USA_files/summary_fixed_stelfi.csv")
# 
# ### --------------------------------------------------------------------------->> CSV value sd (stelfi)
# fixed <- data.frame(matrix(unlist(list_sd_stelfi), ncol = 5, byrow = TRUE))
# colnames(fixed) <- names(stelfi.fit[[1]]$sd)
# fixed$Mesh <- paste("Mesh", valid_meshes_stelfi, sep = "_")
# write.csv(fixed, file = "USA_files/summary_sd_stelfi.csv")
# 
# 
# ### --------------------------------------------------------------------------->> CSV DIC
# dic <- data.frame(matrix(unlist(list_dic), ncol = 1, byrow = TRUE))
# colnames(dic) <- "DIC"
# dic$Mesh <- paste("Mesh", valid_meshes_INLA, sep = "_")
# write.csv(dic, file = "USA_files/DIC.csv")
# 
# 
# ### --------------------------------------------------------------------------->> CSV WAIC
# waic <- data.frame(matrix(unlist(list_waic), ncol = 1, byrow = TRUE))
# colnames(waic) <- "WAIC"
# waic$Mesh <- paste("Mesh", valid_meshes_INLA, sep = "_")
# write.csv(waic, file = "USA_files/WAIC.csv")
# 
# ### --------------------------------------------------------------------------->> CSV INLA Computation time
# times <- data.frame(matrix(unlist(times.INLA), ncol = 1, byrow = TRUE))
# colnames(times) <- "CompTime"
# times$Mesh <- paste("Mesh", 1:nrow(mesh_mat), sep = "_")
# write.csv(times, file = "USA_files/Times_INLA.csv")
# 
# ### --------------------------------------------------------------------------->> CSV stelfi Computation time
# times <- data.frame(matrix(unlist(times.stelfi), ncol = 1, byrow = TRUE))
# colnames(times) <- "CompTime"
# times$Mesh <- paste("Mesh", 1:nrow(mesh_mat), sep = "_")
# write.csv(times, file = "USA_files/Times_stelfi.csv")

