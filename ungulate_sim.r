## read in "index" argument on NeSi
args <- commandArgs(trailingOnly = TRUE)
## index is the first element of this form slurm
## second for bash
k <- as.numeric(args[1])
## packages
library(sp)
library(INLA)
library(inlabru)
library(stelfi)
library(spatstat)
library(maptools) ## need for sp --> owin
## region data
url <- "https://gist.githubusercontent.com/cmjt/9a5e43cce69a3babb129b0a448e65752/raw/8653d11cb6ace667a3f5683982ae82366a3ac5d2/horse.csv"
horse <- read.csv(url)
sp <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(horse)), '0')))
bound <- inla.sp2segment(sp)
win <- as.owin(sp)

## Parameter mu can be regarded as a global mean
## level for the log intensity; i.e. the log-intensity fluctuates about it according to the spatial process
## Smaller values of the practical range will produce a spatial process that changes rapidly over the study window.
## Similarly, very large values of the practical range will produce an almost constant spatial process
## so that the log-intensity will be very close to  at all points of the study window.
## "practical" range is r = rho*sqrt(8) so that rho = r/sqrt(8) i.e., scale (rho) = prac_range (r) / sqrt(8)

######****************************#################
######**8 PARAMETER SETS***********#################
######****************************#################
sets <- expand.grid(data.frame(beta = c(5, 2)/748, ## per unit area
                               var = c(5, 0.1),
                               rho = c(5, 0.2)))
## stop sking to update
RandomFieldsUtils::RFoptions(install = "no")
## iterate over k for seeds
set.seed(4321 + k)
## simulate  8 different point pattern "sets"
points <- list()
for (i in 1:nrow(sets)) {
    X <- rLGCP("matern", mu = sets[i, 1],
               var = sets[i, 2], scale = sets[i, 3], nu = 1,
               win = win)
    Lamda <- attr(X, 'Lambda')
    points[[i]] <- data.frame(x = X$x, y = X$y)
}

## mesh construction (9 for now)
tmp <- lapply(seq(1, 4, length.out = 9), function(x) INLA::inla.mesh.2d(boundary = bound,
                                                                        max.edge = c(x, 2 * x)))
## mesh attributed
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
        t_bru <- system.time(fit_bru <- lgcp(cmp, locs,
                             samplers = sp,
                             domain = list(coordinates = tmp[[i]]),
                             options = list(control.inla = list(int.strategy = "eb")))
                        )
        results[(i + j - 1), 1:4] <- c(fit_bru$summary.fixed[, 1], fit_bru$summary.hyperpar[, 1], as.numeric(t_bru[3]))
        ## Stelfi fit
        t_stelfi <- system.time(fit_stelfi <- stelfi::fit_lgcp(locs = points[[j]], sp = sp, smesh = tmp[[i]]))
        results[(i + j - 1), 5:8] <- c(get_coefs(fit_stelfi)[c(1, 4:5), 1], as.numeric(t_stelfi[3]))
    }
}

write.csv(results, file = paste("res_", k, ".csv", sep = ""))
