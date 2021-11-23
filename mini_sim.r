source("functions.r")
library(sp)
library(INLA)
library(inlabru)
library(spatstat)
library(maptools) ## need for sp --> owi
require(ggplot2) ## plotting
require(patchwork)
## region data
url <- "https://gist.githubusercontent.com/cmjt/9a5e43cce69a3babb129b0a448e65752/raw/8653d11cb6ace667a3f5683982ae82366a3ac5d2/horse.csv"
horse <- read.csv(url)
sp <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(horse)), '0')))
bound <- inla.sp2segment(sp)
win <- as.owin(sp)
## choose resolution
npix <- 300
spatstat.options(npixel = npix)
## simulate one realisation
## Parameter mu can be regarded as a global mean
## level for the log intensity; i.e. the log-intensity fluctuates about it according to the spatial process
## Smaller values of the practical range will produce a spatial process that changes rapidly over the study window.
## Similarly, very large values of the practical range will produce an almost constant spatial process
## so that the log-intensity will be very close to  at all points of the study window.
## "practical" range is r = rho*sqrt(8) so that rho = r/sqrt(8) i.e., scale (rho) = prac_range (r) / sqrt(8)
set.seed(4321)
X <- rLGCP("matern", mu = 2,
           var = 2, scale = 0.7, nu = 1,
           win = win)
## plot simulated rf
Lamda <- attr(X, 'Lambda')
points <- data.frame(x = X$x, y = X$y)
rf <- log(Lamda$v)
rf <- apply(rf, 2, t) ## field flips
ex <- expand.grid(seq(win$xrange[1], win$xrange[2],length.out = npix),
                  seq(win$yrange[1], win$yrange[2],length.out = npix))
dat <- data.frame(x = ex$Var1, y = ex$Var2, rf = as.vector(t(rf)))
ins <- inside.owin(dat$x, dat$y, win)
dat$rf[!ins] <- NA
true <- ggplot(dat , aes(x = x, y = y, fill = rf)) +
    geom_tile() +
    labs(fill = "Expected (log) \nnumber \nof points") +
    geom_point(data = points, aes(x = x, y = y), inherit.aes = FALSE,
               shape = 18, size = 1, alpha = 0.5) +
    scale_fill_viridis_b(na.value = "transparent") +
    theme_void()  +
    geom_sf(data = sf::st_as_sf(sp), size = 1,
            inherit.aes = FALSE, fill = "transparent") +
    ggtitle("Simulated field and points")

## mesh construction
m1 <- INLA::inla.mesh.2d(boundary = bound,
                         max.edge = 0.3)
m2 <- INLA::inla.mesh.2d(boundary = bound,
                         max.edge = 3)
m3 <- INLA::inla.mesh.create(horse, bound = bound)

mesh_list <- list(m1, m2, m3)

## inlabru for model fitting
locs <- points
coordinates(locs) <- c("x", "y")
fits <- lambda <- list()

for (i in 1:length(mesh_list)){
    ## Define SPDE prior
    matern <- INLA::inla.spde2.pcmatern(mesh_list[[i]],
                                        prior.sigma = c(0.1, 0.01),
                                        prior.range = c(0.01, 0.01)
                                        )

    ## Define domain of the LGCP as well as the model components (spatial SPDE
    ## effect and Intercept)
    cmp <- coordinates ~ mySmooth(main = coordinates, model = matern) + Intercept
    ## Fit the model (with int.strategy = "eb" to make the example take less time)
    fits[[i]] <- lgcp(cmp, locs,
                samplers = sp,
                domain = list(coordinates = mesh_list[[i]]),
                options = list(control.inla = list(int.strategy = "eb"))
                )
    ## Predict the spatial intensity surface
    lambda[[i]] <- predict(fits[[i]], pixels(mesh_list[[i]]), ~ mySmooth + Intercept)
}
## parameter estimates

ests <- lapply(fits, function(x) rbind(x$summary.fixed[, 1:2], x$summary.hyperpar[, 1:2]))


## Plot the predicted intensity
l1 <- ggplot() +
    gg(lambda[[1]]) +
    theme_void() +
    scale_fill_viridis_b(na.value = "transparent") +
    coord_fixed()  +
    labs(fill = "Predicted (log) \nnumber \nof points")
mp1 <- ggplot() + gg(m1) + theme_void() + gg(locs, shape = 18, size = 1, alpha = 0.5)
l2 <- ggplot() +
    gg(lambda[[2]]) +
    theme_void() +
    scale_fill_viridis_b(na.value = "transparent") +
    ##gg(points, shape = 18, size = 1, alpha = 0.5) +
    coord_fixed()  +
    labs(fill = "Predicted (log) \nnumber \nof points")
mp2 <- ggplot() + gg(m2) + theme_void() + gg(locs, shape = 18, size = 1, alpha = 0.5)
l3 <- ggplot() +
    gg(lambda[[3]]) +
    theme_void() +
    scale_fill_viridis_b(na.value = "transparent") +
    ##gg(points, shape = 18, size = 1, alpha = 0.5) +
    coord_fixed()  +
    labs(fill = "Predicted (log) \nnumber \nof points")
mp3 <- ggplot() + gg(m3) + theme_void() + gg(locs, shape = 18, size = 1, alpha = 0.5)



## mesh attributes

attrs <- lapply(mesh_list, get_triag_attributes)
sapply(attrs, function (x) mean(x$triangles$rr))

rr1 <- ggplot(attrs[[1]]$sf, aes(fill = attrs[[1]]$triangles$rr)) +
    geom_sf(colour = NA) + theme_void() +
    scale_fill_continuous(name = "Radius Ratio", limits = c(0, 0.5)) +
    labs(title = paste("Mean radius ratio = ", round(mean(attrs[[1]]$triangles$rr), 4)))

rr2 <- ggplot(attrs[[2]]$sf, aes(fill = attrs[[2]]$triangles$rr)) +
    geom_sf(colour = NA) + theme_void() +
    scale_fill_continuous(name = "Radius Ratio", limits = c(0, 0.5))+
    labs(title = paste("Mean radius ratio = ", round(mean(attrs[[2]]$triangles$rr), 4)))

rr3 <- ggplot(attrs[[3]]$sf, aes(fill = attrs[[3]]$triangles$rr)) +
    geom_sf(colour = NA) + theme_void() +
    scale_fill_continuous(name = "Radius Ratio", limits = c(0, 0.5)) +
    labs(title = paste("Mean radius ratio = ", round(mean(attrs[[3]]$triangles$rr), 4)))



re1 <- ggplot(attrs[[1]]$sf, aes(fill = attrs[[1]]$triangles$re)) +
    geom_sf(colour = NA) + theme_void() +
    scale_fill_continuous(name = "Radius-edge Ratio") +
    labs(title = paste("Mean radius-edge ratio = ", round(mean(attrs[[1]]$triangles$re), 4)))

re2 <- ggplot(attrs[[2]]$sf, aes(fill = attrs[[2]]$triangles$re)) +
    geom_sf(colour = NA) + theme_void() +
    scale_fill_continuous(name = "Radius-edge Ratio")+
    labs(title = paste("Mean radius-edge ratio = ", round(mean(attrs[[2]]$triangles$re), 4)))

re3 <- ggplot(attrs[[3]]$sf, aes(fill = attrs[[3]]$triangles$re)) +
    geom_sf(colour = NA) + theme_void() +
    scale_fill_continuous(name = "Radius-edge Ratio") +
    labs(title = paste("Mean radius-edge ratio = ", round(mean(attrs[[3]]$triangles$re), 4)))

mesh_plots <- (mp1 + mp2 + mp3)  + plot_annotation(tag_levels = 'A')

radius_edge <- (re1 / re2 / re3)  + plot_annotation(tag_levels = 'A')

radius_ratio <- (rr1 / rr2 / rr3)  + plot_layout(guides = "collect")  + plot_annotation(tag_levels = 'A')

predictions <- (l1 / l2 / l3)  + plot_annotation(tag_levels = 'A')

tab <-  do.call('rbind', ests)
rownames(tab) <- paste(c("Intercept","Range", "Sigma"), rep(LETTERS[1:3], each = 3))
actual <- data.frame(mu = 2, Sigma = sqrt(2), Range = 0.7*sqrt(8))
rownames(actual) <- "True"
layout <- "
AAAA
AAAA
BB#C
BB##
"
fitted <- (true + gridExtra::tableGrob(round(tab, 3)) +  gridExtra::tableGrob(round(actual, 3))) + 
  plot_layout(design = layout)

mesh_plots
radius_edge
radius_ratio
predictions
fitted
