## plotting packages
library(ggplot2)
library(patchwork)
## model fitting packages
library(INLA)
library(inlabru)
## read in data
bigfoot <- readr::read_csv("data/bigfoot.csv") 

usa <- map_data("state")
ggplot(data = usa, aes(x = long, y = lat, group = group)) + 
  geom_polygon(fill = "grey", alpha = 0.5) + 
  theme_void() +
  coord_fixed(1.3) +
   geom_point(data = bigfoot, aes(y = latitude,  x = longitude), alpha = 0.1,
               size = 1, inherit.aes = FALSE)

## mesh
region <-  sf::st_as_sf(maps::map("usa", fill = TRUE, plot = FALSE))
region_sp <- as(region, "Spatial")
locs <- data.frame(x = bigfoot$longitude, y = bigfoot$latitude)

## mesh examples (first three obvioulsy ridiculous mesh)

mesh01 <- inla.mesh.2d(boundary = inla.sp2segment(region_sp),
                       max.edge = c(7.5,15))

## add locs
mesh02 <- inla.mesh.2d(loc = locs, boundary = inla.sp2segment(region_sp),
                       max.edge = c(7.5,15))

## add cutoff
mesh03 <- inla.mesh.2d(loc = locs, boundary = inla.sp2segment(region_sp),
                       max.edge = c(0.75,3), cutoff = 1)

## model fit inlabru
mesh_list <- list(mesh01, mesh02, mesh03)
coordinates(locs) <- c("x", "y")
proj4string(locs) <- proj4string(region)
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
    fits[[i]] <- lgcp(cmp, data = locs,
                      samplers = data.frame(mesh_list[[i]]$loc[,1:2]),
                      domain = list(coordinates = mesh_list[[i]]),
                      options = list(control.inla = list(int.strategy = "eb"))
                      )
    ## Predict the spatial intensity surface
    lambda[[i]] <- predict(fits[[i]], pixels(mesh_list[[i]]), ~ exp(mySmooth + Intercept))
}
 
## Plot preds
l1 <- ggplot() +
    gg(lambda[[1]]) +
    theme_void() +
    scale_fill_viridis_b(na.value = "transparent") +
    coord_fixed()  +
    labs(fill = "Estimated surface")
l2 <- ggplot() +
    gg(lambda[[2]]) +
    theme_void() +
    scale_fill_viridis_b(na.value = "transparent") +
    ##gg(points, shape = 18, size = 1, alpha = 0.5) +
    coord_fixed()  +
    labs(fill = "Estimated surface")
l3 <- ggplot() +
    gg(lambda[[3]]) +
    theme_void() +
    scale_fill_viridis_b(na.value = "transparent") +
    ##gg(points, shape = 18, size = 1, alpha = 0.5) +
    coord_fixed()  +
    labs(fill = "Estimated surface")

l1 / l2 / l3

## model estimates

ests <- lapply(fits, function(x) rbind(x$summary.fixed[, 1:2], x$summary.hyperpar[, 1:2]))
ests

## mesh attributes
source("functions.r")
attrs <- lapply(mesh_list, get_triag_attributes)
sapply(attrs, function (x) mean(x$triangles$rr))

## estimated vs actual
n_actual <- nrow(bigfoot) ## 972
areakmsq <- sum(raster::area(region_sp))/1000^2 ## whole lot is 9147420 sqkm (inc ALaska and Hawai'i)
## spatial intensity
n_actual/areakmsq
