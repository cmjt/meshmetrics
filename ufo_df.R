library(spatstat)
library(INLA)
library(patchwork)

ufo <- read.csv("data/ufo.csv")
names(ufo)
head(ufo)


library(sp)
## convert lat/long to UTM ...
dpts <- SpatialPoints(ufo[, c("city_longitude", "city_latitude")], proj4string=CRS("+proj=longlat +datum=WGS84"))
# spTransform(dpts, CRS("+proj=utm +zone=13 +ellps=WGS84"))

## create an interactive map using use a basemap given by the OpenStreetmap provider
library(tmap)
tmap_mode("view")
tm_basemap(leaflet::providers$OpenStreetMap) +
  tm_shape(dpts) + tm_dots()



## Model without covariates

library(maps) ## Display of maps
library(ggplot2)
library(mapproj) ## Converts latitude/longitude into projected coordinates
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


### Mesh construction

## source required functions
source("functions.r")
## define region and points
region <-  as(sf::st_as_sf(maps::map("usa", fill = TRUE, plot = FALSE)), "Spatial")
locs <- data.frame(x = ufo$city_longitude, y = ufo$city_latitude)

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

## plot mesh

mesh_list <- list(mesh01, mesh02, mesh03,
                  mesh04, mesh05, mesh06)

lst <- lapply(mesh_list, plt_mesh, xy = locs, domain = sf::st_as_sf(region))

((lst[[1]] + lst[[2]] + lst[[3]]) / (lst[[4]] + lst[[5]] + lst[[6]]))


## Use mesh05 to fit the model
plot(mesh05, asp=1)


## Points
ufo_pts <- as.matrix(ufo[, c("city_longitude", "city_latitude")])
mesh_pts <- as.matrix(mesh05$loc[, 1:2])
allpts <- rbind(mesh_pts, ufo_pts)

## Number of vertices in the mesh05
(nv <- mesh05$n)
## Number of points in the data
(n <- nrow(ufo_pts))

## Create SPDE
ufo_spde <- inla.spde2.pcmatern(mesh = mesh05, alpha = 2,
                                prior.range = c(0.1, 0.01), # Pro(range < 0.1) = 0.01
                                prior.sigma = c(0.01, 0.01)) # Pro(sigma > 0.01) = 0.01


## mesh weights
source("book.mesh.dual.R")

dmesh <- book.mesh.dual(mesh05)
usa_bdy <- cbind(usa[, c("long", "lat")])
domainSP <- SpatialPolygons(list(Polygons(list(Polygon(usa_bdy)), '0')))


library(rgeos)

domainSP <- gBuffer(domainSP, byid=TRUE, width=0) ## self-intersection

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
lmat <- inla.spde.make.A(mesh05, ufo_pts)
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





