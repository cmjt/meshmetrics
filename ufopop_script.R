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

dir.create("UFOPOP_stelfi_maps")
dir.create("UFOPOP_bru_maps")
dir.create("UFOPOP_meshes")


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


tmp <- lapply(seq(250, 1050, length.out = 10), function(x) INLA::inla.mesh.2d(loc.domain = bound$loc,
                                                                              max.edge = c(x, 2 * x)))

# create random mesh nodes
win <- as.owin.SpatialPolygons(usa_utm)
## choose resolution
set.seed(4321)
points <- list()
for (i in 1:10) {
  X <- rLGCP("matern", mu = -9.5,
             var = 0.1, scale = 100, nu = 1,
             win = win)
  Lamda <- attr(X, 'Lambda')
  points[[i]] <- data.frame(x = X$x, y = X$y)
}

tmp2 = lapply(seq(250, 1150, length.out = 10), function(x) INLA::inla.mesh.2d(loc.domain = bound$loc,
                                                                              loc = points[[(as.integer(0.01*x - 1.5))]],
                                                                              cutoff = 0.2*x,
                                                                              max.edge = c(x, 2 * x)))

tmp = append(tmp, tmp2)


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
#conti_US20 <- arrange(conti_US20, NAMELSAD)
#pop20 <- arrange(pop20, COUNTY)
conti_US20 <- conti_US20[order(conti_US20$NAMELSAD, conti_US20$STUSPS), ]
pop20 <- pop20[order(pop20$COUNTY, pop20$STATE_ABBR), ]
conti_US20$TOT_POP <- pop20$TOTAL_POPULATION
### create shape object with state polygons
us_pop20 <- conti_US20
us_pop20 <- as(us_pop20, "Spatial")
us_pop20 <- spTransform(us_pop20,
                        #CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs "))
                        CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs "))
us_pop20_sf <- st_as_sf(us_pop20)

# Population Density covariate for inlabru
# Treat population density as constant within the county
county_pops <- matrix(0, nrow=nrow(us_pop20_sf), ncol=1)
county_pops[,1] <- us_pop20_sf$TOT_POP
county_areas <- matrix(0, nrow=nrow(us_pop20_sf), ncol=1)
for (k in 1:nrow(us_pop20_sf)) {
  county_areas[k,1] = st_area(us_pop20_sf$geometry[[k]])
}

county_den <- county_pops / county_areas

pop_density <- function(x, y) {
  data = us_pop20_sf
  xy <- data.frame(x=x, y=y)
  xy <- sf::st_as_sf(xy, coords = c("x","y"))
  xy <- sf::st_geometry(xy)
  for (k in 1:nrow(data)) {
    county <- data$geometry[[k]]
    if (!all(st_is_valid(county))) { # check for invalid geometries
      county <- st_make_valid(county)
    }
    if (any(st_contains(county, xy,sparse=FALSE) > 0)) {
      return(county_den[k])
    }
  }
  return(0)
}


# Population density covariate for stelfi
# Treat population of each county as concentrated at N random pts
# N <- 8
# K <- nrow(us_pop20_sf)
# county_matrix <- matrix(0, nrow=N*K, ncol=3)
# for (k in 1:K) {
#   pts <- st_sample(us_pop20_sf$geometry[k], size = N)
#   for (n in 1:N) {
#     county_matrix[(N*(k-1) + n),1:2] <- pts[[n]]
#     county_matrix[(N*(k-1) + n),3] <- us_pop20_sf$TOT_POP[k] / N
#   }
# }

county_matrix <- read.csv("county_matrix.csv")[,2:4]

attrs <- lapply(tmp, stelfi:::meshmetrics)
results <- as.data.frame(matrix(0, nrow = length(tmp), ncol = 10))
names(results) <- c(paste("bru_", 1:4, sep = ""),"bru_time", paste("stelfi_", 1:4, sep = ""), "stelfi_time")

## append mesh attributes summaries
results$mean_radius_edge <- sapply(attrs, function(x) mean(x$triangles$radius_edge))
results$sd_radius_edge <- sapply(attrs, function(x) sd(x$triangles$radius_edge))
results$uq_radius_edge <- sapply(attrs, function(x) quantile(x$triangles$radius_edge)[4])
results$mean_radius_ratio <- sapply(attrs, function(x) mean(x$triangles$radius_ratio))
results$sd_radius_ratio <- sapply(attrs, function(x) sd(x$triangles$radius_ratio))
results$uq_radius_ratio <- sapply(attrs, function(x) quantile(x$triangles$radius_ratio)[4])
results$n_triangles <- sapply(tmp, function(x) x$n)

# correct data formats of locs for stelfi and inlabru
points <- locs
coordinates(locs) <- c("x", "y")

for (j in 1:length(tmp)) {
  # plot mesh
  filepath <- paste("UFOPOP_meshes/UFO_mesh_",j, ".jpg", sep="")
  jpeg(filepath, width = 1440, height = 1440)
  plot(tmp[[j]])
  lines(coordinates(usa_utm@polygons[[1]]@Polygons[[1]])[,1], coordinates(usa_utm@polygons[[1]]@Polygons[[1]])[,2], col = "red")
  dev.off()
  
  # plot mesh attributes
  filename <- paste("UFO_mesh_attrs_",j, ".jpg", sep="")
  mesh_plot <- stelfi::plot_mesh(tmp[[j]])
  ggsave(path = "UFOPOP_meshes", filename = filename, plot = mesh_plot,
         units = "in", width = 10, height = 10)
  
  # stelfi
  dmesh <- inla.mesh.dual(tmp[[j]])
  dmesh_areas <- numeric(length=tmp[[j]]$n)
  for (i in 1:tmp[[j]]$n) {
    coord <- sf::st_coordinates(sf::st_as_sf(dmesh[i,]))[,1:2]
    coord <- sf::st_polygon(list(coord))
    dmesh_areas[i] = st_area(coord)
  }
  xy <- data.frame(x=county_matrix[,1], y=county_matrix[,2])
  weights <- stelfi::points.in.mesh(xy, dmesh, county_matrix[,3])
  weights <- weights / dmesh_areas

  weights <- as.matrix(weights)
  out <- tryCatch(
    {
      t_stelfi <- system.time(fit_stelfi <- stelfi::fit_lgcp(locs = points, sp = usa_utm, smesh = tmp[[j]], covariates = weights))
      results[j, 6:10] <- c(get_coefs(fit_stelfi)[c(1, 2, 5:6), 1], as.numeric(t_stelfi[3]))
      res <- TMB::sdreport(fit_stelfi)
      field_plot <- stelfi::show_field(x = res$par.random, smesh = tmp[[j]], border = usa_utm)
      filename <- paste("UFO_field_", j, ".png", sep = "")
      ggsave(path = "UFOPOP_stelfi_maps", filename = filename, plot = field_plot,
             units = "in", width = 10, height = 10)
    },
    error = function(cond) {
      results[j, 6:10] <- rep(NA,5)
    }
  )

  # inlabru
  matern <- inla.spde2.pcmatern(tmp[[j]],
                              prior.sigma = c(0.1, 0.01),
                              prior.range = c(100, 0.01)
  )
 
  ecomp <- coordinates ~ elev(pop_density(x, y), model = "linear") +
    mySmooth(coordinates, model = matern) + Intercept(1)
  out <- tryCatch(
    {
      t_bru <- system.time(fit_bru <- lgcp(ecomp, locs,
                                           samplers = usa_utm,
                                           domain = list(coordinates = tmp[[j]]),
                                           options = list(control.inla = list(int.strategy = "eb")))
      )
      results[j, 1:5] <- c(fit_bru$summary.fixed[, 1], fit_bru$summary.hyperpar[, 1], as.numeric(t_bru[3]))
      field_plot <- stelfi::show_field(x=fit_bru$summary.random$mySmooth$mean, smesh = tmp[[j]], border = usa_utm)
      filename <- paste("UFO_field_", j, ".png", sep = "")
      ggsave(path = "UFOPOP_bru_maps", filename = filename, plot = field_plot,
             units = "in", width = 10, height = 10)
    },
    error = function(cond) {
      results[j, 1:5] <- rep(NA,5)
    }
  )
}

  write.csv(results, file = paste("res_ufopop", ".csv", sep = ""))