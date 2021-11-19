## packages
library(sf)
library(INLA)
require(ggplot2)
require(patchwork) ## just for aligning plots
## read in data
ufo <- readr::read_csv("data/ufo.csv")
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

lst <- lapply(mesh_list,
              plt_mesh, xy = locs, domain = sf::st_as_sf(region))

((lst[[1]] + lst[[2]] + lst[[3]]) / (lst[[4]] + lst[[5]] + lst[[6]])) + plot_annotation(tag_levels = 'A')


## mesh attributes

attributes <- lapply(mesh_list, get_triag_attributes) ## takes a while

all_r <-  do.call(`rbind`, s <- lapply(attributes, function(x) x$triangles))
all_r$mesh <- rep(paste("Mesh", 1:6, sep = " "), times = sapply(s, nrow))
## data manipulation
require(tidyverse)
all_r <- all_r %>%
        group_by(mesh) %>%
        mutate(q_rr = ifelse(rr <= quantile(rr , probs = c(0.05)), 0, 1)) %>%
        mutate(q_re = ifelse(rr <= quantile(rr , probs = c(0.95)), 0, 1))


at <- c(0, 1/4, 1/2)
L <- parse(text = paste("frac( 1,", c(4, 2),")", sep = ""))
r_r <- ggplot(all_r, aes(y = rr, x = mesh)) +
    geom_violin(fill = "white") +
    geom_violin(data = all_r[all_r$q_rr == 1, ],fill = "grey") +
    geom_boxplot(width = 0.1) +
    theme_classic() +
    geom_hline(yintercept = 1/2) +
    scale_y_continuous(breaks = at, labels  = c(0, L)) +
    ggtitle("Radius ratio") + ylab("") + xlab("")

at <- seq(1,3, by = 1)*(1/sqrt(3))
L <- parse( text = paste("frac(", seq(1,3, by = 1), ", sqrt(3))", sep = ""))

r_e <- ggplot(all_r, aes(y = re, x = mesh)) +
    geom_violin(fill = "white") +
    scale_y_continuous(breaks = at, labels  = L) + 
    geom_boxplot(width = 0.1) +
    theme_classic() +
    geom_hline(yintercept = 1/sqrt(3)) + 
    ggtitle("Radius-edge ratio") + ylab("") + xlab("")

r_r + r_e

