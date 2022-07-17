source("mesh_sim.R")
ufo_coo <- coordinates(ufo_sp)
mesh_val <- c(0, "usa_utm", c(750000, 1500000), 0,
              "ufo_coo", "usa_utm", c(750000, 1500000), 0,
              0, "usa_utm", c(750000, 1500000), 1000000,
              "ufo_coo", "usa_utm", c(150000,200000), 300000,
              "ufo_coo", "usa_utm", c(75000,200000), 50000,
              "ufo_coo", "usa_utm", c(75000,100000), 50000)
mesh_mat <- matrix(mesh_val, ncol=5, byrow = TRUE)
mesh_n <- mesh_sim(mesh_mat)


## get triangulation attributes of each mesh

mesh_info <- function(mesh_ls){
  
  mesh_attr <- lapply(mesh_ls, get_triag_attributes)
  angles <- vector("list", length = length(mesh_ls))
  radius <- vector("list", length = length(mesh_ls))
  
  for (x in 1:length(mesh_ls)){
    radius[[x]] <- data.frame(radius_edge=mesh_attr[[x]]$triangles$re, 
                              radius_ratio=mesh_attr[[x]]$triangles$rr)
    angles[[x]] <- mesh_attr[[x]]$angles
  }
  
  allradius <- bind_rows(radius)
  allradius$mesh <- rep(paste0("Mesh", 1:length(mesh_ls)), times=sapply(radius, nrow))
  levels <- paste0("Mesh", 1:length(mesh_ls))
  allradius$meshes <- factor(allradius$mesh, levels = levels)
  return(list(mesh_attr = mesh_attr, allradius = allradius))
}


x <- mesh_info(mesh_n)

## Write the mean and sd of the meshes to CSV
rr_mean <- lapply(x$mesh_attr, function (y) mean(y$triangles$rr))
rr_sd <- lapply(x$mesh_attr, function (y) sd(y$triangles$rr))
re_mean <- lapply(x$mesh_attr, function (y) mean(y$triangles$re))
re_sd <- lapply(x$mesh_attr, function (y) sd(y$triangles$re))
n_triangles <- lapply(mesh_n, function (y) y$n)

mesh_attrs <- vector("list", length = 5)
mesh_attrs[[1]] <- rr_mean
mesh_attrs[[2]] <- rr_sd
mesh_attrs[[3]] <- re_mean
mesh_attrs[[4]] <- re_sd
mesh_attrs[[5]] <- n_triangles

## create a file to save all the plots
dir.create("ufo_sim_files")


## functions to visual the mesh attributes

rr_func <- function(mesh_attr){
  rr <- ggplot(mesh_attr$sf, aes(fill = mesh_attr$triangles$rr)) +
    geom_sf(colour = NA) + 
    theme_void() +
    theme(legend.title = element_text(size = 8),
          legend.key.size = unit(5, "mm")) +
    scale_fill_continuous(name = "Radius Ratio", limits = c(0, 0.5))
}

re_func <- function(mesh_attr){
  re <- ggplot(mesh_attr$sf, aes(fill = mesh_attr$triangles$re)) +
    geom_sf(colour = NA) + 
    theme_void() +  
    theme(legend.title = element_text(size = 8),
          legend.key.size = unit(5, "mm")) +
    scale_fill_continuous(name = "Radius-edge ratio", limits = c(1/sqrt(3), NA))
}



## plot the mesh

pdf(paste0("ufo_sim_files/ufo_mesh.pdf"), width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4") 
par(mfrow = c(3,2), mar = c(2, 2, 2, 2))
lapply(mesh_n, function(x) {plot(x, asp = 1, main=" ")})
dev.off()

## The radius-edge ratio R/l_min, 1/sqrt(3) for an equilateral triangle

p1 <- ggplot(data=x[["allradius"]], aes(x=meshes, y=radius_edge, group=meshes, fill=meshes)) +
  geom_violin(position="dodge", alpha=0.5, trim = FALSE) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  geom_hline(yintercept = 1/sqrt(3), linetype="dashed", size=1) +
  coord_flip() +
  scale_fill_brewer(palette="Paired") +
  labs(title = "Radius-edge") +
  theme_light() +
  theme(legend.position="none", axis.title.y=element_blank())

## The radius ratio r/R, 1/2 for an equilateral triangle.

p2 <- ggplot(data=x[["allradius"]], aes(x=meshes,y=radius_ratio, group=meshes, fill=meshes)) +
  geom_violin(position="dodge", alpha=0.5, trim = FALSE) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  geom_hline(yintercept = 1/2, linetype="dashed",size=1) +
  scale_fill_brewer(palette="Paired") +
  coord_flip() +
  labs(title = "Radius-ratio") +
  theme_light() +
  theme(legend.position="none", axis.title.y=element_blank())

p3 <- p1 / p2
ggsave(path = "ufo_sim_files", filename = "violin_ratio.pdf", plot = p3,
       width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4")

rr <- lapply(x[["mesh_attr"]], rr_func)
rr_pl <- grid.arrange(grobs=rr, ncol=2)
ggsave(path = "ufo_sim_files", filename = "radius_ratio.pdf", plot = rr_pl, 
       width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4") 

re <- lapply(x[["mesh_attr"]], re_func)
re_pl <- grid.arrange(grobs=re, ncol=2)
ggsave(path = "ufo_sim_files", filename = "radius_edge.pdf", plot = re_pl, 
       width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4")

dev.off()


## combine pdf files

setwd("ufo_sim_files/")
pdfs <- list.files(pattern = "\\.pdf$")
qpdf::pdf_combine(input = pdfs, output = "output.pdf")
setwd("..")
