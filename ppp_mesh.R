## ---- Build mesh based on point locations

## Mesh 1 -- Mesh 12, with 8 different point patterns, total 96 meshes.

mesh_ls_1 <- list()
mesh_ls_2 <- list()
mesh_ls_3 <- list()
mesh_ls_4 <- list()
mesh_ls_5 <- list()
mesh_ls_6 <- list()
mesh_ls_7 <- list()
mesh_ls_8 <- list()
mesh_ls_9 <- list()
mesh_ls_10 <- list()
mesh_ls_11 <- list()
mesh_ls_12 <- list()

for (l in 1:length(points_df)){
  
  ## basic mesh created by the point locations and boundary constraint segments
  mesh_ls_1[[l]] <- inla.mesh.create(
    loc = points_df[[l]],
    boundary = inla.mesh.segment(domain),
  )
  
  ## adding parameter max.edge to refine the triangulation, 
  ## the maximum allowed edge length in the triangle is 0.1
  mesh_ls_2[[l]] <- inla.mesh.create(
    loc = points_df[[l]],
    boundary = inla.mesh.segment(domain),
    refine = list(max.edge = 0.3)
  )
  
  ## set the minimum allowed distance between points, to avoid small triangle 
  mesh_ls_3[[l]] <- inla.mesh.create(
    loc = points_df[[l]],
    boundary = inla.mesh.segment(domain),
    refine = list(min.angle = 30)
  )
  
  ## set the smallest triangle angle in the mesh 
  mesh_ls_4[[l]] <- inla.mesh.create(
    loc = points_df[[l]],
    boundary = inla.mesh.segment(domain),
    refine = list(max.edge = 0.5),
    quality.spec = list(loc = 0.1)
  )
  
  ## 
  mesh_ls_5[[l]] <- inla.mesh.create(
    loc = points_df[[l]],
    boundary = inla.mesh.segment(domain),
    refine = list(max.edge = 0.5),
    quality.spec = list(loc = 0.1, segm = 0.3)
  )
  
  ##
  mesh_ls_6[[l]] <- inla.mesh.create(
    loc = points_df[[l]],
    boundary = inla.mesh.segment(domain),
    refine = list(max.edge = 0.2, min.angle = 30)
  )
  
  mesh_ls_7[[l]] <- inla.mesh.create(
    loc = points_df[[l]],
    boundary = inla.mesh.segment(domain),
    refine = list(max.edge = 0.2, min.angle = 30),
    cutoff = 0.01
  )
  
  mesh_ls_8[[l]] <- inla.mesh.create(
    loc = points_df[[l]],
    boundary = inla.mesh.segment(domain),
    refine = list(max.edge = 0.2, min.angle = 30),
    quality.spec = list(loc = 0.2, segm = 0.3),
    cutoff = 0.01
  )
  
  ## mesh created by function `inla.mesh.2d()` 
  
  mesh_ls_9[[l]] <- inla.mesh.2d(loc = points_df[[l]], boundary = inla.mesh.segment(domain), 
                                 max.edge = c(0.5, 0.5))
  
  mesh_ls_10[[l]]  <- inla.mesh.2d(loc = points_df[[l]], boundary = inla.mesh.segment(domain), 
                                   offset = c(0.5, 1), max.edge = c(0.5, 0.5))
  
  mesh_ls_11[[l]]  <- inla.mesh.2d(loc = points_df[[l]], boundary = inla.mesh.segment(domain), 
                                   offset = c(0.5, 1), max.edge = c(0.2, 0.5),cutoff = 0.01)
  
  mesh_ls_12[[l]]  <- inla.mesh.2d(loc = points_df[[l]], boundary = inla.mesh.segment(domain), 
                                   offset = c(0.5, 1), max.edge = c(0.2, 0.5), min.angle = 30)
}

## put all meshes into one list

mesh_list <- list(Mesh_band_1 = mesh_ls_1, 
                  Mesh_band_2 = mesh_ls_2,
                  Mesh_band_3 = mesh_ls_3,
                  Mesh_band_4 = mesh_ls_4,
                  Mesh_band_5 = mesh_ls_5,
                  Mesh_band_6 = mesh_ls_6,
                  Mesh_band_7 = mesh_ls_7,
                  Mesh_band_8 = mesh_ls_8,
                  Mesh_band_9 = mesh_ls_9,
                  Mesh_band_10 = mesh_ls_10,
                  Mesh_band_11 = mesh_ls_11,
                  Mesh_band_12 = mesh_ls_12)


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

x <- lapply(mesh_list, mesh_info)


## CSV file mesh_attr

lapply(1:length(x), function(i) write.csv(x[[i]]$allradius, 
                                          file = paste0(names(x[i]), ".csv"),
                                          row.names = FALSE))



p1 <- list()
p2 <- list()
p3 <- list()
rr <- list()
rr_pl <- list()
re <- list()
re_pl <- list()


for (i in seq_along(mesh_list)) {
  
  ## plot the mesh
  
  pdf(paste0("ppp_pdf_files/mesh_plot_mesh", i, ".pdf"), 
      width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4") 
  par(mfrow = c(4,2), mar = c(2, 2, 2, 2))
  lapply(mesh_list[[i]], 
         function(x) {
           plot(x, asp = 1, main=" ")
         }
  )
  dev.off()
  
  
  ## The radius-edge ratio R/l_min, 1/sqrt(3) for an equilateral triangle
  
  at <- seq(1, 12)*(1/sqrt(3))
  L <- parse(text = paste("frac(", seq(1, 12), ", sqrt(3))", sep = ""))
  
  p1[[i]] <- ggplot(data=x[[i]][["allradius"]], aes(x=meshes, y=radius_edge, group=meshes, fill=meshes)) +
    geom_violin(fill = "grey") +
    geom_boxplot(width=0.1) +
    geom_hline(yintercept = 1/sqrt(3), linetype = "dashed", size = 1) +
    scale_y_continuous(breaks = at, labels = L) + 
    coord_flip() +
    labs(title = paste0("Radius-edge ration, Mesh-band", i)) +
    theme_classic() +
    theme(legend.position="none", axis.title.y=element_blank())
  
  
  ## The radius ratio r/R, 1/2 for an equilateral triangle.
  
  at <- seq(0, 1/2, by = 0.05)
  L <- parse(text = paste("frac(", seq(0, 10), ", 20)", sep = ""))
  
  p2[[i]] <- ggplot(data=x[[i]][["allradius"]], aes(x=meshes,y=radius_ratio, group=meshes, fill=meshes)) +
    geom_violin(fill = "grey") +
    geom_boxplot(width=0.1) +
    geom_hline(yintercept = 1/2, linetype = "dashed", size = 1) +
    scale_y_continuous(breaks = at, labels = L) + 
    coord_flip() +
    labs(title = paste0("Radius ratio, Mesh-band", i)) +
    theme_classic() +
    theme(legend.position="none", axis.title.y=element_blank())
  
  p3[[i]] <- p1[[i]] / p2[[i]]
  ggsave(path = "ppp_pdf_files", filename = paste0("ratio_mesh", i, ".pdf"), plot = p3[[i]],
         width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4")
  
  rr[[i]] <- lapply(x[[i]][["mesh_attr"]], rr_func)
  rr_pl[[i]] <- grid.arrange(grobs=rr[[i]], ncol=2)
  ggsave(path = "ppp_pdf_files", filename = paste0("radius_ratio_mesh", i, ".pdf"), plot = rr_pl[[i]], 
         width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4") 
  
  re[[i]] <- lapply(x[[i]][["mesh_attr"]], re_func)
  re_pl[[i]] <- grid.arrange(grobs=re[[i]], ncol=2)
  ggsave(path = "ppp_pdf_files", filename = paste0("radius_edge_mesh", i, ".pdf"), plot = re_pl[[i]], 
         width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4")
  
  dev.off()
}

## combine pdf files

setwd("ppp_pdf_files/")
pdfs <- list.files(pattern = "\\.pdf$")
qpdf::pdf_combine(input = pdfs, output = "ppp_output.pdf")
setwd("..")

