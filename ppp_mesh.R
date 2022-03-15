## ---- Build mesh based on point locations

## get the point locations for 8 different point patterns in one simulation

points_df <- lapply(
  unlist(unlist(sim[["points"]], recursive=FALSE), recursive=FALSE), 
  as.data.frame
)

## Mesh 1 and Mesh 2, with 8 different point patterns

mesh_ls_1 <- list()
mesh_ls_2 <- list()
mesh_ls_3 <- list()
mesh_ls_4 <- list()

for (l in 1:length(points_df)){
  mesh_ls_1[[l]] <- inla.mesh.create(loc = points_df[[l]], 
                                     boundary = inla.mesh.segment(domain))
  mesh_ls_2[[l]] <- inla.mesh.create(loc = points_df[[l]],
                                     boundary = inla.mesh.segment(domain),
                                     refine = list(max.edge = 0.5))
  mesh_ls_3[[l]] <- inla.mesh.create(loc = points_df[[l]],
                                     boundary = inla.mesh.segment(domain),
                                     refine = list(min.angle = 22, max.edge = 0.5),
                                     cutoff = 0.1)
  mesh_ls_4[[l]] <- inla.mesh.create(loc = points_df[[l]],
                                     boundary = inla.mesh.segment(domain),
                                     refine = list(min.angle = 22),
                                     quality.spec = list(segm = 0.2, loc = 0.05))
}

## put all meshes into one list

mesh_list <- list(mesh_1 = mesh_ls_1, 
                  mesh_2 = mesh_ls_2,
                  mesh_3 = mesh_ls_3,
                  mesh_4 = mesh_ls_4)


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
  allradius$mesh <- rep(paste0("Pattern", 1:length(mesh_ls)), times=sapply(radius, nrow))
  levels <- paste0("Pattern", 1:length(mesh_ls))
  allradius$meshes <- factor(allradius$mesh, levels = levels)
  return(list(mesh_attr = mesh_attr, allradius = allradius))
}

x <- lapply(mesh_list, mesh_info)

## create a file to save all the plots

dir.create("ppp_pdf_files")


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
  par(mfrow = c(3,1), mar = c(2, 2, 2, 2))
  lapply(mesh_list[[i]], 
         function(x) {
           plot(x, asp = 1, main=" ") 
           points(loc, col=2, pch = 16)
         }
  )
  dev.off()
  
  ## The radius-edge ratio R/l_min, 1/sqrt(3) for an equilateral triangle
  
  p1[[i]] <- ggplot(data=x[[i]][["allradius"]], aes(x=meshes, y=radius_edge, group=meshes, fill=meshes)) +
    geom_violin(position="dodge", alpha=0.5, trim = FALSE) +
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    geom_hline(yintercept = 1/sqrt(3), linetype="dashed", size=1) +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    labs(title = "Radius-edge") +
    theme_light() +
    theme(legend.position="none", axis.title.y=element_blank())
  
  ## The radius ratio r/R, 1/2 for an equilateral triangle.
  
  p2[[i]] <- ggplot(data=x[[i]][["allradius"]], aes(x=meshes,y=radius_ratio, group=meshes, fill=meshes)) +
    geom_violin(position="dodge", alpha=0.5, trim = FALSE) +
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    geom_hline(yintercept = 1/2, linetype="dashed",size=1) +
    scale_fill_brewer(palette="Paired") +
    coord_flip() +
    labs(title = "Radius-ratio") +
    theme_light() +
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
qpdf::pdf_combine(input = pdfs, output = "output.pdf")
setwd("..")

