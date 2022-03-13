## ---- Build mesh based on point locations

## get the point locations for 8 different point patterns in one simulation

points_df <- lapply(
  unlist(unlist(sim[["points"]], recursive=FALSE), recursive=FALSE), 
  as.data.frame
)

## Mesh 1 with 8 different point patterns

mesh_ls <- list()
for (l in 1:8){
  mesh_ls[[l]] <- inla.mesh.create(loc = points_df[[l]], 
                                   boundary = inla.mesh.segment(domain))
}

## function to plot mesh attributes

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


mesh_info <- function(mesh_ls){
  
  ## plot the mesh
  
  pdf("mesh_plot.pdf", width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4")
  par(mfrow = c(2,4), mai = c(0, 0, 0, 0))
  plot_ls <- lapply(mesh_ls, 
                    function(x) {
                      plot(x, asp = 1, main="") 
                    }
  )
  dev.off()
  
  ## get the mesh attributes
  
  mesh_attr <- lapply(mesh_ls, get_triag_attributes)
  angles <- vector("list", length = length(mesh_ls))
  radius <- vector("list", length = length(mesh_ls))
  for(x in 1:length(mesh_ls)){
    radius[[x]] <- data.frame(radius_edge=mesh_attr[[x]]$triangles$re, 
                              radius_ratio=mesh_attr[[x]]$triangles$rr)
    angles[[x]] <- mesh_attr[[x]]$angles
  }
  
  allradius <- bind_rows(radius)
  allradius$mesh <- rep(paste0("Mesh", 1:length(mesh_ls)), times=sapply(radius, nrow))
  levels <- paste0("Mesh", 1:length(mesh_ls))
  allradius$meshes <- factor(allradius$mesh, levels = levels)
  
  ## plot the attributes
  
  ## The radius-edge ratio R/l_min, 1/sqrt(3) for an equilateral triangle
  
  p1 <- ggplot(data=allradius, aes(x=meshes, y=radius_edge, group=meshes, fill=meshes)) +
    geom_violin(position="dodge", alpha=0.5, trim = FALSE) +
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    geom_hline(yintercept = 1/sqrt(3), linetype="dashed", size=1) +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    labs(title = "Radius-edge") +
    theme_light() +
    theme(legend.position="none", axis.title.y=element_blank())
  
  
  ## The radius ratio r/R, 1/2 for an equilateral triangle.
  
  p2 <- ggplot(data=allradius, aes(x=meshes,y=radius_ratio, group=meshes, fill=meshes)) +
    geom_violin(position="dodge", alpha=0.5, trim = FALSE) +
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    geom_hline(yintercept = 1/2, linetype="dashed",size=1) +
    scale_fill_brewer(palette="Paired") +
    coord_flip() +
    labs(title = "Radius-ratio") +
    theme_light() +
    theme(legend.position="none", axis.title.y=element_blank())
  
  p3 <- p1/p2
  ggsave(path = "ppp_pdf_files", filename = "ratio.pdf", plot = p3, 
         width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4")
  
  rr <- lapply(mesh_attr, rr_func)
  rr_pl <- grid.arrange(grobs=rr, ncol=2)
  ggsave(path = "ppp_pdf_files", filename = "radius_ratio.pdf", plot = rr_pl, 
         width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4") 
  
  re <- lapply(mesh_attr, re_func)
  re_pl <- grid.arrange(grobs=re, ncol=2)
  ggsave(path = "ppp_pdf_files", filename = "radius_edge.pdf", plot = re_pl, 
         width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4") 
  dev.off()
}

mesh_info(mesh_ls)


