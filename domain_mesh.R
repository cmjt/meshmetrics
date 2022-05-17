
dir.create("domain_pdf_files")

## ========================================================================= Mesh Construction

mesh1 <- inla.mesh.2d(loc.domain = domain, max.edge = 0.5)

mesh2 <- inla.mesh.2d(loc.domain = domain, max.edge = 0.5, min.angle = 30)

mesh3 <- inla.mesh.2d(loc.domain = domain, max.edge = 0.5, min.angle = 35)

mesh4 <- inla.mesh.2d(loc.domain = domain, max.edge = 0.2)

mesh5 <- inla.mesh.2d(loc.domain = domain, max.edge = 0.2, offset = 0.1)

mesh6 <- inla.mesh.2d(loc.domain = domain, max.edge = c(0.2,0.5), offset = c(0.1, 0.3))

mesh7 <- inla.mesh.2d(loc.domain = domain, max.edge = c(0.2,0.5), offset = c(0.1, 0.3), min.angle = c(32, 26)) 

mesh8 <- inla.mesh.2d(loc.domain = domain, max.edge = c(0.2, 0.5), offset = c(0.05, 0.3))

mesh9 <- inla.mesh.2d(loc.domain = domain, max.edge = c(0.2, 0.5), offset = c(0.05, 0.3), min.angle = c(32, 26))


## put all meshes into one list
mesh_list_1 <- list(
  m1 = mesh1,
  m2 = mesh2,
  m3 = mesh3,
  m4 = mesh4,
  m5 = mesh5,
  m6 = mesh6,
  m7 = mesh7,
  m8 = mesh8,
  m9 = mesh9)


mesh10 <- inla.mesh.2d(loc.domain = domain, max.edge = c(0.1, 0.2), min.angle = c(32, 26))

mesh11 <- inla.mesh.2d(loc.domain = domain, max.edge = c(0.1, 0.2), offset = c(0.1, 0.5))

mesh12 <- inla.mesh.2d(loc.domain = domain, max.edge = c(0.1, 0.2), offset = c(0.05, 0.2))

mesh13 <- inla.mesh.2d(loc.domain = domain, max.edge = c(0.1, 0.2), offset = c(0.05, 0.2), min.angle = c(32, 26))

mesh14 <- inla.mesh.2d(loc.domain = domain, max.edge = c(0.1, 0.2), offset = c(0.3, 0.8))

mesh15 <- inla.mesh.2d(loc.domain = domain, max.edge = c(0.1, 0.2), min.angle = c(34, 30))

mesh16 <- inla.mesh.2d(loc.domain = domain, max.edge = c(0.05, 0.1), offset = c(0.05, 0.2), min.angle = c(32, 26))

mesh17 <- inla.mesh.2d(loc.domain = domain, max.edge = c(0.05, 0.1), offset = c(0.05, 0.3), min.angle = c(34, 26))

mesh18 <- inla.mesh.2d(loc.domain = domain, max.edge = c(0.05, 0.1), offset = c(0.05, 0.3), min.angle = c(34, 30))


mesh_list_2 <- list(
  m10 = mesh10,
  m11 = mesh11,
  m12 = mesh12,
  m13 = mesh13,
  m14 = mesh14,
  m15 = mesh15,
  m16 = mesh16,
  m17 = mesh17,
  m18 = mesh18)

## ========================================================================= Mesh Assessment

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
    scale_fill_continuous(name = "Radius-edge ratio", limits = c(1/sqrt(3), 1.4))
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

x <- mesh_info(mesh_list_1)

y <- mesh_info(mesh_list_2)


## plot the mesh
## group 1
pdf("domain_pdf_files/mesh_plot_1.pdf", 
  width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4")
par(mfrow = c(3,3), mar = c(2, 2, 2, 2))
lapply(mesh_list_1, 
  function(x) {
    plot(x, asp = 1, main=" ")
  }
)

## group 2
pdf("domain_pdf_files/mesh_plot_2.pdf", 
  width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4")
par(mfrow = c(3,3), mar = c(2, 2, 2, 2))
lapply(mesh_list_2, 
  function(x) {
    plot(x, asp = 1, main=" ")
  }
)


## the violin plot 

## The radius-edge ratio R/l_min, 1/sqrt(3) for an equilateral triangle
at <- seq(1, 12)*(1/sqrt(3))
L <- parse(text = paste("frac(", seq(1, 12), ", sqrt(3))", sep = ""))
x1 <- ggplot(data=x[["allradius"]], aes(x=meshes, y=radius_edge, group=meshes, fill=meshes)) +
  geom_violin(fill = "grey") +
  geom_boxplot(width=0.1) +
  geom_hline(yintercept = 1/sqrt(3), linetype = "dashed", size = 1) +
  scale_y_continuous(breaks = at, labels = L) + 
  coord_flip() +
  labs(title = "Radius-edge ratio") +
  theme_classic() +
  theme(legend.position="none", axis.title.y=element_blank())

y1 <- ggplot(data=y[["allradius"]], aes(x=meshes, y=radius_edge, group=meshes, fill=meshes)) +
  geom_violin(fill = "grey") +
  geom_boxplot(width=0.1) +
  geom_hline(yintercept = 1/sqrt(3), linetype = "dashed", size = 1) +
  scale_y_continuous(breaks = at, labels = L) + 
  coord_flip() +
  labs(title = "Radius-edge ratio") +
  theme_classic() +
  theme(legend.position="none", axis.title.y=element_blank())


## The radius ratio r/R, 1/2 for an equilateral triangle.
at <- seq(0, 1/2, by = 0.05)
L <- parse(text = paste("frac(", seq(0, 10), ", 20)", sep = ""))
x2 <- ggplot(data=x[["allradius"]], aes(x=meshes,y=radius_ratio, group=meshes, fill=meshes)) +
  geom_violin(fill = "grey") +
  geom_boxplot(width=0.1) +
  geom_hline(yintercept = 1/2, linetype = "dashed", size = 1) +
  scale_y_continuous(breaks = at, labels = L) + 
  coord_flip() +
  labs(title = "Radius ratio") +
  theme_classic() +
  theme(legend.position="none", axis.title.y=element_blank())

y2 <- ggplot(data=x[["allradius"]], aes(x=meshes,y=radius_ratio, group=meshes, fill=meshes)) +
  geom_violin(fill = "grey") +
  geom_boxplot(width=0.1) +
  geom_hline(yintercept = 1/2, linetype = "dashed", size = 1) +
  scale_y_continuous(breaks = at, labels = L) + 
  coord_flip() +
  labs(title = "Radius ratio") +
  theme_classic() +
  theme(legend.position="none", axis.title.y=element_blank())

violin_x <- x1 / x2
ggsave(path = "domain_pdf_files", filename = "violin_mesh_1.pdf", plot = violin_x,
  width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4")

violin_y <- y1 / y2
ggsave(path = "domain_pdf_files", filename = "violin_mesh_2.pdf", plot = violin_y,
  width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4")


re_1 <- lapply(x[["mesh_attr"]], re_func)
re_x <- grid.arrange(grobs=re_1, ncol=3)
ggsave(path = "domain_pdf_files", filename = "radius_edge_1.pdf", plot = re_x, 
  width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4")

rr_1 <- lapply(x[["mesh_attr"]], rr_func)
rr_x <- grid.arrange(grobs=rr_1, ncol=3)
ggsave(path = "domain_pdf_files", filename = "radius_ratio_1.pdf", plot = rr_x, 
  width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4") 



re_2 <- lapply(y[["mesh_attr"]], re_func)
re_y <- grid.arrange(grobs=re_2, ncol=3)
ggsave(path = "domain_pdf_files", filename = "radius_edge_2.pdf", plot = re_y, 
  width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4")


rr_2 <- lapply(y[["mesh_attr"]], rr_func)
rr_y <- grid.arrange(grobs=rr_2, ncol=3)
ggsave(path = "domain_pdf_files", filename = "radius_ratio_2.pdf", plot = rr_y, 
  width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4") 

dev.off()



## combine pdf files

setwd("domain_pdf_files/")
pdfs <- list.files(pattern = "\\.pdf$")
qpdf::pdf_combine(input = pdfs, output = "domain_output.pdf")
setwd("..")

