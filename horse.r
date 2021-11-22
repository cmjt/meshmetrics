url <- "https://gist.githubusercontent.com/cmjt/9a5e43cce69a3babb129b0a448e65752/raw/8653d11cb6ace667a3f5683982ae82366a3ac5d2/horse.csv"
horse <- read.csv(url)
source("functions.r")
library(sp)
library(INLA)
library(ggplot2)
library(patchwork)
sp <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(horse)), '0')))
bound <- inla.sp2segment(sp)

m1 <- INLA::inla.mesh.2d(boundary = bound,
                         max.edge = 2)

m2 <- INLA::inla.mesh.2d(boundary = bound,
                         max.edge = 2, cutoff = 1.5)
m3 <- INLA::inla.mesh.create(horse, bound = bound)
m4 <- INLA::inla.mesh.2d(boundary = bound,
                         max.edge = c(0.5,2))

m4 <- INLA::inla.mesh.2d(boundary = bound,
                         max.edge = c(0.5,2))

m5 <- INLA::inla.mesh.2d(loc = horse,
                         max.edge = 0.5)

png("horse_mesh_example.png", width = 1000, height = 500)
par(mfrow = c(1,2), mar = c(5,0,0,0))
plot(m3, asp = 1, main = "")
plot(sp, add = TRUE, lwd = 3)
mtext(1, line = -3, text = "A")
plot(m1, asp = 1, main = "")
plot(sp, add = TRUE, lwd = 3)
mtext(1, line = -3, text = "B")
dev.off()

## *********** ##
## mesh 1 vs mesh 3

m1_tris <- get_triag_attributes(m1)
m3_tris <- get_triag_attributes(m3)

## mesh attribute plots

rr3 <- ggplot(m3_tris$sf, aes(fill = m3_tris$triangles$rr)) +
    geom_sf(colour = NA) + theme_void() +
    scale_fill_continuous(name = "Radius Ratio", limits = c(0, 0.5))


rr1 <- ggplot(m1_tris$sf, aes(fill = m1_tris$triangles$rr)) +
    geom_sf(colour = NA) + theme_void() +
    scale_fill_continuous(name = "Radius Ratio", limits = c(0, 0.5))

(rr3 / rr1) + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect")


re3 <- ggplot(m3_tris$sf, aes(fill = m3_tris$triangles$re)) +
    geom_sf(colour = NA) + theme_void() +
    scale_fill_continuous(name = "Radius-edge ratio" , limits = c(1/sqrt(3), NA))


re1 <- ggplot(m1_tris$sf, aes(fill = m1_tris$triangles$re)) +
    geom_sf(colour = NA) + theme_void() +
    scale_fill_continuous(name = "Radius-edge ratio", limits = c(1/sqrt(3), NA))

(re3 / re1) + plot_annotation(tag_levels = 'A') 

## circumcircle plots
png("horse_circums.png", width = 1000, height = 500)
par(xpd = TRUE, mfrow = c(1, 2), mar = c(5,0,0,0))
## mesh 3
plot(m3_tris$triangles$c_Ox, m3_tris$triangles$c_Oy, ,axes = FALSE, ylab = "",
     xlab = "", pch = 20, cex = 0.5, asp = 1)
mtext(1, line = -3, text = "A")
for(i in 1:nrow(m3_tris$triangles)){
    plotrix::draw.circle(m3_tris$triangles$c_Ox[i], m3_tris$triangles$c_Oy[i],
                         m3_tris$triangles$circumcircle_R[i])
}
for(i in 1:nrow(m3_tris$triangles)){
    plotrix::draw.circle(m3_tris$triangles$i_Ox[i], m3_tris$triangles$i_Oy[i],
                         m3_tris$triangles$incircle_r[i], border = "red")
}
## mesh 1
plot(m1_tris$triangles$c_Ox, m1_tris$triangles$c_Oy, ,axes = FALSE, ylab = "",
     xlab = "", pch = 20, cex = 0.5, asp = 1)
mtext(1, line = -3, text = "B")
for(i in 1:nrow(m1_tris$triangles)){
    plotrix::draw.circle(m1_tris$triangles$c_Ox[i], m1_tris$triangles$c_Oy[i],
                         m1_tris$triangles$circumcircle_R[i])
}
for(i in 1:nrow(m1_tris$triangles)){
    plotrix::draw.circle(m1_tris$triangles$i_Ox[i], m1_tris$triangles$i_Oy[i],
                         m1_tris$triangles$incircle_r[i], border = "red")
}
dev.off()

## density plots
all <- data.frame(radius_edge = c(m1_tris$triangles$re, m3_tris$triangles$re),
                  radius_ratio = c(m1_tris$triangles$rr, m3_tris$triangles$rr),
                  mesh = c(rep("B", nrow(m1_tris$triangles)), rep("A", nrow(m3_tris$triangles))))
all_es <- data.frame(edge_length = c(m1_tris$edges$length, m3_tris$edges$length),
                     angles = c(c(m1_tris$angles), c(m3_tris$angles)),
                     mesh = c(rep("B", nrow(m1_tris$edges)), rep("A", nrow(m3_tris$edges))))

at <- seq(1,12, by = 2)*(1/sqrt(3))
L <- parse( text = paste("frac(", seq(1,12, by = 2), ", sqrt(3))", sep = ""))
p1 <- ggplot(all, aes(y = radius_edge, x = mesh)) +
    theme_classic() + geom_violin(fill = "grey") +
    geom_boxplot(width = 0.07) +
    scale_y_continuous(breaks = at, labels  = L) + 
    coord_flip() +
    geom_hline(yintercept = 1/sqrt(3)) + 
    theme(legend.position = "none") +
    ggtitle("Radius-edge ratio") + ylab("") + xlab("")
at <- c(0, 1/4, 1/2)
L <- parse( text = paste("frac( 1,", c(4,2),")", sep = ""))

p2 <- ggplot(all, aes(y = radius_ratio, x = mesh)) +
    theme_classic() + geom_violin(fill = "grey") +
    geom_boxplot(width = 0.07) +
    scale_y_continuous(breaks = at, labels  = c(0, L)) + 
    coord_flip() +
    geom_hline(yintercept = 1/2) + 
    theme(legend.position = "none") +
    ggtitle("Radius ratio") + ylab("") + xlab("")

p1 + p2
## edges

p3 <- ggplot(all_es, aes(y = edge_length, x = mesh)) +
    theme_classic() + geom_violin(fill = "grey") +
    geom_boxplot(width = 0.07) +
    coord_flip() +
    theme(legend.position = "none") +
    ggtitle("Edgelengths") + ylab("") + xlab("")


at <- seq(0,150, by = 30)
L <- parse(text = paste(seq(0,150, by = 30), "*degree ", sep = ""))
p4 <- ggplot(all_es, aes(y = angles, x = mesh)) +
    theme_classic() + geom_violin(fill = "grey") +
    geom_boxplot(width = 0.07) +
    scale_y_continuous(breaks = at, labels  = L) +
    geom_hline(yintercept = 60) + 
    coord_flip() +
    theme(legend.position = "none") +
    ggtitle("Interior angles") + ylab("") + xlab("")

(p1 + p2)/(p3 + p4)

ggsave("horse_mesh_attributes.png", width = 7, height = 7)
