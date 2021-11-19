#' Calculate the distance between two points
#' @param A vector of length 2 specifying (first) location
#' @param B vector of length 2 specifying (second) location
dist <- function(A, B){
    sqrt((A[1] - B[1])^2 + (A[2] - B[2])^2)
}
#' cosine triangle function to calculate angle C (in degrees)
#' given edge lengths a, b, and c
#' @param a triangle edgelength (opposits angle A)
#' @param b triangle edgelength (opposits angle B)
#' @param c triangle edgelength (opposits angle C)
ang <- function(a, b, c){
    cosC <- (a^2 + b^2 -c^2)/(2*a*b)
    angC <- acos(cosC)
    return(angC*180/pi)
}

#' Calculate cicumcircle radius given vertex A, B, C
#' locations of a triangle
#' @param A vector of length 2 specifying vertex location "A"
#' @param B vector of length 2 specifying vertex location "B"
#' @param C vector of length 2 specifying vertex location "C"
circum_R <- function(A, B, C){
    a <- dist(B, C)
    b <- dist(A, C)
    c <- dist(B, A)
    abc <- a*b*c
    d1 <- a + b + c
    d2 <- b + c - a
    d3 <- c + a - b
    d4 <- a + b - c
    return((abc)/(sqrt(d1 * d2 * d3 * d4)))
}
#' Calculate cicumcircle centroid given vertex A, B, C
#' locations of a triangle
#' @param A vector of length 2 specifying vertex location "A"
#' @param B vector of length 2 specifying vertex location "B"
#' @param C vector of length 2 specifying vertex location "C"
circum_O <- function(A, B, C){
    a <- dist(B, C)
    b <- dist(A, C)
    c <- dist(B, A)
    angA <- ang(b, c, a)*pi/180
    angB <- ang(a, c, b)*pi/180
    angC <- ang(b, a, c)*pi/180
    sumsins <- sin(2*angA) + sin(2*angB) + sin(2*angC)
    xo <- (A[1]*sin(2*angA) + B[1]*sin(2*angB) + C[1]*sin(2*angC))/sumsins
    yo <- (A[2]*sin(2*angA) + B[2]*sin(2*angB) + C[2]*sin(2*angC))/sumsins
    return(c(xo, yo))
}
#' Calculate incumcircle radius given vertex A, B, C
#' locations of a triangle
#' @param A vector of length 2 specifying vertex location "A"
#' @param B vector of length 2 specifying vertex location "B"
#' @param C vector of length 2 specifying vertex location "C"
incircle_O <- function(A, B, C){
    a <- dist(B, C)
    b <- dist(A, C)
    c <- dist(B, A)
    abc <- a + b + c
    xc <- (a*A[1] + b*B[1] + c* C[1])/abc
    yc <- (a*A[2] + b*B[2] + c* C[2])/abc
    return(c(xc,yc))
}
#' Calculate incumcircle centroid given vertex A, B, C
#' locations of a triangle
#' @param A vector of length 2 specifying vertex location "A"
#' @param B vector of length 2 specifying vertex location "B"
#' @param C vector of length 2 specifying vertex location "C"
incircle_r <- function(A, B, C){
    a <- dist(B, C)
    b <- dist(A, C)
    c <- dist(B, A)
    s <- (a + b + c)/2
    sqrt(((s - a)*(s - b)*(s - c))/s)
}
#' Extract a dataframe of mesh triangle segments (start and end locations)
#' @param mesh an \code{inla.mesh.2d()} object
segs <- function(mesh){
    df <- rbind(data.frame(a = mesh$loc[mesh$graph$tv[, 1], c(1, 2)],
                           b = mesh$loc[mesh$graph$tv[, 2], c(1, 2)]),
                data.frame(a = mesh$loc[mesh$graph$tv[, 2], c(1, 2)],
                           b = mesh$loc[mesh$graph$tv[, 3], c(1, 2)]),
                data.frame(a = mesh$loc[mesh$graph$tv[, 1], c(1, 2)],
                           b = mesh$loc[mesh$graph$tv[, 3], c(1, 2)]))
    colnames(df) <- c("x", "y", "xend", "yend")
    df$length <- raster::pointDistance(df[,1:2], df[,3:4], FALSE, allpairs = FALSE)
    return(df)
}
#' Calculate all interior mesh triangle angles
#' @param  mesh an \code{inla.mesh.2d()} object
#' @param s a dataframe of mesh triangle segments as returned by
#' \code(segs())
mesh_ang <- function(mesh, s){
    tv <- mesh$graph$tv
    angs <- matrix(numeric(3*nrow(tv)), ncol = 3)
    for(i in 1:nrow(tv)){
        ## the three verts of one triangle
        vs <- data.frame(x = rep(mesh$loc[tv[i,],1], each = 2),
                         y = rep(mesh$loc[tv[i,],2], each = 2))
        vs$xend <- vs$x[c(3, 5, 1, 5, 1, 3)];vs$yend <- vs$y[c(3, 5, 1, 5, 1, 3)]
        vs <- vs[c(1,2,4),]
        dists <- unique(plyr::match_df(s, vs, on = NULL))$length
        angs[i, 1] <- ang(dists[1], dists[2], dists[3])
        angs[i, 2] <- ang(dists[2], dists[3], dists[1])
        angs[i, 3] <- ang(dists[3], dists[1], dists[2])
    }
    return(angs)
}
#' Calculate minimum edge length for each triangle in a triangulation
#' @param mesh an \code{inla.mesh.2d()} object
lmin <- function(mesh){
    tv <- mesh$graph$tv
    lmin <- numeric(nrow(tv))
    for (i in 1:nrow(tv)){
        A <- mesh$loc[tv[i,1], 1:2]
        B <- mesh$loc[tv[i,2], 1:2]
        C <- mesh$loc[tv[i,3], 1:2]
        a <- dist(B, C)
        b <- dist(A, C)
        c <- dist(A, B)
        lmin[i] <- min(a, b, c)
    }
    return(lmin)
}
#' Calculate a number of different mesh attributes
#' @param mesh an \code{inla.mesh.2d()} object
#' @return a list of length three: \code{edges} a
#' dataframe of mesh segments as returned by \code{segs};
#' \code{triangles} a dataframe with, for each triangle
#' in the triangulation, the circumcircle and incircle radii
#' (\code{circumcircle_R}, \code{incircle_r}), the assocoated centroid
#' locations (\code{c_Ox, cOy} and \code{i_Ox, iOy}), and
#' the radius-edge ratio \code{re} and radius ratio \code{rr};
#' the third element \code{angles} is a 3 x nvert dataframe of all
#' triangle interior angles

get_triag_attributes <- function(mesh){
    verts <- segs(mesh = mesh)
    angles <- mesh_ang(mesh = mesh, s = verts)
    tv <- mesh$graph$tv
    c_R <- i_R <- numeric(nrow(tv))
    c_O <- i_O <- matrix(rep(0, 2*nrow(tv)), ncol = 2)
    for (i in 1:nrow(tv)){
        A <- mesh$loc[tv[i,1], 1:2]
        B <- mesh$loc[tv[i,2], 1:2]
        C <- mesh$loc[tv[i,3], 1:2]
        c_R[i] <- circum_R(A, B , C)
        i_R[i] <- incircle_r(A, B , C)
        c_O[i, ] <- circum_O(A, B, C)
        i_O[i, ] <- incircle_O(A, B, C)
    }
    mn <- lmin(mesh)
    df <- data.frame(incircle_r = i_R, circumcircle_R = c_R,
                     c_Ox = c_O[, 1],  c_Oy = c_O[, 2],
                     i_Ox = i_O[, 1],  i_Oy = i_O[, 2],
                     re = c_R/mn, rr = i_R/c_R)
    return(list(edges = verts, triangles = df, angles = angles))
}
#' Wrapper function to plot mesh as displayed in ms (using ggplot2)
#' @param mesh an \code{inla.mesh.2d()} object
#' @param xy xy coordinates of points
#' @param domain spatial polygon of region/domain
plt_mesh <- function(mesh, xy, domain){
    points <- data.frame(x = mesh$loc[, 1], y = mesh$loc[, 2])
    tmp <-  ggplot2::ggplot(points, ggplot2::aes(x,y)) +
        ggforce::geom_delaunay_tile(alpha = 0.3, colour = 'black',fill = "transparent") +
        ggplot2::theme_void() +
        ggplot2::geom_sf(data = domain, inherit.aes = FALSE, fill = NA, color = "black", size = 1.5) +
        ggplot2::geom_point(data = xy, ggplot2::aes(x, y), shape = 18, size = 2)
    return(tmp)
}
