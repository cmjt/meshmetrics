#' construct mesh with different shape and conditions
#' a wrapper function for inla.mesh.2d
#' @param mesh_mat a matrix that contains all the conditions for the mesh construction.


# mesh_val <- c(loc, boundary, max.edge, cutoff)
mesh_val <- c(0, "region", c(7.5,15), 0,
              "locs", "region", c(7.5,15), 0,
              0, "region", c(7.5,15), 10,
              "locs", "region", c(1.5,2), 3,
              "locs", "region", c(0.75,2), 0.5,
              "locs", "region", c(0.75,1), 0.5)
mesh_mat <- matrix(mesh_val, ncol=5, byrow = TRUE)

mesh_sim <- function(mesh_mat){
  mesh_list <- vector(mode = "list", length = nrow(mesh_mat))
  idx <- which(mesh_mat == 0, arr.ind = TRUE)
  i <- idx[, 1]
  j <- idx[, 2]
  
  if (any(i == 1)){
    mesh_mat <- mesh_mat[unique(idx[,1]), ]
    for (x in 1:nrow(mesh_mat)) {
      param <- mesh_mat[x, ]
      mesh_list_1[[x]] <- inla.mesh.2d(boundary = inla.sp2segment(get(param[2])), 
                                       max.edge = c(param[3], param[4]), cutoff = param[5])
    }
    return(mesh_list_1)
  }
  
  
  if (any(i == 1)) {
    mesh_mat <- mesh_mat[-unique(idx[,1]), ]
    for (x in 1:nrow(mesh_mat)) {
      param <- mesh_mat[x, ]
      mesh_list_4[[x]] <- inla.mesh.2d(loc=get(param[1]), boundary = inla.sp2segment(get(param[2])), 
                                       max.edge = c(param[3], param[4]), cutoff = param[5])
    }
    return(mesh_list_4)
  }
  return(c(mesh_list_1, mesh_list_4))
}

mesh_n <- mesh_sim(mesh_mat)