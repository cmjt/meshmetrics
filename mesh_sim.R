#' construct mesh with different shape and conditions
#' a wrapper function for inla.mesh.2d
#' @param mesh_mat a matrix that contains all the conditions for the mesh construction.

mesh_sim <- function(mesh_mat){
  mesh_list <- list()
  for (x in 1:nrow(mesh_mat)){
    param <- mesh_mat[x, ]
    if (!sum(param == 0) > 0) {
      mesh_list_1 <- inla.mesh.2d(loc=get(param[1]), boundary = inla.sp2segment(get(param[2])),
                                  max.edge = c(param[3], param[4]), cutoff = param[5])
      mesh_list <- append(mesh_list, list(mesh_list_1))
    } 
    else if (sum(param == 0) > 0) {
      
      loc <- param[1]
      boundary <- param[2]
      inner <- param[3]
      outer <- param[4]
      cutoff <- param[5]
      
      if (sum(param == 0) == 1) {
        if (loc == 0) {
          mesh_list_2 <- inla.mesh.2d(boundary = inla.sp2segment(get(boundary)),
                                      max.edge = c(inner, outer), cutoff = cutoff)
          mesh_list <- append(mesh_list, list(mesh_list_2))
        }
        if (boundary == 0) {
          mesh_list_3 <- inla.mesh.2d(loc = get(loc), max.edge = c(inner, outer), cutoff = cutoff)
          mesh_list <- append(mesh_list, list(mesh_list_3))
        }
        if (outer == 0) { 
          outer <- NULL
          mesh_list_4 <- inla.mesh.2d(loc=get(loc), boundary = inla.sp2segment(get(boundary)),
                                      max.edge = c(inner, outer), cutoff = cutoff)
          mesh_list <- append(mesh_list, list(mesh_list_4))
        }
        if (cutoff == 0) {
          cutoff <- NULL
          mesh_list_5 <- inla.mesh.2d(loc=get(loc), boundary = inla.sp2segment(get(boundary)),
                                      max.edge = c(inner, outer), cutoff = cutoff)
          mesh_list <- append(mesh_list, list(mesh_list_5))
        }
      }
      
      if (sum(param == 0) == 2) {
        if (loc == 0) {
          if (outer == 0) { 
            outer <- NULL
          }
          if (cutoff == 0) {
            cutoff <- NULL
          }
          mesh_list_6 <- inla.mesh.2d(boundary = inla.sp2segment(get(boundary)),
                                      max.edge = c(inner, outer), cutoff = cutoff)
          mesh_list <- append(mesh_list, list(mesh_list_6))
        }
        else if (boundary == 0){
          if (outer == 0) { 
            outer <- NULL
          }
          if (cutoff == 0) {
            cutoff <- NULL
          }
          mesh_list_7 <- inla.mesh.2d(loc = get(loc),
                                      max.edge = c(inner, outer), cutoff = cutoff)
          mesh_list <- append(mesh_list, list(mesh_list_7))
        }
        else {
          outer <- NULL
          cutoff <- NULL
          mesh_list_8 <- inla.mesh.2d(loc=get(loc), boundary = inla.sp2segment(get(boundary)),
                                      max.edge = c(inner, outer), cutoff = cutoff)
          mesh_list <- append(mesh_list, list(mesh_list_8))
        }
      }
      
      if (sum(param == 0) > 2) {
        outer <- NULL
        mesh_list_9 <- inla.mesh.2d(boundary = inla.sp2segment(get(boundary)),
                                    max.edge = c(inner, outer))
        mesh_list <- append(mesh_list, list(mesh_list_9))
      }
    } else {
      stop("Reached the limit of function!")
    }
  }
  return(mesh_list)
}








