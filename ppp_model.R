
## ---- packages and functions
library(sp)
library(rgeos)
library(INLA)
library(spatstat)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork)
library(maptools)
source("functions.r")
source("ppp_sim.R") ## generate 8 different point patterns
source("ppp_mesh.R") ## one mesh was made using point locations from 8 different point patterns


## ---- fitting the model, 8 models in total

n <- list()
nv <- list()

spde <- list()
dmesh <- list()
w <- list()

y.pp <- list()
e.pp <- list()
imat <- list()
lmat <- list()
A.pp <- list()

stk.pp <- list()
pp.res <- list()
pp.res.est <- list()

for (m in 1:length(mesh_ls)) {
  
  ## number of simulated points
  
  n[[m]] <- nrow(points_df[[m]]) 
  
  ## number of vertices in the mesh
  
  nv[[m]] <- mesh_ls[[m]]$n
  
  ## define the SPDE model
  
  spde[[m]] <- inla.spde2.pcmatern(mesh = mesh_ls[[m]],
                                   alpha = 2,
                                   prior.range = c(0.1, 0.01), # P(range < 0.1) = 0.01
                                   prior.sigma = c(0.01, 0.01)) # P(sigma > 0.01) = 0.01)
  
  ## mesh weights
  
  dmesh[[m]] <- book.mesh.dual(mesh_ls[[m]])
  domainSP <- SpatialPolygons(list(Polygons(list(Polygon(domain)), '0')))
  
  w[[m]] <- sapply(1:length(dmesh[[m]]), function(i) {
    if (gIntersects(dmesh[[m]][i, ], domainSP)) 
      return(gArea(gIntersection(dmesh[[m]][i, ], domainSP)))
    else 
      return(0) 
  })
  
  ## projector matrix 
  
  y.pp[[m]] <- rep(0:1, c(nv[[m]], n[[m]]))
  e.pp[[m]] <- c(w[[m]], rep(0, n[[m]]))
  imat[[m]] <- Diagonal(nv[[m]], rep(1, nv[[m]]))
  lmat[[m]] <- inla.spde.make.A(mesh = mesh_ls[[m]], loc = coordinates(points_df[[m]]))
  A.pp[[m]] <- rbind(imat[[m]], lmat[[m]])
  
  ## the data stack
  
  stk.pp[[m]] <- inla.stack(
    data = list(y = y.pp[[m]], e = e.pp[[m]]),
    A = list(1, A.pp[[m]]),
    effects = list(list(b0 = rep(1, nv[[m]] + n[[m]])), list(i = 1:nv[[m]])),
    tag = 'pp')
  
  ## fit the model
  
  pp.res[[m]] <- inla(y ~ 0 + b0 + f(i, model = spde[[m]]),
                      family = 'poisson', 
                      data = inla.stack.data(stk.pp[[m]]),
                      control.predictor = list(A = inla.stack.A(stk.pp[[m]])),
                      control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
                      E = inla.stack.data(stk.pp[[m]])$e)
  
  ## transform to user scale
  
  pp.res.est[[m]] <- inla.spde.result(inla = pp.res[[m]], name = "i", spde = spde[[m]], do.transf = TRUE)
  
}


### --------------------------------------------------------------------------->> Final output CSV file

list_fixed <- list()
marg_fixed <- list()
marg_kappa <- list()
marg_variance <- list()
marg_range <- list()
list_dic <- list()
list_waic <- list()


for (f in 1:nrow(mesh_mat)){
  list_fixed[[f]] <- pp.res[[f]]$summary.fixed
  marg_fixed[[f]] <- pp.res[[f]]$marginals.fixed[[1]]
  marg_kappa[[f]] <- pp.res.est[[f]]$marginals.kappa[[1]]
  marg_variance[[f]] <- pp.res.est[[f]]$marginals.variance.nominal[[1]]
  marg_range[[f]] <- pp.res.est[[f]]$marginals.range.nominal[[1]]
  list_dic[[f]] <- pp.res[[f]]$dic$dic
  list_waic[[f]] <- pp.res[[f]]$waic$waic
}

### --------------------------------------------------------------------------->> CSV summary.fixed
fixed <- data.frame(matrix(unlist(list_fixed), ncol = 7, byrow = TRUE))
colnames(fixed) <- names(pp.res[[1]]$summary.fixed)
fixed$Mesh <- paste("Mesh", 1:nrow(mesh_mat), sep = "_")
write.csv(fixed, file = "summary_fixed.csv")


### --------------------------------------------------------------------------->> CSV marginals.fixed
margfixed <- data.frame(do.call('rbind', marg_fixed))
margfixed$Mesh <- rep(paste("Mesh", 1:nrow(mesh_mat), sep = "_"), each = nrow(pp.res[[1]]$marginals.fixed[[1]]))
write.csv(margfixed, file = "marginals_fixed.csv")


### --------------------------------------------------------------------------->> CSV marginals.log.kappa
margkappa <- data.frame(do.call('rbind', marg_kappa))
margkappa$Mesh <- rep(paste("Mesh", 1:nrow(mesh_mat), sep = "_"), each = nrow(pp.res.est[[1]]$marginals.kappa[[1]]))
write.csv(margkappa, file = "marginals_kappa.csv")


### --------------------------------------------------------------------------->> CSV marginals.log.variance.nominal
margvar <- data.frame(do.call('rbind', marg_variance))
margvar$Mesh <- rep(paste("Mesh", 1:nrow(mesh_mat), sep = "_"), each = nrow(pp.res.est[[1]]$marginals.variance.nominal[[1]]))
write.csv(margvar, file = "marginals_variance.csv")


### --------------------------------------------------------------------------->> CSV marginals.log.range.nominal
margrange <- data.frame(do.call('rbind', marg_range))
margrange$Mesh <- rep(paste("Mesh", 1:nrow(mesh_mat), sep = "_"), each = nrow(pp.res.est[[1]]$marginals.range.nominal[[1]]))
write.csv(margrange, file = "marginals_range.csv")


### --------------------------------------------------------------------------->> CSV DIC
dic <- data.frame(matrix(unlist(list_dic), ncol = 1, byrow = TRUE))
colnames(dic) <- "DIC"
dic$Mesh <- paste("Mesh", 1:nrow(mesh_mat), sep = "_")
write.csv(dic, file = "DIC.csv")


### --------------------------------------------------------------------------->> CSV WAIC
waic <- data.frame(matrix(unlist(list_waic), ncol = 1, byrow = TRUE))
colnames(waic) <- "WAIC"
waic$Mesh <- paste("Mesh", 1:nrow(mesh_mat), sep = "_")
write.csv(waic, file = "WAIC.csv")




