
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
## generate 8 different point patterns
source("ppp_sim.R")
## 4 mesh were made using point locations from 8 different point patterns
source("ppp_mesh.R")


## ---- fitting the model, 32 models in total

mesh_flatten <- unlist(mesh_list, recursive = FALSE)

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

for (p in 1:(length(mesh_flatten)/length(mesh_list))) {
  
  ## number of simulated points
  
  n[[p]] <- nrow(points_df[[p]])
}

n <- rep(unlist(n), length(mesh_list))


for (q in 1:length(mesh_flatten)) {
  
  ## number of vertices in the mesh
  
  nv[[q]] <- mesh_flatten[[q]]$n
  
  ## define the SPDE model
  
  spde[[q]] <- inla.spde2.pcmatern(mesh = mesh_flatten[[q]],
                                   alpha = 2,
                                   prior.range = c(0.1, 0.01), # P(range < 0.1) = 0.01
                                   prior.sigma = c(0.01, 0.01)) # P(sigma > 0.01) = 0.01)
  
  ## mesh weights
  
  dmesh[[q]] <- book.mesh.dual(mesh_flatten[[q]])
  domainSP <- SpatialPolygons(list(Polygons(list(Polygon(domain)), '0')))
  
  w[[q]] <- sapply(1:length(dmesh[[q]]), function(i) {
    if (gIntersects(dmesh[[q]][i, ], domainSP)) 
      return(gArea(gIntersection(dmesh[[q]][i, ], domainSP)))
    else 
      return(0) 
  })
  
  ## projector matrix 
  
  y.pp[[q]] <- rep(0:1, c(nv[[q]], n[q]))
  e.pp[[q]] <- c(w[[q]], rep(0, n[q]))
  imat[[q]] <- Diagonal(nv[[q]], rep(1, nv[[q]]))
  lmat[[q]] <- inla.spde.make.A(mesh = mesh_flatten[[q]], 
                                loc = coordinates(rep(points_df, length(mesh_list))[[q]]))
  A.pp[[q]] <- rbind(imat[[q]], lmat[[q]])
  
  ## the data stack
  
  stk.pp[[q]] <- inla.stack(
    data = list(y = y.pp[[q]], e = e.pp[[q]]),
    A = list(1, A.pp[[q]]),
    effects = list(list(b0 = rep(1, nv[[q]] + n[q])), list(i = 1:nv[[q]])),
    tag = 'pp')
  
  ## fit the model
  
  pp.res[[q]] <- inla(y ~ 0 + b0 + f(i, model = spde[[q]]),
                      family = 'poisson', 
                      data = inla.stack.data(stk.pp[[q]]),
                      control.predictor = list(A = inla.stack.A(stk.pp[[q]])),
                      control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
                      E = inla.stack.data(stk.pp[[q]])$e)
  
  ## transform to user scale
  
  pp.res.est[[q]] <- inla.spde.result(inla = pp.res[[q]], name = "i", spde = spde[[q]], do.transf = TRUE)
  
}


### --------------------------------------------------------------------------->> Final output CSV file

list_fixed <- list()
marg_fixed <- list()
marg_kappa <- list()
marg_variance <- list()
marg_range <- list()
list_dic <- list()
list_waic <- list()


for (f in 1:length(mesh_flatten)){
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
fixed$Mesh <- paste("Mesh", 1:length(mesh_list), sep = "_")
write.csv(fixed, file = "summary_fixed.csv")


### --------------------------------------------------------------------------->> CSV marginals.fixed
margfixed <- data.frame(do.call('rbind', marg_fixed))
margfixed$Mesh <- rep(paste("Mesh", 1:length(mesh_list), sep = "_"), each = nrow(pp.res[[1]]$marginals.fixed[[1]]))
write.csv(margfixed, file = "marginals_fixed.csv")


### --------------------------------------------------------------------------->> CSV marginals.log.kappa
margkappa <- data.frame(do.call('rbind', marg_kappa))
margkappa$Mesh <- rep(paste("Mesh", 1:length(mesh_list), sep = "_"), each = nrow(pp.res.est[[1]]$marginals.kappa[[1]]))
write.csv(margkappa, file = "marginals_kappa.csv")


### --------------------------------------------------------------------------->> CSV marginals.log.variance.nominal
margvar <- data.frame(do.call('rbind', marg_variance))
margvar$Mesh <- rep(paste("Mesh", 1:length(mesh_list), sep = "_"), each = nrow(pp.res.est[[1]]$marginals.variance.nominal[[1]]))
write.csv(margvar, file = "marginals_variance.csv")


### --------------------------------------------------------------------------->> CSV marginals.log.range.nominal
margrange <- data.frame(do.call('rbind', marg_range))
margrange$Mesh <- rep(paste("Mesh", 1:length(mesh_list), sep = "_"), each = nrow(pp.res.est[[1]]$marginals.range.nominal[[1]]))
write.csv(margrange, file = "marginals_range.csv")


### --------------------------------------------------------------------------->> CSV DIC
dic <- data.frame(matrix(unlist(list_dic), ncol = 1, byrow = TRUE))
colnames(dic) <- "DIC"
dic$Mesh <- paste("Mesh", 1:length(mesh_list), sep = "_")
write.csv(dic, file = "DIC.csv")


### --------------------------------------------------------------------------->> CSV WAIC
waic <- data.frame(matrix(unlist(list_waic), ncol = 1, byrow = TRUE))
colnames(waic) <- "WAIC"
waic$Mesh <- paste("Mesh", 1:length(mesh_list), sep = "_")
write.csv(waic, file = "WAIC.csv")



