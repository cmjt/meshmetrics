
## ======================================================================== packages and functions
library(sp)
library(rgeos)
library(INLA)
inla.setOption(mkl=TRUE)
library(spatstat)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork)
library(maptools)
source("functions.r")

## generate 8 different point patterns
source("ppp_sim.R")


## get the point locations for 8 different point patterns
points_df <- lapply(
  unlist(unlist(sim[["points"]], recursive=FALSE), recursive=FALSE), 
  as.data.frame
)

## 12 mesh were made using point locations from 8 different point patterns
## total 96 mesh
source("ppp_mesh.R")

## ======================================================================== fitting the model with 96 mesh

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

## number of simulated points 
for (p in 1:(length(mesh_flatten)/length(mesh_list))) {
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
                      E = inla.stack.data(stk.pp[[q]])$e,
                      control.inla = list(diagonal = 100, strategy = "gaussian", int.strategy = "eb"),
                      verbose = TRUE)
  
  ## transform to user scale
  
  pp.res.est[[q]] <- inla.spde.result(inla = pp.res[[q]], name = "i", spde = spde[[q]], do.transf = TRUE)
  
}


### ======================================================================== Final output CSV file

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

### ------------------------------------------------------------------------->> CSV summary.fixed
fixed <- data.frame(matrix(unlist(list_fixed), ncol = 7, byrow = TRUE))
colnames(fixed) <- names(pp.res[[1]]$summary.fixed)
fixed$Pattern <- rep(paste("Pattern", 1:(length(mesh_flatten)/length(mesh_list)), sep = "_"), 
                     times = length(mesh_list))
fixed$Mesh <- rep(paste("Mesh_band", 1:length(mesh_list), sep = "_"), each = length(mesh_flatten)/length(mesh_list))
write.csv(fixed, file = "loc_summary_fixed.csv")


### ------------------------------------------------------------------------->> CSV marginals.fixed
margfixed <- data.frame(do.call('rbind', marg_fixed))
margfixed$Pattern <- rep(paste("Pattern", 1:(length(mesh_flatten)/length(mesh_list)), sep = "_"), 
                         each = sapply(marg_fixed, nrow)[1], times = length(mesh_list))
margfixed$Mesh <- rep(paste("Mesh_band", 1:length(mesh_list), sep = "_"), 
                      each = sapply(marg_fixed, nrow)[1] * (length(mesh_flatten)/length(mesh_list)))
write.csv(margfixed, file = "loc_marginals_fixed.csv")


### ------------------------------------------------------------------------->> CSV marginals.log.kappa
margkappa <- data.frame(do.call('rbind', marg_kappa))
margkappa$Pattern <- rep(paste("Pattern", 1:(length(mesh_flatten)/length(mesh_list)), sep = "_"), 
                         each = sapply(marg_kappa, nrow)[1], times = length(mesh_list))
margkappa$Mesh <- rep(paste("Mesh_band", 1:length(mesh_list), sep = "_"), 
                      each = sapply(marg_kappa, nrow)[1] * (length(mesh_flatten)/length(mesh_list)))
write.csv(margkappa, file = "loc_marginals_kappa.csv")


### ------------------------------------------------------------------------->> CSV marginals.log.variance.nominal
margvar <- data.frame(do.call('rbind', marg_variance))
margvar$Pattern <- rep(paste("Pattern", 1:(length(mesh_flatten)/length(mesh_list)), sep = "_"), 
                       each = sapply(marg_variance, nrow)[1], times = length(mesh_list))
margvar$Mesh <- rep(paste("Mesh_band", 1:length(mesh_list), sep = "_"), 
                    each = sapply(marg_variance, nrow)[1] * (length(mesh_flatten)/length(mesh_list)))
write.csv(margvar, file = "loc_marginals_variance.csv")


### ------------------------------------------------------------------------->> CSV marginals.log.range.nominal
margrange <- data.frame(do.call('rbind', marg_range))
margrange$Pattern <- rep(paste("Pattern", 1:(length(mesh_flatten)/length(mesh_list)), sep = "_"), 
                         each = sapply(marg_range, nrow)[1], times = length(mesh_list))
margrange$Mesh <- rep(paste("Mesh_band", 1:length(mesh_list), sep = "_"), 
                      each = sapply(marg_range, nrow)[1] * (length(mesh_flatten)/length(mesh_list)))
write.csv(margrange, file = "loc_marginals_range.csv")


### ------------------------------------------------------------------------->> CSV DIC
dic <- data.frame(matrix(unlist(list_dic), ncol = 1, byrow = TRUE))
colnames(dic) <- "DIC"
dic$Pattern <- rep(paste("Pattern", 1:(length(mesh_flatten)/length(mesh_list)), sep = "_"), 
                   times = length(mesh_list))
dic$Mesh <- rep(paste("Mesh_band", 1:length(mesh_list), sep = "_"), 
                each = length(mesh_flatten)/length(mesh_list))
write.csv(dic, file = "loc_DIC.csv")


### ------------------------------------------------------------------------->> CSV WAIC
waic <- data.frame(matrix(unlist(list_waic), ncol = 1, byrow = TRUE))
colnames(waic) <- "WAIC"
waic$Pattern <- rep(paste("Pattern", 1:(length(mesh_flatten)/length(mesh_list)), sep = "_"), 
                    times = length(mesh_list))
waic$Mesh <- rep(paste("Mesh_band", 1:length(mesh_list), sep = "_"), 
                 each = length(mesh_flatten)/length(mesh_list))
write.csv(waic, file = "loc_WAIC.csv")



## ======================================================================== fitting the model with 18*8 mesh

## 18 mesh were made using the domain area.
source("domain_mesh.R")

## number of simulated points
n <- lapply(points_df, nrow)

mesh_list <- c(mesh_list_1, mesh_list_2)

nv <- c()
dual <- list()
weights <- list()
spde <- list()

tmp <- list()
for (l in 1:length(mesh_list)) {
  tmp[[l]] <- vector("list", length = length(points_df))
}

y.pp <- tmp
e.pp <- tmp
imat <- tmp
lmat <- tmp
A.pp <- tmp
stk <- tmp
res <- tmp
res.est <- tmp

for (k in 1:length(mesh_list)) {
  for (l in 1:length(points_df)){
    
    ## number of vertices
    nv[k] <- mesh_list[[k]]$n
    
    ## the dual mesh
    dual[[k]] <- book.mesh.dual(mesh_list[[k]])
    
    ## compute intersection between each polygon
    weights[[k]] <- sapply(1:length(dual[[k]]),
      function(i) {
        if (gIntersects(dual[[k]][i, ], domainSP))
          return(gArea(gIntersection(dual[[k]][i, ], domainSP)))
        else return(0)
      })
    
    ## set up the spde model
    spde[[k]] <- inla.spde2.pcmatern(mesh_list[[k]],
      prior.range = c(0.05, 0.01),
      prior.sigma = c(1, 0.01))
    
    ## projection matrices
    y.pp[[k]][[l]] <- rep(0:1, c(nv[k], n[[l]]))
    e.pp[[k]][[l]] <- c(weights[[k]], rep(0,n[[l]]))
    imat[[k]][[l]] <- Diagonal(nv[k], rep(1, nv[k]))
    lmat[[k]][[l]] <- inla.spde.make.A(mesh = mesh_list[[k]], 
      loc = coordinates(points_df[[l]]))
    A.pp[[k]][[l]] <- rbind(imat[[k]][[l]], lmat[[k]][[l]])
    
    ## set up the data stack
    stk[[k]][[l]] <- inla.stack(
      data = list(y = y.pp[[k]][[l]], e = e.pp[[k]][[l]]),
      A = list(1, A.pp[[k]][[l]]),
      effects = list(list(b0 = rep(1, nv[k] + n[[l]])), list(i = 1:nv[k])))
    
    ## fitting the model
    res[[k]][[l]] <- inla(y ~ 0 + b0 + f(i, model = spde[[k]]),
      family = 'poisson', 
      data = inla.stack.data(stk[[k]][[l]]),
      control.predictor = list(A = inla.stack.A(stk[[k]][[l]])),
      control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
      E = inla.stack.data(stk[[k]][[l]])$e,
      control.inla = list(diagonal = 100, strategy = "gaussian", int.strategy = "eb"),
      verbose = TRUE)
    
    ## transform to user scale
    res.est[[k]][[l]] <- inla.spde.result(inla = res[[k]][[l]], name = "i", spde = spde[[k]], do.transf = TRUE)
    
  }
}


##========================================================================= Model Results

fixed_list <- tmp
marg_fixed <- tmp
marg_kappa <- tmp
marg_variance <- tmp
marg_range <- tmp
dic_list <- tmp
waic_list <- tmp

for (f in 1:length(mesh_list)){
  for (g in 1:length(points_df)){
    fixed_list[[f]][[g]] <- res[[f]][[g]]$summary.fixed
    marg_fixed[[f]][[g]] <- res[[f]][[g]]$marginals.fixed[[1]]
    marg_kappa[[f]][[g]] <- res.est[[f]][[g]]$marginals.kappa[[1]]
    marg_variance[[f]][[g]] <- res.est[[f]][[g]]$marginals.variance.nominal[[1]]
    marg_range[[f]][[g]] <- res.est[[f]][[g]]$marginals.range.nominal[[1]]
    dic_list[[f]][[g]] <- res[[f]][[g]]$dic$dic
    waic_list[[f]][[g]] <- res[[f]][[g]]$waic$waic
  }
}

### ------------------------------------------------------------------------->> CSV summary.fixed
fixed <- data.frame(matrix(unlist(fixed_list), ncol = 7, byrow = TRUE))
colnames(fixed) <- names(res[[1]][[1]]$summary.fixed)
fixed$True <- rep(beta0, each=4)
fixed$Pattern <- rep(paste0("Pattern", 1:length(points_df)), length(mesh_list))
fixed$Mesh <- rep(paste0("Mesh", 1:length(mesh_list)), each = length(points_df))
write.csv(fixed, file = "domain_summary_fixed.csv")

### ------------------------------------------------------------------------->> CSV marginals.fixed
mf <- lapply(rapply(marg_fixed, enquote, how="unlist"), eval)
margfixed <- data.frame(do.call('rbind', mf))
margfixed$Pattern <- rep(paste0("Pattern", 1:length(points_df)), each = sapply(mf,nrow)[1])
margfixed$Mesh <- rep(paste0("Mesh", 1:length(mesh_list)), each = sapply(mf,nrow)[1]*length(points_df))
write.csv(margfixed, file = "domain_marginal_fixed.csv")

### -------------------------------------------------------------------------->> CSV marginals.log.kappa
mk <- lapply(rapply(marg_kappa, enquote, how="unlist"), eval)
margkappa <- data.frame(do.call('rbind', mk))
margkappa$Pattern <- rep(paste0("Pattern", 1:length(points_df)), each = sapply(mk,nrow)[1])
margkappa$Mesh <- rep(paste0("Mesh", 1:length(mesh_list)), each = sapply(mk,nrow)[1]*length(points_df))
write.csv(margkappa, file = "domain_marginal_kappa.csv")


### --------------------------------------------------------------------------->> CSV marginals.log.variance.nominal
mv <- lapply(rapply(marg_variance, enquote, how="unlist"), eval)
margvar <- data.frame(do.call('rbind', mv))
margvar$Pattern <- rep(paste0("Pattern", 1:length(points_df)), each = sapply(mv,nrow)[1])
margvar$Mesh <- rep(paste0("Mesh", 1:length(mesh_list)), each = sapply(mv,nrow)[1]*length(points_df))
write.csv(margvar, file = "domain_marginal_variance.csv")


### --------------------------------------------------------------------------->> CSV marginals.log.range.nominal
mr <- lapply(rapply(marg_range, enquote, how="unlist"), eval)
margrange <- data.frame(do.call('rbind', mr))
margrange$Pattern <- rep(paste0("Pattern", 1:length(points_df)), each = sapply(mr,nrow)[1])
margrange$Mesh <- rep(paste0("Mesh", 1:length(mesh_list)), each = sapply(mr,nrow)[1]*length(points_df))
write.csv(margrange, file = "domain_marginal_range.csv")


### --------------------------------------------------------------------------->> CSV DIC
dic <- data.frame(matrix(unlist(dic_list), ncol = 1, byrow = TRUE))
colnames(dic) <- "DIC"
dic$Pattern <- rep(paste0("Pattern", 1:length(points_df)), length(mesh_list))
dic$Mesh <- rep(paste0("Mesh", 1:length(mesh_list)), each = length(points_df))
write.csv(dic, file = "domain_DIC.csv")


### --------------------------------------------------------------------------->> CSV WAIC
waic <- data.frame(matrix(unlist(waic_list), ncol = 1, byrow = TRUE))
colnames(waic) <- "WAIC"
waic$Pattern <- rep(paste0("Pattern", 1:length(points_df)), length(mesh_list))
waic$Mesh <- rep(paste0("Mesh", 1:length(mesh_list)), each = length(points_df))
write.csv(waic, file = "domain_WAIC.csv")





