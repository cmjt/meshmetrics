
ppp_sim <- function(nsim, beta0, sigma2x, kappa, win, seed = 1234, plot = TRUE) {
  
  argg <- c(as.list(environment()), list())
  
  temp <- list()
  for (x in 1:2) {
    temp[[x]] <- vector("list", length = 2)
    for (y in 1:2) {
      temp[[x]][[y]] <- vector("list", length = 2)
    }
  }
  
  sim <- temp
  npoints <- temp
  
  for (j in 1:length(beta0)) {
    for (k in 1:length(sigma2x)) {
      for (l in 1:length(kappa)) {
        set.seed(seed)
        sim[[j]][[k]][[l]] <- rLGCP('matern', nsim = nsim, mu = beta0[j], var = sigma2x[k], 
                                    scale = 1 / kappa[l], nu = 1, win = as.owin(domain))
      }
    }
  }
  
  for (a in 1:length(beta0)) {
    for (b in 1:length(sigma2x)) {
      for (c in 1:length(kappa)) {
        npoints[[a]][[b]][[c]] <- sim[[a]][[b]][[c]]$n
      }
    }
  }
  
  return(list(simulation = sim, points = npoints))
}


