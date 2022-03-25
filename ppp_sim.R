
ppp_sim <- function(nsim, beta0, sigma2x, prac_range, win) {
  argg <- c(as.list(environment()), list())
  xy <- points <- rf <- dat <- lam_attr <- plts <- list()
  for(i in 1:length(beta0)){
    xy[[i]] <- points[[i]] <- rf[[i]] <- dat[[i]] <- lam_attr[[i]] <- plts[[i]] <- list()
    for(j in 1:length(sigma2x)){
      xy[[i]][[j]] <- points[[i]][[j]] <- rf[[i]][[j]] <- dat[[i]][[j]] <- lam_attr[[i]][[j]] <- plts[[i]][[j]] <-  list()
      for(k in 1:length(prac_range)){
        xy[[i]][[j]][[k]] <- rLGCP('matern', nsim = nsim, mu = beta0[i], var = sigma2x[j],
                                   scale = prac_range[k] / sqrt(8), nu = 1, win = as.owin(win))
        points[[i]][[j]][[k]] <- data.frame(x = xy[[i]][[j]][[k]]$x, y = xy[[i]][[j]][[k]]$y)
        lam_attr[[i]][[j]][[k]] <- summary(as.vector(log(attr(xy[[i]][[j]][[k]], 'Lambda')$v)))
        rf[[i]][[j]][[k]] <- apply(log(attr(xy[[i]][[j]][[k]], 'Lambda')$v), 2, t) ## field flips
        ex <- expand.grid(seq(0,2,length.out = 300), seq(0,2,length.out = 300))
        dat[[i]][[j]][[k]] <- data.frame(x = ex$Var1, y = ex$Var2, rf = as.vector(t(rf[[i]][[j]][[k]])))
        plts[[i]][[j]][[k]] <- 
          ggplot(dat[[i]][[j]][[k]], aes(x = x, y = y, fill = rf)) +
          geom_tile() +
          labs(fill = "Expected number \n of points") +
          geom_point(data = points[[i]][[j]][[k]], aes(x = x, y = y), inherit.aes = FALSE,
                     shape = 18, size = 2, alpha = 0.5) +
          scale_fill_viridis_b() +
          theme_void()  +
          geom_path(data = domain, aes(x, y), size = 3, inherit.aes = FALSE) +
          ggtitle(bquote(paste("Parameter values: ", beta[0] , " = ", .(beta0[i]), "; ",
                               sigma[z]^2 , " = ",  .(sigma2x[j]), "; ",
                               rho , " = ",  .(prac_range[k]), ".", sep = "")))
      }
    }
  }
  return(list(sim = xy, points = points, dat = dat, lam = lam_attr, plot = plts))
}


## ---- data_and_domain

domain <- data.frame(x = c(0, 2, 2, 0, 0),y = c(0, 0, 2, 2, 0))
domainSP <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(domain)), '0')))
npix <- 300
spatstat.options(npixel = npix)
beta0 <- c(5, 2)
sigma2x <- c(5, 0.1) 
prac_range <- c(5, 0.2)
sim <- ppp_sim(1, beta0, sigma2x, prac_range, domainSP)

points_plts <- (sim[["plot"]][[1]][[1]][[1]] | sim[["plot"]][[1]][[1]][[2]])/
  (sim[["plot"]][[1]][[2]][[1]] | sim[["plot"]][[1]][[2]][[2]])/
  (sim[["plot"]][[2]][[1]][[1]] | sim[["plot"]][[2]][[1]][[2]])/
  (sim[["plot"]][[2]][[2]][[1]] | sim[["plot"]][[2]][[2]][[2]])


## create a new file to save all the plots

dir.create("ppp_pdf_files")

ggsave(path = "ppp_pdf_files", filename = "point_patterns.pdf", plot = points_plts, 
       width = 8, height = 12, bg = "white", colormodel = "cmyk", paper = "A4")



