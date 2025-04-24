require(ks)
require(mclust)
require(spatstat)
require(RandomFields) # 4b


# 1b, oracle

N <- 100
set.seed(887330)
verbose = FALSE
nm.1b = 0

ISE.1b.oracle <- rep(NA, N)

for (n in 1:N){
  cat(n, "  ")
  im.1b.128 <- as.im(lambda1b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.1b <- rpoispp(im.1b.128)
  pp.1b.df <- data.frame(x=pp.1b$x, y=pp.1b$y)
  ise.1b <- function(vh) {
    mat <- invvech(vh)
    dens <- density.ppp(pp.1b, dimyx=c(128,128), varcov=mat)
    dens.dim <- nrow(im.1b.128$v)*ncol(im.1b.128$v)
    ans <- dens.dim * sum((im.1b.128$v - dens$v)^2) / sum(im.1b.128$v)^2
    return(ans)
  }
  res.1b <- optim(vech(Hpi(pp.1b.df)), ise.1b, method = "Nelder-Mead", 
              control = list(trace = as.numeric(verbose), REPORT = 1))
  H.oracle.1b <- invvech(res.1b$par)
  dens.oracle.1b <- density.ppp(pp.1b, dimyx=c(128,128), varcov=H.oracle.1b)
  
  ISE.1b.oracle[n] <- F.computeISE(im.1b.128, dens.oracle.1b)
}



# 1b Thomas, oracle

N <- 100
set.seed(887330)
verbose = FALSE

ISE.1b.thomas.oracle <- rep(NA, N)

for (n in 1:N){
  cat(n, "  ")
  im.1b.128 <- as.im(lambda1b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.1b <- rThomas(kappa=200, scale=0.1, mu=5*im.1b.128*128*128/sum(im.1b.128))
  pp.1b.df <- data.frame(x=pp.1b$x, y=pp.1b$y)
  ise.1b <- function(vh) {
    mat <- invvech(vh)
    dens <- density.ppp(pp.1b, dimyx=c(128,128), varcov=mat)
    dens.dim <- nrow(im.1b.128$v)*ncol(im.1b.128$v)
    ans <- dens.dim * sum((im.1b.128$v - dens$v)^2) / sum(im.1b.128$v)^2
    return(ans)
  }
  res.1b.thomas <- optim(vech(Hpi(pp.1b.df)), ise.1b, method = "Nelder-Mead", 
                  control = list(trace = as.numeric(verbose), REPORT = 1))
  H.oracle.1b.thomas <- invvech(res.1b.thomas$par)
  dens.oracle.1b.thomas <- density.ppp(pp.1b, dimyx=c(128,128), varcov=H.oracle.1b.thomas)
  
  ISE.1b.thomas.oracle[n] <- F.computeISE(im.1b.128, dens.oracle.1b.thomas)
}


# 2b, oracle

N <- 100
set.seed(887330)
verbose = FALSE
nm.2b = 0

ISE.2b.oracle <- rep(NA, N)

for (n in 1:N){
  cat(n, "  ")
  im.2b.128 <- as.im(lambda2b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.2b <- rpoispp(im.2b.128)
  pp.2b.df <- data.frame(x=pp.2b$x, y=pp.2b$y)
  ise.2b <- function(vh) {
    mat <- invvech(vh)
    dens <- density.ppp(pp.2b, dimyx=c(128,128), varcov=mat)
    dens.dim <- nrow(im.2b.128$v)*ncol(im.2b.128$v)
    ans <- dens.dim * sum((im.2b.128$v - dens$v)^2) / sum(im.2b.128$v)^2
    return(ans)
  }
  res.2b <- optim(vech(Hpi(pp.2b.df)), ise.2b, method = "Nelder-Mead", 
                  control = list(trace = as.numeric(verbose), REPORT = 1))
  H.oracle.2b <- invvech(res.2b$par)
  dens.oracle.2b <- density.ppp(pp.2b, dimyx=c(128,128), varcov=H.oracle.2b)
  
  ISE.2b.oracle[n] <- F.computeISE(im.2b.128, dens.oracle.2b)
}




# 2b Thomas, oracle

N <- 100
set.seed(887330)
verbose = FALSE

ISE.2b.thomas.oracle <- rep(NA, N)

for (n in 1:N){
  cat(n, "  ")
  im.2b.128 <- as.im(lambda2b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.2b <- rThomas(kappa=200, scale=0.1, mu=5*im.2b.128*128*128/sum(im.2b.128))
  pp.2b.df <- data.frame(x=pp.2b$x, y=pp.2b$y)
  ise.2b <- function(vh) {
    mat <- invvech(vh)
    dens <- density.ppp(pp.2b, dimyx=c(128,128), varcov=mat)
    dens.dim <- nrow(im.2b.128$v)*ncol(im.2b.128$v)
    ans <- dens.dim * sum((im.2b.128$v - dens$v)^2) / sum(im.2b.128$v)^2
    return(ans)
  }
  res.2b.thomas <- optim(vech(Hpi(pp.2b.df)), ise.2b, method = "Nelder-Mead", 
                         control = list(trace = as.numeric(verbose), REPORT = 1))
  H.oracle.2b.thomas <- invvech(res.2b.thomas$par)
  dens.oracle.2b.thomas <- density.ppp(pp.2b, dimyx=c(128,128), varcov=H.oracle.2b.thomas)
  
  ISE.2b.thomas.oracle[n] <- F.computeISE(im.2b.128, dens.oracle.2b.thomas)
}



# 3b, oracle

N <- 100
set.seed(887330)
verbose = FALSE
nm.3b = 0

ISE.3b.oracle <- rep(NA, N)

for (n in 1:N){
  cat(n, "  ")
  im.3b.128 <- as.im(lambda3b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.3b <- rpoispp(im.3b.128)
  pp.3b.df <- data.frame(x=pp.3b$x, y=pp.3b$y)
  ise.3b <- function(vh) {
    mat <- invvech(vh)
    dens <- density.ppp(pp.3b, dimyx=c(128,128), varcov=mat)
    dens.dim <- nrow(im.3b.128$v)*ncol(im.3b.128$v)
    ans <- dens.dim * sum((im.3b.128$v - dens$v)^2) / sum(im.3b.128$v)^2
    return(ans)
  }
  res.3b <- optim(vech(Hpi(pp.3b.df)), ise.3b, method = "Nelder-Mead", 
                  control = list(trace = as.numeric(verbose), REPORT = 1))
  H.oracle.3b <- invvech(res.3b$par)
  dens.oracle.3b <- density.ppp(pp.3b, dimyx=c(128,128), varcov=H.oracle.3b)
  
  ISE.3b.oracle[n] <- F.computeISE(im.3b.128, dens.oracle.3b)
}



# 4b, oracle

N <- 100
set.seed(887330)
verbose = FALSE
nm.4b = 0

ISE.4b.oracle <- rep(NA, N)

for (n in 1:N){
  cat(n, "  ")
  im.4b.128 <- as.im(lambda4b(my.sim.RF), W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.4b <- rpoispp(im.4b.128)
  pp.4b.df <- data.frame(x=pp.4b$x, y=pp.4b$y)
  ise.4b <- function(vh) {
    mat <- invvech(vh)
    dens <- density.ppp(pp.4b, dimyx=c(128,128), varcov=mat)
    dens.dim <- nrow(im.4b.128$v)*ncol(im.4b.128$v)
    ans <- dens.dim * sum((im.4b.128$v - dens$v)^2) / sum(im.4b.128$v)^2
    return(ans)
  }
  res.4b <- optim(vech(Hpi(pp.4b.df)), ise.4b, method = "Nelder-Mead", 
                  control = list(trace = as.numeric(verbose), REPORT = 1))
  H.oracle.4b <- invvech(res.4b$par)
  dens.oracle.4b <- density.ppp(pp.4b, dimyx=c(128,128), varcov=H.oracle.4b)
  
  ISE.4b.oracle[n] <- F.computeISE(im.4b.128, dens.oracle.4b)
}



# 5b, oracle

N <- 100
set.seed(887330)
verbose = FALSE
nm.5b = 0

ISE.5b.oracle <- rep(NA, N)

for (n in 1:N){
  cat(n, "  ")
  im.5b.128 <- as.im(lambda5b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.5b <- rpoispp(im.5b.128)
  pp.5b.df <- data.frame(x=pp.5b$x, y=pp.5b$y)
  ise.5b <- function(vh) {
    mat <- invvech(vh)
    dens <- density.ppp(pp.5b, dimyx=c(128,128), varcov=mat)
    dens.dim <- nrow(im.5b.128$v)*ncol(im.5b.128$v)
    ans <- dens.dim * sum((im.5b.128$v - dens$v)^2) / sum(im.5b.128$v)^2
    return(ans)
  }
  res.5b <- optim(vech(Hpi(pp.5b.df)), ise.5b, method = "Nelder-Mead", 
                  control = list(trace = as.numeric(verbose), REPORT = 1))
  H.oracle.5b <- invvech(res.5b$par)
  dens.oracle.5b <- density.ppp(pp.5b, dimyx=c(128,128), varcov=H.oracle.5b)
  
  ISE.5b.oracle[n] <- F.computeISE(im.5b.128, dens.oracle.5b)
}



# 6b, oracle

N <- 100
set.seed(887330)
verbose = FALSE
nm.6b = 0

ISE.6b.oracle <- rep(NA, N)

for (n in 1:N){
  cat(n, "  ")
  im.6b.128 <- as.im(lambda6b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.6b <- rpoispp(im.6b.128)
  pp.6b.df <- data.frame(x=pp.6b$x, y=pp.6b$y)
  ise.6b <- function(vh) {
    mat <- invvech(vh)
    dens <- density.ppp(pp.6b, dimyx=c(128,128), varcov=mat)
    dens.dim <- nrow(im.6b.128$v)*ncol(im.6b.128$v)
    ans <- dens.dim * sum((im.6b.128$v - dens$v)^2) / sum(im.6b.128$v)^2
    return(ans)
  }
  res.6b <- optim(vech(Hpi(pp.6b.df)), ise.6b, method = "Nelder-Mead", 
                  control = list(trace = as.numeric(verbose), REPORT = 1))
  H.oracle.6b <- invvech(res.6b$par)
  dens.oracle.6b <- density.ppp(pp.6b, dimyx=c(128,128), varcov=H.oracle.6b)
  
  ISE.6b.oracle[n] <- F.computeISE(im.6b.128, dens.oracle.6b)
}





save(ISE.1b.oracle, ISE.1b.thomas.oracle, ISE.2b.oracle, ISE.2b.thomas.oracle,
     ISE.3b.oracle, ISE.4b.oracle, ISE.5b.oracle, ISE.6b.oracle,
     file = "ISEoracle.RData")



