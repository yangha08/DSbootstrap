
require(spatstat)
require(ks)
require(mvtnorm)

F.computeISE <- function(true.dens, est.dens){
  dens.dim <- nrow(true.dens$v)*ncol(true.dens$v)
  #ans <- (sum((true.dens$v - est.dens$v)^2) / dens.dim) / (sum(true.dens$v)/dens.dim)^2
  ans <- dens.dim * sum((true.dens$v - est.dens$v)^2) / sum(true.dens$v)^2
  return(ans)
}

pilotg <- function (x, binned, deriv.order = 0, verbose = FALSE) 
{
  n <- nrow(x)
  d <- ncol(x)
  r <- deriv.order
  S <- var(x)
  binned <- default.bflag(d = d, n = n)
  
  g6 <- gsamse(S, n=n, modr=6)
  psihat6.star <- kfe(x = x, G = g6^2 * diag(d), 
                      deriv.order = 6, deriv.vec = TRUE, binned = binned, 
                      add.index = FALSE, verbose = verbose)
  g <- gsamse(S, n = n, modr = 4, nstage = 2, 
              psihat = psihat6.star)
  G <- g^2 * diag(d)
  return(G)
}


# 1b, empirical bootstrap

N <- 100
set.seed(887330)
verbose = FALSE

ISE.1b.empboot <- rep(NA, N)

im.1b.128 <- as.im(lambda1b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
dens.dim <- nrow(im.1b.128$v)*ncol(im.1b.128$v)

for(j in 1:N) {
  cat(j, '\t')
  pp.1b <- rpoispp(im.1b.128)
  pp.1b.df <- data.frame(x=pp.1b$x, y=pp.1b$y)
  G <- pilotg(pp.1b.df)
  dens.pilot <- density.ppp(pp.1b, dimyx=c(128,128), varcov=G)
  m <- sum(dens.pilot/128^2)
  bmise.1b <- function(vh) {
    mat <- invvech(vh)
    ans <- 0
    for(i in 1:20) {
      N.star <- rpois(1, m)
      boot.index <- sample(1:nrow(pp.1b.df), N.star, replace = T)
      pp.1b.SB.df <- pp.1b.df[boot.index,] + as.data.frame(rmvnorm.mixt(n = N.star, Sigmas = G))
      pp.1b.SB <- as.ppp(pp.1b.SB.df, owin(c(0,1),c(0,1)))
      dens <- density.ppp(pp.1b.SB, dimyx=c(128,128), varcov=mat)
      ans <- ans + dens.dim * sum((dens.pilot$v - dens$v)^2) / sum(dens.pilot$v)^2
    }
    ans <- ans/N
    return(ans)
  }
  res.1b <- optim(vech(Hpi(pp.1b.df)), bmise.1b, method = "Nelder-Mead", 
                  control = list(trace = as.numeric(verbose), REPORT = 1))
  H.empboot.1b <- invvech(res.1b$par)
  dens.empboot.1b <- density.ppp(pp.1b, dimyx=c(128,128), varcov=H.empboot.1b)
  
  ISE.1b.empboot[j] <- F.computeISE(im.1b.128, dens.empboot.1b)
}




# 1b Thomas, empirical bootstrap

N <- 100
set.seed(887330)
verbose = FALSE

ISE.1b.thomas.empboot <- rep(NA, N)

im.1b.128 <- as.im(lambda1b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
dens.dim <- nrow(im.1b.128$v)*ncol(im.1b.128$v)

for(j in 1:N) {
  cat(j, '\t')
  pp.1b <- rThomas(kappa=200, scale=0.1, mu=5*im.1b.128*128*128/sum(im.1b.128))
  pp.1b.df <- data.frame(x=pp.1b$x, y=pp.1b$y)
  G <- pilotg(pp.1b.df)
  dens.pilot <- density.ppp(pp.1b, dimyx=c(128,128), varcov=G)
  m <- sum(dens.pilot/128^2)
  bmise.1b <- function(vh) {
    mat <- invvech(vh)
    ans <- 0
    for(i in 1:20) {
      N.star <- rpois(1, m)
      boot.index <- sample(1:nrow(pp.1b.df), N.star, replace = T)
      pp.1b.SB.df <- pp.1b.df[boot.index,] + as.data.frame(rmvnorm.mixt(n = N.star, Sigmas = G))
      pp.1b.SB <- as.ppp(pp.1b.SB.df, owin(c(0,1),c(0,1)))
      dens <- density.ppp(pp.1b.SB, dimyx=c(128,128), varcov=mat)
      ans <- ans + dens.dim * sum((dens.pilot$v - dens$v)^2) / sum(dens.pilot$v)^2
    }
    ans <- ans/N
    return(ans)
  }
  res.1b <- optim(vech(Hpi(pp.1b.df)), bmise.1b, method = "Nelder-Mead", 
                  control = list(trace = as.numeric(verbose), REPORT = 1))
  H.empboot.1b <- invvech(res.1b$par)
  dens.empboot.1b <- density.ppp(pp.1b, dimyx=c(128,128), varcov=H.empboot.1b)
  
  ISE.1b.thomas.empboot[j] <- F.computeISE(im.1b.128, dens.empboot.1b)
}




# 2b, empirical bootstrap

N <- 100
set.seed(887330)
verbose = FALSE

ISE.2b.empboot <- rep(NA, N)

im.2b.128 <- as.im(lambda2b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
dens.dim <- nrow(im.2b.128$v)*ncol(im.2b.128$v)

for(j in 1:N) {
  cat(j, '\t')
  pp.2b <- rpoispp(im.2b.128)
  pp.2b.df <- data.frame(x=pp.2b$x, y=pp.2b$y)
  G <- pilotg(pp.2b.df)
  dens.pilot <- density.ppp(pp.2b, dimyx=c(128,128), varcov=G)
  m <- sum(dens.pilot/128^2)
  bmise.2b <- function(vh) {
    mat <- invvech(vh)
    ans <- 0
    for(i in 1:20) {
      N.star <- rpois(1, m)
      boot.index <- sample(1:nrow(pp.2b.df), N.star, replace = T)
      pp.2b.SB.df <- pp.2b.df[boot.index,] + as.data.frame(rmvnorm.mixt(n = N.star, Sigmas = G))
      pp.2b.SB <- as.ppp(pp.2b.SB.df, owin(c(0,1),c(0,1)))
      dens <- density.ppp(pp.2b.SB, dimyx=c(128,128), varcov=mat)
      ans <- ans + dens.dim * sum((dens.pilot$v - dens$v)^2) / sum(dens.pilot$v)^2
    }
    ans <- ans/N
    return(ans)
  }
  res.2b <- optim(vech(Hpi(pp.2b.df)), bmise.2b, method = "Nelder-Mead", 
                  control = list(trace = as.numeric(verbose), REPORT = 1))
  H.empboot.2b <- invvech(res.2b$par)
  dens.empboot.2b <- density.ppp(pp.2b, dimyx=c(128,128), varcov=H.empboot.2b)
  
  ISE.2b.empboot[j] <- F.computeISE(im.2b.128, dens.empboot.2b)
}




# 2b thomas, empirical bootstrap

N <- 100
set.seed(887330)
verbose = FALSE

ISE.2b.thomas.empboot <- rep(NA, N)

im.2b.128 <- as.im(lambda2b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
dens.dim <- nrow(im.2b.128$v)*ncol(im.2b.128$v)

for(j in 1:N) {
  cat(j, '\t')
  pp.2b <- rThomas(kappa=200, scale=0.1, mu=5*im.2b.128*128*128/sum(im.2b.128))
  pp.2b.df <- data.frame(x=pp.2b$x, y=pp.2b$y)
  G <- pilotg(pp.2b.df)
  dens.pilot <- density.ppp(pp.2b, dimyx=c(128,128), varcov=G)
  m <- sum(dens.pilot/128^2)
  bmise.2b <- function(vh) {
    mat <- invvech(vh)
    ans <- 0
    for(i in 1:20) {
      N.star <- rpois(1, m)
      boot.index <- sample(1:nrow(pp.2b.df), N.star, replace = T)
      pp.2b.SB.df <- pp.2b.df[boot.index,] + as.data.frame(rmvnorm.mixt(n = N.star, Sigmas = G))
      pp.2b.SB <- as.ppp(pp.2b.SB.df, owin(c(0,1),c(0,1)))
      dens <- density.ppp(pp.2b.SB, dimyx=c(128,128), varcov=mat)
      ans <- ans + dens.dim * sum((dens.pilot$v - dens$v)^2) / sum(dens.pilot$v)^2
    }
    ans <- ans/N
    return(ans)
  }
  res.2b <- optim(vech(Hpi(pp.2b.df)), bmise.2b, method = "Nelder-Mead", 
                  control = list(trace = as.numeric(verbose), REPORT = 1))
  H.empboot.2b <- invvech(res.2b$par)
  dens.empboot.2b <- density.ppp(pp.2b, dimyx=c(128,128), varcov=H.empboot.2b)
  
  ISE.2b.thomas.empboot[j] <- F.computeISE(im.2b.128, dens.empboot.2b)
}



# 3b, empirical bootstrap

N <- 100
set.seed(887330)
verbose = FALSE

ISE.3b.empboot <- rep(NA, N)

im.3b.128 <- as.im(lambda3b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
dens.dim <- nrow(im.3b.128$v)*ncol(im.3b.128$v)

for(j in 1:N) {
  cat(j, '\t')
  pp.3b <- rpoispp(im.3b.128)
  pp.3b.df <- data.frame(x=pp.3b$x, y=pp.3b$y)
  G <- pilotg(pp.3b.df)
  dens.pilot <- density.ppp(pp.3b, dimyx=c(128,128), varcov=G)
  m <- sum(dens.pilot/128^2)
  bmise.3b <- function(vh) {
    mat <- invvech(vh)
    ans <- 0
    for(i in 1:20) {
      N.star <- rpois(1, m)
      boot.index <- sample(1:nrow(pp.3b.df), N.star, replace = T)
      pp.3b.SB.df <- pp.3b.df[boot.index,] + as.data.frame(rmvnorm.mixt(n = N.star, Sigmas = G))
      pp.3b.SB <- as.ppp(pp.3b.SB.df, owin(c(0,1),c(0,1)))
      dens <- density.ppp(pp.3b.SB, dimyx=c(128,128), varcov=mat)
      ans <- ans + dens.dim * sum((dens.pilot$v - dens$v)^2) / sum(dens.pilot$v)^2
    }
    ans <- ans/N
    return(ans)
  }
  res.3b <- optim(vech(Hpi(pp.3b.df)), bmise.3b, method = "Nelder-Mead", 
                  control = list(trace = as.numeric(verbose), REPORT = 1))
  H.empboot.3b <- invvech(res.3b$par)
  dens.empboot.3b <- density.ppp(pp.3b, dimyx=c(128,128), varcov=H.empboot.3b)
  
  ISE.3b.empboot[j] <- F.computeISE(im.3b.128, dens.empboot.3b)
}


# 4b, empirical bootstrap

N <- 100
set.seed(887330)
verbose = FALSE

ISE.4b.empboot <- rep(NA, N)

im.4b.128 <- as.im(lambda4b(my.sim.RF), W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
dens.dim <- nrow(im.4b.128$v)*ncol(im.4b.128$v)

for(j in 1:N) {
  cat(j, '\t')
  pp.4b <- rpoispp(im.4b.128)
  pp.4b.df <- data.frame(x=pp.4b$x, y=pp.4b$y)
  G <- pilotg(pp.4b.df)
  dens.pilot <- density.ppp(pp.4b, dimyx=c(128,128), varcov=G)
  m <- sum(dens.pilot/128^2)
  bmise.4b <- function(vh) {
    mat <- invvech(vh)
    ans <- 0
    for(i in 1:20) {
      N.star <- rpois(1, m)
      boot.index <- sample(1:nrow(pp.4b.df), N.star, replace = T)
      pp.4b.SB.df <- pp.4b.df[boot.index,] + as.data.frame(rmvnorm.mixt(n = N.star, Sigmas = G))
      pp.4b.SB <- as.ppp(pp.4b.SB.df, owin(c(0,1),c(0,1)))
      dens <- density.ppp(pp.4b.SB, dimyx=c(128,128), varcov=mat)
      ans <- ans + dens.dim * sum((dens.pilot$v - dens$v)^2) / sum(dens.pilot$v)^2
    }
    ans <- ans/N
    return(ans)
  }
  res.4b <- optim(vech(Hpi(pp.4b.df)), bmise.4b, method = "Nelder-Mead", 
                  control = list(trace = as.numeric(verbose), REPORT = 1))
  H.empboot.4b <- invvech(res.4b$par)
  dens.empboot.4b <- density.ppp(pp.4b, dimyx=c(128,128), varcov=H.empboot.4b)
  
  ISE.4b.empboot[j] <- F.computeISE(im.4b.128, dens.empboot.4b)
}





# 5b, empirical bootstrap

N <- 100
set.seed(887330)
verbose = FALSE

ISE.5b.empboot <- rep(NA, N)

im.5b.128 <- as.im(lambda5b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
dens.dim <- nrow(im.5b.128$v)*ncol(im.5b.128$v)

for(j in 1:N) {
  cat(j, '\t')
  pp.5b <- rpoispp(im.5b.128)
  pp.5b.df <- data.frame(x=pp.5b$x, y=pp.5b$y)
  G <- pilotg(pp.5b.df)
  dens.pilot <- density.ppp(pp.5b, dimyx=c(128,128), varcov=G)
  m <- sum(dens.pilot/128^2)
  bmise.5b <- function(vh) {
    mat <- invvech(vh)
    ans <- 0
    for(i in 1:20) {
      N.star <- rpois(1, m)
      boot.index <- sample(1:nrow(pp.5b.df), N.star, replace = T)
      pp.5b.SB.df <- pp.5b.df[boot.index,] + as.data.frame(rmvnorm.mixt(n = N.star, Sigmas = G))
      pp.5b.SB <- as.ppp(pp.5b.SB.df, owin(c(0,1),c(0,1)))
      dens <- density.ppp(pp.5b.SB, dimyx=c(128,128), varcov=mat)
      ans <- ans + dens.dim * sum((dens.pilot$v - dens$v)^2) / sum(dens.pilot$v)^2
    }
    ans <- ans/N
    return(ans)
  }
  res.5b <- optim(vech(Hpi(pp.5b.df)), bmise.5b, method = "Nelder-Mead", 
                  control = list(trace = as.numeric(verbose), REPORT = 1))
  H.empboot.5b <- invvech(res.5b$par)
  dens.empboot.5b <- density.ppp(pp.5b, dimyx=c(128,128), varcov=H.empboot.5b)
  
  ISE.5b.empboot[j] <- F.computeISE(im.5b.128, dens.empboot.5b)
}



# 6b, empirical bootstrap

N <- 100
set.seed(887330)
verbose = FALSE

ISE.6b.empboot <- rep(NA, N)

im.6b.128 <- as.im(lambda6b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
dens.dim <- nrow(im.6b.128$v)*ncol(im.6b.128$v)

for(j in 1:N) {
  cat(j, '\t')
  pp.6b <- rpoispp(im.6b.128)
  pp.6b.df <- data.frame(x=pp.6b$x, y=pp.6b$y)
  G <- pilotg(pp.6b.df)
  dens.pilot <- density.ppp(pp.6b, dimyx=c(128,128), varcov=G)
  m <- sum(dens.pilot/128^2)
  bmise.6b <- function(vh) {
    mat <- invvech(vh)
    ans <- 0
    for(i in 1:20) {
      N.star <- rpois(1, m)
      boot.index <- sample(1:nrow(pp.6b.df), N.star, replace = T)
      pp.6b.SB.df <- pp.6b.df[boot.index,] + as.data.frame(rmvnorm.mixt(n = N.star, Sigmas = G))
      pp.6b.SB <- as.ppp(pp.6b.SB.df, owin(c(0,1),c(0,1)))
      dens <- density.ppp(pp.6b.SB, dimyx=c(128,128), varcov=mat)
      ans <- ans + dens.dim * sum((dens.pilot$v - dens$v)^2) / sum(dens.pilot$v)^2
    }
    ans <- ans/N
    return(ans)
  }
  res.6b <- optim(vech(Hpi(pp.6b.df)), bmise.6b, method = "Nelder-Mead", 
                  control = list(trace = as.numeric(verbose), REPORT = 1))
  H.empboot.6b <- invvech(res.6b$par)
  dens.empboot.6b <- density.ppp(pp.6b, dimyx=c(128,128), varcov=H.empboot.6b)
  
  ISE.6b.empboot[j] <- F.computeISE(im.6b.128, dens.empboot.6b)
}


save(ISE.1b.empboot,
     ISE.1b.thomas.empboot, 
     ISE.2b.empboot,
     ISE.2b.thomas.empboot,
     ISE.3b.empboot, 
     ISE.4b.empboot, 
     ISE.5b.empboot, 
     ISE.6b.empboot, 
     file = "ISEempboot.RData")
