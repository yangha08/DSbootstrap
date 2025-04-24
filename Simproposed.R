require(ks)
require(mclust)
require(spatstat)
require(RandomFields) # 4b

### Functions
### density.ppp : kernel intensity estimation
### F.computeISE : ISE between estimated density and true density


# 1b, scv

N <- 100

set.seed(887330)

MISE.mat.1b <- 0
MISE.diag.1b <- 0
MISE.mat.1b.1 <- 0
MISE.diag.1b.1 <- 0
MISE.mat.1b.1 <- 0
MISE.diag.1b.1 <- 0
MISE.mat.1b.2 <- 0
MISE.diag.1b.2 <- 0
MISE.mat.1b.3 <- 0
MISE.diag.1b.3 <- 0
MISE.mat.1b.4 <- 0
MISE.diag.1b.4 <- 0
MISE.mat.1b.5 <- 0
MISE.diag.1b.5 <- 0
MISE.mat.1b.10 <- 0
MISE.diag.1b.10 <- 0
mean.Bw.mat.1b <- diag(0,2,2)
mean.Bw.diag.1b <- diag(0,2,2)
mean.Bw.mat.1b.1 <- diag(0,2,2)
mean.Bw.diag.1b.1 <- diag(0,2,2)
mean.Bw.mat.1b.2 <- diag(0,2,2)
mean.Bw.diag.1b.2 <- diag(0,2,2)
mean.Bw.mat.1b.3 <- diag(0,2,2)
mean.Bw.diag.1b.3 <- diag(0,2,2)
mean.Bw.mat.1b.4 <- diag(0,2,2)
mean.Bw.diag.1b.4 <- diag(0,2,2)
mean.Bw.mat.1b.5 <- diag(0,2,2)
mean.Bw.diag.1b.5 <- diag(0,2,2)
mean.Bw.mat.1b.10 <- diag(0,2,2)
mean.Bw.diag.1b.10 <- diag(0,2,2)
time.1b <- 0

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.1b.128 <- as.im(lambda1b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.1b <- rpoispp(im.1b.128)
  pp.1b.df <- data.frame(x=pp.1b$x, y=pp.1b$y)
  Bw.mat <- Hscv(pp.1b.df)
  Bw.diag <- Hscv.diag(pp.1b.df)
  Bw.mat.1 <- Hscv.dir(pp.1b.df, alpha=1)
  Bw.diag.1 <- Hscv.diag.dir(pp.1b.df, alpha=1)
  Bw.mat.2 <- Hscv.dir(pp.1b.df, alpha=2)
  Bw.diag.2 <- Hscv.diag.dir(pp.1b.df, alpha=2)
  Bw.mat.3 <- Hscv.dir(pp.1b.df, alpha=3)
  Bw.diag.3 <- Hscv.diag.dir(pp.1b.df, alpha=3)
  Bw.mat.4 <- Hscv.dir(pp.1b.df, alpha=4)
  Bw.diag.4 <- Hscv.diag.dir(pp.1b.df, alpha=4)
  Bw.mat.5 <- Hscv.dir(pp.1b.df, alpha=5)
  Bw.diag.5 <- Hscv.diag.dir(pp.1b.df, alpha=5)
  Bw.mat.10 <- Hscv.dir(pp.1b.df, alpha=10)
  Bw.diag.10 <- Hscv.diag.dir(pp.1b.df, alpha=10)
  
  dens.128 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat)
  dens.diag.128 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag)
  dens.128.1 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat.1)
  dens.diag.128.1 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag.1)
  dens.128.2 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat.2)
  dens.diag.128.2 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag.2)
  dens.128.3 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat.3)
  dens.diag.128.3 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag.3)
  dens.128.4 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat.4)
  dens.diag.128.4 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag.4)
  dens.128.5 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat.5)
  dens.diag.128.5 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag.5)
  dens.128.10 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat.10)
  dens.diag.128.10 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag.10)
  
  mean.Bw.mat.1b <- mean.Bw.mat.1b + Bw.mat
  mean.Bw.diag.1b <- mean.Bw.diag.1b + Bw.diag
  MISE.mat.1b <- MISE.mat.1b + F.computeISE(im.1b.128, dens.128)
  MISE.diag.1b <- MISE.diag.1b + F.computeISE(im.1b.128, dens.diag.128)
  mean.Bw.mat.1b.1 <- mean.Bw.mat.1b.1 + Bw.mat.1
  mean.Bw.diag.1b.1 <- mean.Bw.diag.1b.1 + Bw.diag.1
  MISE.mat.1b.1 <- MISE.mat.1b.1 + F.computeISE(im.1b.128, dens.128.1)
  MISE.diag.1b.1 <- MISE.diag.1b.1 + F.computeISE(im.1b.128, dens.diag.128.1)
  mean.Bw.mat.1b.2 <- mean.Bw.mat.1b.2 + Bw.mat.2
  mean.Bw.diag.1b.2 <- mean.Bw.diag.1b.2 + Bw.diag.2
  MISE.mat.1b.2 <- MISE.mat.1b.2 + F.computeISE(im.1b.128, dens.128.2)
  MISE.diag.1b.2 <- MISE.diag.1b.2 + F.computeISE(im.1b.128, dens.diag.128.2)
  mean.Bw.mat.1b.3 <- mean.Bw.mat.1b.3 + Bw.mat.3
  mean.Bw.diag.1b.3 <- mean.Bw.diag.1b.3 + Bw.diag.3
  MISE.mat.1b.3 <- MISE.mat.1b.3 + F.computeISE(im.1b.128, dens.128.3)
  MISE.diag.1b.3 <- MISE.diag.1b.3 + F.computeISE(im.1b.128, dens.diag.128.3)
  mean.Bw.mat.1b.4 <- mean.Bw.mat.1b.4 + Bw.mat.4
  mean.Bw.diag.1b.4 <- mean.Bw.diag.1b.4 + Bw.diag.4
  MISE.mat.1b.4 <- MISE.mat.1b.4 + F.computeISE(im.1b.128, dens.128.4)
  MISE.diag.1b.4 <- MISE.diag.1b.4 + F.computeISE(im.1b.128, dens.diag.128.4)
  mean.Bw.mat.1b.5 <- mean.Bw.mat.1b.5 + Bw.mat.5
  mean.Bw.diag.1b.5 <- mean.Bw.diag.1b.5 + Bw.diag.5
  MISE.mat.1b.5 <- MISE.mat.1b.5 + F.computeISE(im.1b.128, dens.128.5)
  MISE.diag.1b.5 <- MISE.diag.1b.5 + F.computeISE(im.1b.128, dens.diag.128.5)
  mean.Bw.mat.1b.10 <- mean.Bw.mat.1b.10 + Bw.mat.10
  mean.Bw.diag.1b.10 <- mean.Bw.diag.1b.10 + Bw.diag.10
  MISE.mat.1b.10 <- MISE.mat.1b.10 + F.computeISE(im.1b.128, dens.128.10)
  MISE.diag.1b.10 <- MISE.diag.1b.10 + F.computeISE(im.1b.128, dens.diag.128.10)
  time.1b <- time.1b + Sys.time() - st

  # ise = F.computeISE(im.1b.128, dens.128)
  # MISE.mat.1b <- MISE.mat.1b + ise
  # ise.diag = F.computeISE(im.1b.128, dens.diag.128)
  # MISE.diag.1b <- MISE.diag.1b + ise.diag
  # cat(ise, ise.diag, "\n")
}

MISE.mat.1b <- MISE.mat.1b/N
MISE.diag.1b <- MISE.diag.1b/N
mean.Bw.mat.1b <- mean.Bw.mat.1b/N
mean.Bw.diag.1b <- mean.Bw.diag.1b/N
MISE.mat.1b.1 <- MISE.mat.1b.1/N
MISE.diag.1b.1 <- MISE.diag.1b.1/N
mean.Bw.mat.1b.1 <- mean.Bw.mat.1b.1/N
mean.Bw.diag.1b.1 <- mean.Bw.diag.1b.1/N
MISE.mat.1b.2 <- MISE.mat.1b.2/N
MISE.diag.1b.2 <- MISE.diag.1b.2/N
mean.Bw.mat.1b.2 <- mean.Bw.mat.1b.2/N
mean.Bw.diag.1b.2 <- mean.Bw.diag.1b.2/N
MISE.mat.1b.3 <- MISE.mat.1b.3/N
MISE.diag.1b.3 <- MISE.diag.1b.3/N
mean.Bw.mat.1b.3 <- mean.Bw.mat.1b.3/N
mean.Bw.diag.1b.3 <- mean.Bw.diag.1b.3/N
MISE.mat.1b.4 <- MISE.mat.1b.4/N
MISE.diag.1b.4 <- MISE.diag.1b.4/N
mean.Bw.mat.1b.4 <- mean.Bw.mat.1b.4/N
mean.Bw.diag.1b.4 <- mean.Bw.diag.1b.4/N
MISE.mat.1b.5 <- MISE.mat.1b.5/N
MISE.diag.1b.5 <- MISE.diag.1b.5/N
mean.Bw.mat.1b.5 <- mean.Bw.mat.1b.5/N
mean.Bw.diag.1b.5 <- mean.Bw.diag.1b.5/N
MISE.mat.1b.10 <- MISE.mat.1b.10/N
MISE.diag.1b.10 <- MISE.diag.1b.10/N
mean.Bw.mat.1b.10 <- mean.Bw.mat.1b.10/N
mean.Bw.diag.1b.10 <- mean.Bw.diag.1b.10/N



# 1b, pi

N <- 100

set.seed(887330)

MISE.mat.1b.pi <- 0
MISE.diag.1b.pi <- 0
MISE.mat.1b.pi.1 <- 0
MISE.diag.1b.pi.1 <- 0
MISE.mat.1b.pi.2 <- 0
MISE.diag.1b.pi.2 <- 0
MISE.mat.1b.pi.3 <- 0
MISE.diag.1b.pi.3 <- 0
MISE.mat.1b.pi.4 <- 0
MISE.diag.1b.pi.4 <- 0
MISE.mat.1b.pi.5 <- 0
MISE.diag.1b.pi.5 <- 0
MISE.mat.1b.pi.10 <- 0
MISE.diag.1b.pi.10 <- 0
mean.Bw.mat.1b.pi <- diag(0,2,2)
mean.Bw.diag.1b.pi <- diag(0,2,2)
mean.Bw.mat.1b.pi.1 <- diag(0,2,2)
mean.Bw.diag.1b.pi.1 <- diag(0,2,2)
mean.Bw.mat.1b.pi.2 <- diag(0,2,2)
mean.Bw.diag.1b.pi.2 <- diag(0,2,2)
mean.Bw.mat.1b.pi.3 <- diag(0,2,2)
mean.Bw.diag.1b.pi.3 <- diag(0,2,2)
mean.Bw.mat.1b.pi.4 <- diag(0,2,2)
mean.Bw.diag.1b.pi.4 <- diag(0,2,2)
mean.Bw.mat.1b.pi.5 <- diag(0,2,2)
mean.Bw.diag.1b.pi.5 <- diag(0,2,2)
mean.Bw.mat.1b.pi.10 <- diag(0,2,2)
mean.Bw.diag.1b.pi.10 <- diag(0,2,2)
time.1b.pi <- 0

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.1b.128 <- as.im(lambda1b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.1b <- rpoispp(im.1b.128)
  pp.1b.df <- data.frame(x=pp.1b$x, y=pp.1b$y)
  dens.pilot <- density.ppp(pp.1b, dimyx=c(128,128), varcov=pilot.G(pp.1b.df))
  Bw.pi.mat <- Hpi(pp.1b.df)
  Bw.pi.diag <- Hpi.diag(pp.1b.df)
  Bw.pi.mat.1 <- Hpi.dir(pp.1b.df, alpha=1)
  Bw.pi.diag.1 <- Hpi.diag.dir(pp.1b.df, alpha=1)
  Bw.pi.mat.2 <- Hpi.dir(pp.1b.df, alpha=2)
  Bw.pi.diag.2 <- Hpi.diag.dir(pp.1b.df, alpha=2)
  Bw.pi.mat.3 <- Hpi.dir(pp.1b.df, alpha=3)
  Bw.pi.diag.3 <- Hpi.diag.dir(pp.1b.df, alpha=3)
  Bw.pi.mat.4 <- Hpi.dir(pp.1b.df, alpha=4)
  Bw.pi.diag.4 <- Hpi.diag.dir(pp.1b.df, alpha=4)
  Bw.pi.mat.5 <- Hpi.dir(pp.1b.df, alpha=5)
  Bw.pi.diag.5 <- Hpi.diag.dir(pp.1b.df, alpha=5)
  Bw.pi.mat.10 <- Hpi.dir(pp.1b.df, alpha=10)
  Bw.pi.diag.10 <- Hpi.diag.dir(pp.1b.df, alpha=10)
  
  dens.pi.128 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.mat)
  dens.pi.diag.128 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.diag)
  dens.pi.128.1 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.mat.1)
  dens.pi.diag.128.1 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.diag.1)
  dens.pi.128.2 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.mat.2)
  dens.pi.diag.128.2 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.diag.2)
  dens.pi.128.3 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.mat.3)
  dens.pi.diag.128.3 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.diag.3)
  dens.pi.128.4 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.mat.4)
  dens.pi.diag.128.4 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.diag.4)
  dens.pi.128.5 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.mat.5)
  dens.pi.diag.128.5 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.diag.5)
  dens.pi.128.10 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.mat.10)
  dens.pi.diag.128.10 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.diag.10)
  
  mean.Bw.mat.1b.pi <- mean.Bw.mat.1b.pi + Bw.pi.mat
  mean.Bw.diag.1b.pi <- mean.Bw.diag.1b.pi + Bw.pi.diag
  MISE.mat.1b.pi <- MISE.mat.1b.pi + F.computeISE(im.1b.128, dens.pi.128)
  MISE.diag.1b.pi <- MISE.diag.1b.pi + F.computeISE(im.1b.128, dens.pi.diag.128)
  mean.Bw.mat.1b.pi.1 <- mean.Bw.mat.1b.pi.1 + Bw.pi.mat.1
  mean.Bw.diag.1b.pi.1 <- mean.Bw.diag.1b.pi.1 + Bw.pi.diag.1
  MISE.mat.1b.pi.1 <- MISE.mat.1b.pi.1 + F.computeISE(im.1b.128, dens.pi.128.1)
  MISE.diag.1b.pi.1 <- MISE.diag.1b.pi.1 + F.computeISE(im.1b.128, dens.pi.diag.128.1)
  mean.Bw.mat.1b.pi.2 <- mean.Bw.mat.1b.pi.2 + Bw.pi.mat.2
  mean.Bw.diag.1b.pi.2 <- mean.Bw.diag.1b.pi.2 + Bw.pi.diag.2
  MISE.mat.1b.pi.2 <- MISE.mat.1b.pi.2 + F.computeISE(im.1b.128, dens.pi.128.2)
  MISE.diag.1b.pi.2 <- MISE.diag.1b.pi.2 + F.computeISE(im.1b.128, dens.pi.diag.128.2)
  mean.Bw.mat.1b.pi.3 <- mean.Bw.mat.1b.pi.3 + Bw.pi.mat.3
  mean.Bw.diag.1b.pi.3 <- mean.Bw.diag.1b.pi.3 + Bw.pi.diag.3
  MISE.mat.1b.pi.3 <- MISE.mat.1b.pi.3 + F.computeISE(im.1b.128, dens.pi.128.3)
  MISE.diag.1b.pi.3 <- MISE.diag.1b.pi.3 + F.computeISE(im.1b.128, dens.pi.diag.128.3)
  mean.Bw.mat.1b.pi.4 <- mean.Bw.mat.1b.pi.4 + Bw.pi.mat.4
  mean.Bw.diag.1b.pi.4 <- mean.Bw.diag.1b.pi.4 + Bw.pi.diag.4
  MISE.mat.1b.pi.4 <- MISE.mat.1b.pi.4 + F.computeISE(im.1b.128, dens.pi.128.4)
  MISE.diag.1b.pi.4 <- MISE.diag.1b.pi.4 + F.computeISE(im.1b.128, dens.pi.diag.128.4)
  mean.Bw.mat.1b.pi.5 <- mean.Bw.mat.1b.pi.5 + Bw.pi.mat.5
  mean.Bw.diag.1b.pi.5 <- mean.Bw.diag.1b.pi.5 + Bw.pi.diag.5
  MISE.mat.1b.pi.5 <- MISE.mat.1b.pi.5 + F.computeISE(im.1b.128, dens.pi.128.5)
  MISE.diag.1b.pi.5 <- MISE.diag.1b.pi.5 + F.computeISE(im.1b.128, dens.pi.diag.128.5)
  mean.Bw.mat.1b.pi.10 <- mean.Bw.mat.1b.pi.10 + Bw.pi.mat.10
  mean.Bw.diag.1b.pi.10 <- mean.Bw.diag.1b.pi.10 + Bw.pi.diag.10
  MISE.mat.1b.pi.10 <- MISE.mat.1b.pi.10 + F.computeISE(im.1b.128, dens.pi.128.10)
  MISE.diag.1b.pi.10 <- MISE.diag.1b.pi.10 + F.computeISE(im.1b.128, dens.pi.diag.128.10)
  time.1b.pi.FS <- time.1b.pi.FS + Sys.time() - st
  
  # ise = F.computeISE(im.1b.128, dens.128)
  # MISE.mat.1b <- MISE.mat.1b + ise
  # ise.diag = F.computeISE(im.1b.128, dens.diag.128)
  # MISE.diag.1b <- MISE.diag.1b + ise.diag
  # cat(ise, ise.diag, "\n")
}

time.1b.pi <- time.1b.pi/N
MISE.mat.1b.pi <- MISE.mat.1b.pi/N
MISE.diag.1b.pi <- MISE.diag.1b.pi/N
mean.Bw.mat.1b.pi <- mean.Bw.mat.1b.pi/N
mean.Bw.diag.1b.pi <- mean.Bw.diag.1b.pi/N
MISE.mat.1b.pi.1 <- MISE.mat.1b.pi.1/N
MISE.diag.1b.pi.1 <- MISE.diag.1b.pi.1/N
mean.Bw.mat.1b.pi.1 <- mean.Bw.mat.1b.pi.1/N
mean.Bw.diag.1b.pi.1 <- mean.Bw.diag.1b.pi.1/N
MISE.mat.1b.pi.2 <- MISE.mat.1b.pi.2/N
MISE.diag.1b.pi.2 <- MISE.diag.1b.pi.2/N
mean.Bw.mat.1b.pi.2 <- mean.Bw.mat.1b.pi.2/N
mean.Bw.diag.1b.pi.2 <- mean.Bw.diag.1b.pi.2/N
MISE.mat.1b.pi.3 <- MISE.mat.1b.pi.3/N
MISE.diag.1b.pi.3 <- MISE.diag.1b.pi.3/N
mean.Bw.mat.1b.pi.3 <- mean.Bw.mat.1b.pi.3/N
mean.Bw.diag.1b.pi.3 <- mean.Bw.diag.1b.pi.3/N
MISE.mat.1b.pi.4 <- MISE.mat.1b.pi.4/N
MISE.diag.1b.pi.4 <- MISE.diag.1b.pi.4/N
mean.Bw.mat.1b.pi.4 <- mean.Bw.mat.1b.pi.4/N
mean.Bw.diag.1b.pi.4 <- mean.Bw.diag.1b.pi.4/N
MISE.mat.1b.pi.5 <- MISE.mat.1b.pi.5/N
MISE.diag.1b.pi.5 <- MISE.diag.1b.pi.5/N
mean.Bw.mat.1b.pi.5 <- mean.Bw.mat.1b.pi.5/N
mean.Bw.diag.1b.pi.5 <- mean.Bw.diag.1b.pi.5/N
MISE.mat.1b.pi.10 <- MISE.mat.1b.pi.10/N
MISE.diag.1b.pi.10 <- MISE.diag.1b.pi.10/N
mean.Bw.mat.1b.pi.10 <- mean.Bw.mat.1b.pi.10/N
mean.Bw.diag.1b.pi.10 <- mean.Bw.diag.1b.pi.10/N


# 1b Thomas, scv

N <- 100

set.seed(887330)

MISE.mat.1b.thomas <- 0
MISE.diag.1b.thomas <- 0
MISE.mat.1b.thomas.1 <- 0
MISE.diag.1b.thomas.1 <- 0
MISE.mat.1b.thomas.2 <- 0
MISE.diag.1b.thomas.2 <- 0
MISE.mat.1b.thomas.3 <- 0
MISE.diag.1b.thomas.3 <- 0
MISE.mat.1b.thomas.4 <- 0
MISE.diag.1b.thomas.4 <- 0
MISE.mat.1b.thomas.5 <- 0
MISE.diag.1b.thomas.5 <- 0
MISE.mat.1b.thomas.10 <- 0
MISE.diag.1b.thomas.10 <- 0
mean.Bw.mat.1b.thomas <- diag(0,2,2)
mean.Bw.diag.1b.thomas <- diag(0,2,2)
mean.Bw.mat.1b.thomas.1 <- diag(0,2,2)
mean.Bw.diag.1b.thomas.1 <- diag(0,2,2)
mean.Bw.mat.1b.thomas.2 <- diag(0,2,2)
mean.Bw.diag.1b.thomas.2 <- diag(0,2,2)
mean.Bw.mat.1b.thomas.3 <- diag(0,2,2)
mean.Bw.diag.1b.thomas.3 <- diag(0,2,2)
mean.Bw.mat.1b.thomas.4 <- diag(0,2,2)
mean.Bw.diag.1b.thomas.4 <- diag(0,2,2)
mean.Bw.mat.1b.thomas.5 <- diag(0,2,2)
mean.Bw.diag.1b.thomas.5 <- diag(0,2,2)
mean.Bw.mat.1b.thomas.10 <- diag(0,2,2)
mean.Bw.diag.1b.thomas.10 <- diag(0,2,2)

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.1b.128 <- as.im(lambda1b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.1b <- rThomas(kappa=200, scale=0.1, mu=5*im.1b.128*128*128/sum(im.1b.128))
  pp.1b.df <- data.frame(x=pp.1b$x, y=pp.1b$y)
  Bw.mat <- Hscv(pp.1b.df)
  Bw.diag <- Hscv.diag(pp.1b.df)
  Bw.mat.1 <- Hscv.dir(pp.1b.df, alpha=1)
  Bw.diag.1 <- Hscv.diag.dir(pp.1b.df, alpha=1)
  Bw.mat.2 <- Hscv.dir(pp.1b.df, alpha=2)
  Bw.diag.2 <- Hscv.diag.dir(pp.1b.df, alpha=2)
  Bw.mat.3 <- Hscv.dir(pp.1b.df, alpha=3)
  Bw.diag.3 <- Hscv.diag.dir(pp.1b.df, alpha=3)
  Bw.mat.4 <- Hscv.dir(pp.1b.df, alpha=4)
  Bw.diag.4 <- Hscv.diag.dir(pp.1b.df, alpha=4)
  Bw.mat.5 <- Hscv.dir(pp.1b.df, alpha=5)
  Bw.diag.5 <- Hscv.diag.dir(pp.1b.df, alpha=5)
  Bw.mat.10 <- Hscv.dir(pp.1b.df, alpha=10)
  Bw.diag.10 <- Hscv.diag.dir(pp.1b.df, alpha=10)
  
  dens.128 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat)
  dens.diag.128 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag)
  dens.128.1 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat.1)
  dens.diag.128.1 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag.1)
  dens.128.2 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat.2)
  dens.diag.128.2 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag.2)
  dens.128.3 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat.3)
  dens.diag.128.3 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag.3)
  dens.128.4 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat.4)
  dens.diag.128.4 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag.4)
  dens.128.5 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat.5)
  dens.diag.128.5 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag.5)
  dens.128.10 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat.10)
  dens.diag.128.10 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag.10)
  
  mean.Bw.mat.1b.thomas <- mean.Bw.mat.1b.thomas + Bw.mat
  mean.Bw.diag.1b.thomas <- mean.Bw.diag.1b.thomas + Bw.diag
  MISE.mat.1b.thomas <- MISE.mat.1b.thomas + F.computeISE(im.1b.128, dens.128)
  MISE.diag.1b.thomas <- MISE.diag.1b.thomas + F.computeISE(im.1b.128, dens.diag.128)
  mean.Bw.mat.1b.thomas.1 <- mean.Bw.mat.1b.thomas.1 + Bw.mat.1
  mean.Bw.diag.1b.thomas.1 <- mean.Bw.diag.1b.thomas.1 + Bw.diag.1
  MISE.mat.1b.thomas.1 <- MISE.mat.1b.thomas.1 + F.computeISE(im.1b.128, dens.128.1)
  MISE.diag.1b.thomas.1 <- MISE.diag.1b.thomas.1 + F.computeISE(im.1b.128, dens.diag.128.1)
  mean.Bw.mat.1b.thomas.2 <- mean.Bw.mat.1b.thomas.2 + Bw.mat.2
  mean.Bw.diag.1b.thomas.2 <- mean.Bw.diag.1b.thomas.2 + Bw.diag.2
  MISE.mat.1b.thomas.2 <- MISE.mat.1b.thomas.2 + F.computeISE(im.1b.128, dens.128.2)
  MISE.diag.1b.thomas.2 <- MISE.diag.1b.thomas.2 + F.computeISE(im.1b.128, dens.diag.128.2)
  mean.Bw.mat.1b.thomas.3 <- mean.Bw.mat.1b.thomas.3 + Bw.mat.3
  mean.Bw.diag.1b.thomas.3 <- mean.Bw.diag.1b.thomas.3 + Bw.diag.3
  MISE.mat.1b.thomas.3 <- MISE.mat.1b.thomas.3 + F.computeISE(im.1b.128, dens.128.3)
  MISE.diag.1b.thomas.3 <- MISE.diag.1b.thomas.3 + F.computeISE(im.1b.128, dens.diag.128.3)
  mean.Bw.mat.1b.thomas.4 <- mean.Bw.mat.1b.thomas.4 + Bw.mat.4
  mean.Bw.diag.1b.thomas.4 <- mean.Bw.diag.1b.thomas.4 + Bw.diag.4
  MISE.mat.1b.thomas.4 <- MISE.mat.1b.thomas.4 + F.computeISE(im.1b.128, dens.128.4)
  MISE.diag.1b.thomas.4 <- MISE.diag.1b.thomas.4 + F.computeISE(im.1b.128, dens.diag.128.4)
  mean.Bw.mat.1b.thomas.5 <- mean.Bw.mat.1b.thomas.5 + Bw.mat.5
  mean.Bw.diag.1b.thomas.5 <- mean.Bw.diag.1b.thomas.5 + Bw.diag.5
  MISE.mat.1b.thomas.5 <- MISE.mat.1b.thomas.5 + F.computeISE(im.1b.128, dens.128.5)
  MISE.diag.1b.thomas.5 <- MISE.diag.1b.thomas.5 + F.computeISE(im.1b.128, dens.diag.128.5)
  mean.Bw.mat.1b.thomas.10 <- mean.Bw.mat.1b.thomas.10 + Bw.mat.10
  mean.Bw.diag.1b.thomas.10 <- mean.Bw.diag.1b.thomas.10 + Bw.diag.10
  MISE.mat.1b.thomas.10 <- MISE.mat.1b.thomas.10 + F.computeISE(im.1b.128, dens.128.10)
  MISE.diag.1b.thomas.10 <- MISE.diag.1b.thomas.10 + F.computeISE(im.1b.128, dens.diag.128.10)
  time.1b.thomas <- time.1b.thomas + Sys.time() - st
  # ise = F.computeISE(im.1b.128, dens.128)
  # MISE.mat.1b <- MISE.mat.1b + ise
  # ise.diag = F.computeISE(im.1b.128, dens.diag.128)
  # MISE.diag.1b <- MISE.diag.1b + ise.diag
  # cat(ise, ise.diag, "\n")
}

time.1b.thomas <- time.1b.thomas/N
MISE.mat.1b.thomas <- MISE.mat.1b.thomas/N
MISE.diag.1b.thomas <- MISE.diag.1b.thomas/N
mean.Bw.mat.1b.thomas <- mean.Bw.mat.1b.thomas/N
mean.Bw.diag.1b.thomas <- mean.Bw.diag.1b.thomas/N
MISE.mat.1b.thomas.1 <- MISE.mat.1b.thomas.1/N
MISE.diag.1b.thomas.1 <- MISE.diag.1b.thomas.1/N
mean.Bw.mat.1b.thomas.1 <- mean.Bw.mat.1b.thomas.1/N
mean.Bw.diag.1b.thomas.1 <- mean.Bw.diag.1b.thomas.1/N
MISE.mat.1b.thomas.2 <- MISE.mat.1b.thomas.2/N
MISE.diag.1b.thomas.2 <- MISE.diag.1b.thomas.2/N
mean.Bw.mat.1b.thomas.2 <- mean.Bw.mat.1b.thomas.2/N
mean.Bw.diag.1b.thomas.2 <- mean.Bw.diag.1b.thomas.2/N
MISE.mat.1b.thomas.3 <- MISE.mat.1b.thomas.3/N
MISE.diag.1b.thomas.3 <- MISE.diag.1b.thomas.3/N
mean.Bw.mat.1b.thomas.3 <- mean.Bw.mat.1b.thomas.3/N
mean.Bw.diag.1b.thomas.3 <- mean.Bw.diag.1b.thomas.3/N
MISE.mat.1b.thomas.4 <- MISE.mat.1b.thomas.4/N
MISE.diag.1b.thomas.4 <- MISE.diag.1b.thomas.4/N
mean.Bw.mat.1b.thomas.4 <- mean.Bw.mat.1b.thomas.4/N
mean.Bw.diag.1b.thomas.4 <- mean.Bw.diag.1b.thomas.4/N
MISE.mat.1b.thomas.5 <- MISE.mat.1b.thomas.5/N
MISE.diag.1b.thomas.5 <- MISE.diag.1b.thomas.5/N
mean.Bw.mat.1b.thomas.5 <- mean.Bw.mat.1b.thomas.5/N
mean.Bw.diag.1b.thomas.5 <- mean.Bw.diag.1b.thomas.5/N
MISE.mat.1b.thomas.10 <- MISE.mat.1b.thomas.10/N
MISE.diag.1b.thomas.10 <- MISE.diag.1b.thomas.10/N
mean.Bw.mat.1b.thomas.10 <- mean.Bw.mat.1b.thomas.10/N
mean.Bw.diag.1b.thomas.10 <- mean.Bw.diag.1b.thomas.10/N




# 1b Thomas, pi

N <- 100

set.seed(887330)

MISE.mat.1b.thomas.pi <- 0
MISE.diag.1b.thomas.pi <- 0
MISE.mat.1b.thomas.pi.1 <- 0
MISE.diag.1b.thomas.pi.1 <- 0
MISE.mat.1b.thomas.pi.2 <- 0
MISE.diag.1b.thomas.pi.2 <- 0
MISE.mat.1b.thomas.pi.3 <- 0
MISE.diag.1b.thomas.pi.3 <- 0
MISE.mat.1b.thomas.pi.4 <- 0
MISE.diag.1b.thomas.pi.4 <- 0
MISE.mat.1b.thomas.pi.5 <- 0
MISE.diag.1b.thomas.pi.5 <- 0
MISE.mat.1b.thomas.pi.10 <- 0
MISE.diag.1b.thomas.pi.10 <- 0
mean.Bw.mat.1b.thomas.pi <- diag(0,2,2)
mean.Bw.diag.1b.thomas.pi <- diag(0,2,2)
mean.Bw.mat.1b.thomas.pi.1 <- diag(0,2,2)
mean.Bw.diag.1b.thomas.pi.1 <- diag(0,2,2)
mean.Bw.mat.1b.thomas.pi.2 <- diag(0,2,2)
mean.Bw.diag.1b.thomas.pi.2 <- diag(0,2,2)
mean.Bw.mat.1b.thomas.pi.3 <- diag(0,2,2)
mean.Bw.diag.1b.thomas.pi.3 <- diag(0,2,2)
mean.Bw.mat.1b.thomas.pi.4 <- diag(0,2,2)
mean.Bw.diag.1b.thomas.pi.4 <- diag(0,2,2)
mean.Bw.mat.1b.thomas.pi.5 <- diag(0,2,2)
mean.Bw.diag.1b.thomas.pi.5 <- diag(0,2,2)
mean.Bw.mat.1b.thomas.pi.10 <- diag(0,2,2)
mean.Bw.diag.1b.thomas.pi.10 <- diag(0,2,2)
time.1b.thomas.pi <- 0

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.1b.128 <- as.im(lambda1b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.1b <- rThomas(kappa=200, scale=0.1, mu=5*im.1b.128*128*128/sum(im.1b.128))
  pp.1b.df <- data.frame(x=pp.1b$x, y=pp.1b$y)
  Bw.pi.mat <- Hpi(pp.1b.df)
  Bw.pi.diag <- Hpi.diag(pp.1b.df)
  Bw.pi.mat.1 <- Hpi.dir(pp.1b.df, alpha=1)
  Bw.pi.diag.1 <- Hpi.diag.dir(pp.1b.df, alpha=1)
  Bw.pi.mat.2 <- Hpi.dir(pp.1b.df, alpha=2)
  Bw.pi.diag.2 <- Hpi.diag.dir(pp.1b.df, alpha=2)
  Bw.pi.mat.3 <- Hpi.dir(pp.1b.df, alpha=3)
  Bw.pi.diag.3 <- Hpi.diag.dir(pp.1b.df, alpha=3)
  Bw.pi.mat.4 <- Hpi.dir(pp.1b.df, alpha=4)
  Bw.pi.diag.4 <- Hpi.diag.dir(pp.1b.df, alpha=4)
  Bw.pi.mat.5 <- Hpi.dir(pp.1b.df, alpha=5)
  Bw.pi.diag.5 <- Hpi.diag.dir(pp.1b.df, alpha=5)
  Bw.pi.mat.10 <- Hpi.dir(pp.1b.df, alpha=10)
  Bw.pi.diag.10 <- Hpi.diag.dir(pp.1b.df, alpha=10)
  
  dens.pi.128 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.mat)
  dens.pi.diag.128 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.diag)
  dens.pi.128.1 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.mat.1)
  dens.pi.diag.128.1 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.diag.1)
  dens.pi.128.2 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.mat.2)
  dens.pi.diag.128.2 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.diag.2)
  dens.pi.128.3 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.mat.3)
  dens.pi.diag.128.3 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.diag.3)
  dens.pi.128.4 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.mat.4)
  dens.pi.diag.128.4 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.diag.4)
  dens.pi.128.5 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.mat.5)
  dens.pi.diag.128.5 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.diag.5)
  dens.pi.128.10 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.mat.10)
  dens.pi.diag.128.10 <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.pi.diag.10)
  
  mean.Bw.mat.1b.thomas.pi <- mean.Bw.mat.1b.thomas.pi + Bw.pi.mat
  mean.Bw.diag.1b.thomas.pi <- mean.Bw.diag.1b.thomas.pi + Bw.pi.diag
  MISE.mat.1b.thomas.pi <- MISE.mat.1b.thomas.pi + F.computeISE(im.1b.128, dens.pi.128)
  MISE.diag.1b.thomas.pi <- MISE.diag.1b.thomas.pi + F.computeISE(im.1b.128, dens.pi.diag.128)
  mean.Bw.mat.1b.thomas.pi.1 <- mean.Bw.mat.1b.thomas.pi.1 + Bw.pi.mat.1
  mean.Bw.diag.1b.thomas.pi.1 <- mean.Bw.diag.1b.thomas.pi.1 + Bw.pi.diag.1
  MISE.mat.1b.thomas.pi.1 <- MISE.mat.1b.thomas.pi.1 + F.computeISE(im.1b.128, dens.pi.128.1)
  MISE.diag.1b.thomas.pi.1 <- MISE.diag.1b.thomas.pi.1 + F.computeISE(im.1b.128, dens.pi.diag.128.1)
  mean.Bw.mat.1b.thomas.pi.2 <- mean.Bw.mat.1b.thomas.pi.2 + Bw.pi.mat.2
  mean.Bw.diag.1b.thomas.pi.2 <- mean.Bw.diag.1b.thomas.pi.2 + Bw.pi.diag.2
  MISE.mat.1b.thomas.pi.2 <- MISE.mat.1b.thomas.pi.2 + F.computeISE(im.1b.128, dens.pi.128.2)
  MISE.diag.1b.thomas.pi.2 <- MISE.diag.1b.thomas.pi.2 + F.computeISE(im.1b.128, dens.pi.diag.128.2)
  mean.Bw.mat.1b.thomas.pi.3 <- mean.Bw.mat.1b.thomas.pi.3 + Bw.pi.mat.3
  mean.Bw.diag.1b.thomas.pi.3 <- mean.Bw.diag.1b.thomas.pi.3 + Bw.pi.diag.3
  MISE.mat.1b.thomas.pi.3 <- MISE.mat.1b.thomas.pi.3 + F.computeISE(im.1b.128, dens.pi.128.3)
  MISE.diag.1b.thomas.pi.3 <- MISE.diag.1b.thomas.pi.3 + F.computeISE(im.1b.128, dens.pi.diag.128.3)
  mean.Bw.mat.1b.thomas.pi.4 <- mean.Bw.mat.1b.thomas.pi.4 + Bw.pi.mat.4
  mean.Bw.diag.1b.thomas.pi.4 <- mean.Bw.diag.1b.thomas.pi.4 + Bw.pi.diag.4
  MISE.mat.1b.thomas.pi.4 <- MISE.mat.1b.thomas.pi.4 + F.computeISE(im.1b.128, dens.pi.128.4)
  MISE.diag.1b.thomas.pi.4 <- MISE.diag.1b.thomas.pi.4 + F.computeISE(im.1b.128, dens.pi.diag.128.4)
  mean.Bw.mat.1b.thomas.pi.5 <- mean.Bw.mat.1b.thomas.pi.5 + Bw.pi.mat.5
  mean.Bw.diag.1b.thomas.pi.5 <- mean.Bw.diag.1b.thomas.pi.5 + Bw.pi.diag.5
  MISE.mat.1b.thomas.pi.5 <- MISE.mat.1b.thomas.pi.5 + F.computeISE(im.1b.128, dens.pi.128.5)
  MISE.diag.1b.thomas.pi.5 <- MISE.diag.1b.thomas.pi.5 + F.computeISE(im.1b.128, dens.pi.diag.128.5)
  mean.Bw.mat.1b.thomas.pi.10 <- mean.Bw.mat.1b.thomas.pi.10 + Bw.pi.mat.10
  mean.Bw.diag.1b.thomas.pi.10 <- mean.Bw.diag.1b.thomas.pi.10 + Bw.pi.diag.10
  MISE.mat.1b.thomas.pi.10 <- MISE.mat.1b.thomas.pi.10 + F.computeISE(im.1b.128, dens.pi.128.10)
  MISE.diag.1b.thomas.pi.10 <- MISE.diag.1b.thomas.pi.10 + F.computeISE(im.1b.128, dens.pi.diag.128.10)
  time.1b.thomas.pi <- time.1b.thomas.pi + Sys.time() - st
  
  # ise = F.computeISE(im.1b.128, dens.128)
  # MISE.mat.1b <- MISE.mat.1b + ise
  # ise.diag = F.computeISE(im.1b.128, dens.diag.128)
  # MISE.diag.1b <- MISE.diag.1b + ise.diag
  # cat(ise, ise.diag, "\n")
}

time.1b.thomas.pi <- time.1b.thomas.pi/N
MISE.mat.1b.thomas.pi <- MISE.mat.1b.thomas.pi/N
MISE.diag.1b.thomas.pi <- MISE.diag.1b.thomas.pi/N
mean.Bw.mat.1b.thomas.pi <- mean.Bw.mat.1b.thomas.pi/N
mean.Bw.diag.1b.thomas.pi <- mean.Bw.diag.1b.thomas.pi/N
MISE.mat.1b.thomas.pi.1 <- MISE.mat.1b.thomas.pi.1/N
MISE.diag.1b.thomas.pi.1 <- MISE.diag.1b.thomas.pi.1/N
mean.Bw.mat.1b.thomas.pi.1 <- mean.Bw.mat.1b.thomas.pi.1/N
mean.Bw.diag.1b.thomas.pi.1 <- mean.Bw.diag.1b.thomas.pi.1/N
MISE.mat.1b.thomas.pi.2 <- MISE.mat.1b.thomas.pi.2/N
MISE.diag.1b.thomas.pi.2 <- MISE.diag.1b.thomas.pi.2/N
mean.Bw.mat.1b.thomas.pi.2 <- mean.Bw.mat.1b.thomas.pi.2/N
mean.Bw.diag.1b.thomas.pi.2 <- mean.Bw.diag.1b.thomas.pi.2/N
MISE.mat.1b.thomas.pi.3 <- MISE.mat.1b.thomas.pi.3/N
MISE.diag.1b.thomas.pi.3 <- MISE.diag.1b.thomas.pi.3/N
mean.Bw.mat.1b.thomas.pi.3 <- mean.Bw.mat.1b.thomas.pi.3/N
mean.Bw.diag.1b.thomas.pi.3 <- mean.Bw.diag.1b.thomas.pi.3/N
MISE.mat.1b.thomas.pi.4 <- MISE.mat.1b.thomas.pi.4/N
MISE.diag.1b.thomas.pi.4 <- MISE.diag.1b.thomas.pi.4/N
mean.Bw.mat.1b.thomas.pi.4 <- mean.Bw.mat.1b.thomas.pi.4/N
mean.Bw.diag.1b.thomas.pi.4 <- mean.Bw.diag.1b.thomas.pi.4/N
MISE.mat.1b.thomas.pi.5 <- MISE.mat.1b.thomas.pi.5/N
MISE.diag.1b.thomas.pi.5 <- MISE.diag.1b.thomas.pi.5/N
mean.Bw.mat.1b.thomas.pi.5 <- mean.Bw.mat.1b.thomas.pi.5/N
mean.Bw.diag.1b.thomas.pi.5 <- mean.Bw.diag.1b.thomas.pi.5/N
MISE.mat.1b.thomas.pi.10 <- MISE.mat.1b.thomas.pi.10/N
MISE.diag.1b.thomas.pi.10 <- MISE.diag.1b.thomas.pi.10/N
mean.Bw.mat.1b.thomas.pi.10 <- mean.Bw.mat.1b.thomas.pi.10/N
mean.Bw.diag.1b.thomas.pi.10 <- mean.Bw.diag.1b.thomas.pi.10/N



# 2b, scv

N <- 100

set.seed(887330)

MISE.mat.2b <- 0
MISE.diag.2b <- 0
MISE.mat.2b.1 <- 0
MISE.diag.2b.1 <- 0
MISE.mat.2b.2 <- 0
MISE.diag.2b.2 <- 0
MISE.mat.2b.3 <- 0
MISE.diag.2b.3 <- 0
MISE.mat.2b.4 <- 0
MISE.diag.2b.4 <- 0
MISE.mat.2b.5 <- 0
MISE.diag.2b.5 <- 0
MISE.mat.2b.10 <- 0
MISE.diag.2b.10 <- 0
mean.Bw.mat.2b <- diag(0,2,2)
mean.Bw.diag.2b <- diag(0,2,2)
mean.Bw.mat.2b.1 <- diag(0,2,2)
mean.Bw.diag.2b.1 <- diag(0,2,2)
mean.Bw.mat.2b.2 <- diag(0,2,2)
mean.Bw.diag.2b.2 <- diag(0,2,2)
mean.Bw.mat.2b.3 <- diag(0,2,2)
mean.Bw.diag.2b.3 <- diag(0,2,2)
mean.Bw.mat.2b.4 <- diag(0,2,2)
mean.Bw.diag.2b.4 <- diag(0,2,2)
mean.Bw.mat.2b.5 <- diag(0,2,2)
mean.Bw.diag.2b.5 <- diag(0,2,2)
mean.Bw.mat.2b.10 <- diag(0,2,2)
mean.Bw.diag.2b.10 <- diag(0,2,2)

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.2b.128 <- as.im(lambda2b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.2b <- rpoispp(im.2b.128)
  pp.2b.df <- data.frame(x=pp.2b$x, y=pp.2b$y)
  Bw.mat <- Hscv(pp.2b.df)
  Bw.diag <- Hscv.diag(pp.2b.df)
  Bw.mat.1 <- Hscv.dir(pp.2b.df, alpha=1)
  Bw.diag.1 <- Hscv.diag.dir(pp.2b.df, alpha=1)
  Bw.mat.2 <- Hscv.dir(pp.2b.df, alpha=2)
  Bw.diag.2 <- Hscv.diag.dir(pp.2b.df, alpha=2)
  Bw.mat.3 <- Hscv.dir(pp.2b.df, alpha=3)
  Bw.diag.3 <- Hscv.diag.dir(pp.2b.df, alpha=3)
  Bw.mat.4 <- Hscv.dir(pp.2b.df, alpha=4)
  Bw.diag.4 <- Hscv.diag.dir(pp.2b.df, alpha=4)
  Bw.mat.5 <- Hscv.dir(pp.2b.df, alpha=5)
  Bw.diag.5 <- Hscv.diag.dir(pp.2b.df, alpha=5)
  Bw.mat.10 <- Hscv.dir(pp.2b.df, alpha=10)
  Bw.diag.10 <- Hscv.diag.dir(pp.2b.df, alpha=10)
  
  dens.128 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat)
  dens.diag.128 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag)
  dens.128.1 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat.1)
  dens.diag.128.1 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag.1)
  dens.128.2 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat.2)
  dens.diag.128.2 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag.2)
  dens.128.3 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat.3)
  dens.diag.128.3 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag.3)
  dens.128.4 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat.4)
  dens.diag.128.4 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag.4)
  dens.128.5 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat.5)
  dens.diag.128.5 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag.5)
  dens.128.10 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat.10)
  dens.diag.128.10 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag.10)
  
  mean.Bw.mat.2b <- mean.Bw.mat.2b + Bw.mat
  mean.Bw.diag.2b <- mean.Bw.diag.2b + Bw.diag
  MISE.mat.2b <- MISE.mat.2b + F.computeISE(im.2b.128, dens.128)
  MISE.diag.2b <- MISE.diag.2b + F.computeISE(im.2b.128, dens.diag.128)
  mean.Bw.mat.2b.1 <- mean.Bw.mat.2b.1 + Bw.mat.1
  mean.Bw.diag.2b.1 <- mean.Bw.diag.2b.1 + Bw.diag.1
  MISE.mat.2b.1 <- MISE.mat.2b.1 + F.computeISE(im.2b.128, dens.128.1)
  MISE.diag.2b.1 <- MISE.diag.2b.1 + F.computeISE(im.2b.128, dens.diag.128.1)
  mean.Bw.mat.2b.2 <- mean.Bw.mat.2b.2 + Bw.mat.2
  mean.Bw.diag.2b.2 <- mean.Bw.diag.2b.2 + Bw.diag.2
  MISE.mat.2b.2 <- MISE.mat.2b.2 + F.computeISE(im.2b.128, dens.128.2)
  MISE.diag.2b.2 <- MISE.diag.2b.2 + F.computeISE(im.2b.128, dens.diag.128.2)
  mean.Bw.mat.2b.3 <- mean.Bw.mat.2b.3 + Bw.mat.3
  mean.Bw.diag.2b.3 <- mean.Bw.diag.2b.3 + Bw.diag.3
  MISE.mat.2b.3 <- MISE.mat.2b.3 + F.computeISE(im.2b.128, dens.128.3)
  MISE.diag.2b.3 <- MISE.diag.2b.3 + F.computeISE(im.2b.128, dens.diag.128.3)
  mean.Bw.mat.2b.4 <- mean.Bw.mat.2b.4 + Bw.mat.4
  mean.Bw.diag.2b.4 <- mean.Bw.diag.2b.4 + Bw.diag.4
  MISE.mat.2b.4 <- MISE.mat.2b.4 + F.computeISE(im.2b.128, dens.128.4)
  MISE.diag.2b.4 <- MISE.diag.2b.4 + F.computeISE(im.2b.128, dens.diag.128.4)
  mean.Bw.mat.2b.5 <- mean.Bw.mat.2b.5 + Bw.mat.5
  mean.Bw.diag.2b.5 <- mean.Bw.diag.2b.5 + Bw.diag.5
  MISE.mat.2b.5 <- MISE.mat.2b.5 + F.computeISE(im.2b.128, dens.128.5)
  MISE.diag.2b.5 <- MISE.diag.2b.5 + F.computeISE(im.2b.128, dens.diag.128.5)
  mean.Bw.mat.2b.10 <- mean.Bw.mat.2b.10 + Bw.mat.10
  mean.Bw.diag.2b.10 <- mean.Bw.diag.2b.10 + Bw.diag.10
  MISE.mat.2b.10 <- MISE.mat.2b.10 + F.computeISE(im.2b.128, dens.128.10)
  MISE.diag.2b.10 <- MISE.diag.2b.10 + F.computeISE(im.2b.128, dens.diag.128.10)
  time.2b <- time.2b + Sys.time() - st
  
  # ise = F.computeISE(im.2b.128, dens.128)
  # MISE.mat.2b <- MISE.mat.2b + ise
  # ise.diag = F.computeISE(im.2b.128, dens.diag.128)
  # MISE.diag.2b <- MISE.diag.2b + ise.diag
  # cat(ise, ise.diag, "\n")
}

time.2b <- time.2b/N
MISE.mat.2b <- MISE.mat.2b/N
MISE.diag.2b <- MISE.diag.2b/N
mean.Bw.mat.2b <- mean.Bw.mat.2b/N
mean.Bw.diag.2b <- mean.Bw.diag.2b/N
MISE.mat.2b.1 <- MISE.mat.2b.1/N
MISE.diag.2b.1 <- MISE.diag.2b.1/N
mean.Bw.mat.2b.1 <- mean.Bw.mat.2b.1/N
mean.Bw.diag.2b.1 <- mean.Bw.diag.2b.1/N
MISE.mat.2b.2 <- MISE.mat.2b.2/N
MISE.diag.2b.2 <- MISE.diag.2b.2/N
mean.Bw.mat.2b.2 <- mean.Bw.mat.2b.2/N
mean.Bw.diag.2b.2 <- mean.Bw.diag.2b.2/N
MISE.mat.2b.3 <- MISE.mat.2b.3/N
MISE.diag.2b.3 <- MISE.diag.2b.3/N
mean.Bw.mat.2b.3 <- mean.Bw.mat.2b.3/N
mean.Bw.diag.2b.3 <- mean.Bw.diag.2b.3/N
MISE.mat.2b.4 <- MISE.mat.2b.4/N
MISE.diag.2b.4 <- MISE.diag.2b.4/N
mean.Bw.mat.2b.4 <- mean.Bw.mat.2b.4/N
mean.Bw.diag.2b.4 <- mean.Bw.diag.2b.4/N
MISE.mat.2b.5 <- MISE.mat.2b.5/N
MISE.diag.2b.5 <- MISE.diag.2b.5/N
mean.Bw.mat.2b.5 <- mean.Bw.mat.2b.5/N
mean.Bw.diag.2b.5 <- mean.Bw.diag.2b.5/N
MISE.mat.2b.10 <- MISE.mat.2b.10/N
MISE.diag.2b.10 <- MISE.diag.2b.10/N
mean.Bw.mat.2b.10 <- mean.Bw.mat.2b.10/N
mean.Bw.diag.2b.10 <- mean.Bw.diag.2b.10/N


# 2b, pi

N <- 100

set.seed(887330)

MISE.mat.2b.pi <- 0
MISE.diag.2b.pi <- 0
MISE.mat.2b.pi.1 <- 0
MISE.diag.2b.pi.1 <- 0
MISE.mat.2b.pi.2 <- 0
MISE.diag.2b.pi.2 <- 0
MISE.mat.2b.pi.3 <- 0
MISE.diag.2b.pi.3 <- 0
MISE.mat.2b.pi.4 <- 0
MISE.diag.2b.pi.4 <- 0
MISE.mat.2b.pi.5 <- 0
MISE.diag.2b.pi.5 <- 0
MISE.mat.2b.pi.10 <- 0
MISE.diag.2b.pi.10 <- 0
mean.Bw.mat.2b.pi <- diag(0,2,2)
mean.Bw.diag.2b.pi <- diag(0,2,2)
mean.Bw.mat.2b.pi.1 <- diag(0,2,2)
mean.Bw.diag.2b.pi.1 <- diag(0,2,2)
mean.Bw.mat.2b.pi.2 <- diag(0,2,2)
mean.Bw.diag.2b.pi.2 <- diag(0,2,2)
mean.Bw.mat.2b.pi.3 <- diag(0,2,2)
mean.Bw.diag.2b.pi.3 <- diag(0,2,2)
mean.Bw.mat.2b.pi.4 <- diag(0,2,2)
mean.Bw.diag.2b.pi.4 <- diag(0,2,2)
mean.Bw.mat.2b.pi.5 <- diag(0,2,2)
mean.Bw.diag.2b.pi.5 <- diag(0,2,2)
mean.Bw.mat.2b.pi.10 <- diag(0,2,2)
mean.Bw.diag.2b.pi.10 <- diag(0,2,2)
time.2b.pi <- 0

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.2b.128 <- as.im(lambda2b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.2b <- rpoispp(im.2b.128)
  pp.2b.df <- data.frame(x=pp.2b$x, y=pp.2b$y)
  Bw.pi.mat <- Hpi(pp.2b.df)
  Bw.pi.diag <- Hpi.diag(pp.2b.df)
  Bw.pi.mat.1 <- Hpi.dir(pp.2b.df, alpha=1)
  Bw.pi.diag.1 <- Hpi.diag.dir(pp.2b.df, alpha=1)
  Bw.pi.mat.2 <- Hpi.dir(pp.2b.df, alpha=2)
  Bw.pi.diag.2 <- Hpi.diag.dir(pp.2b.df, alpha=2)
  Bw.pi.mat.3 <- Hpi.dir(pp.2b.df, alpha=3)
  Bw.pi.diag.3 <- Hpi.diag.dir(pp.2b.df, alpha=3)
  Bw.pi.mat.4 <- Hpi.dir(pp.2b.df, alpha=4)
  Bw.pi.diag.4 <- Hpi.diag.dir(pp.2b.df, alpha=4)
  Bw.pi.mat.5 <- Hpi.dir(pp.2b.df, alpha=5)
  Bw.pi.diag.5 <- Hpi.diag.dir(pp.2b.df, alpha=5)
  Bw.pi.mat.10 <- Hpi.dir(pp.2b.df, alpha=10)
  Bw.pi.diag.10 <- Hpi.diag.dir(pp.2b.df, alpha=10)
  
  dens.pi.128 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.mat)
  dens.pi.diag.128 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.diag)
  dens.pi.128.1 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.mat.1)
  dens.pi.diag.128.1 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.diag.1)
  dens.pi.128.2 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.mat.2)
  dens.pi.diag.128.2 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.diag.2)
  dens.pi.128.3 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.mat.3)
  dens.pi.diag.128.3 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.diag.3)
  dens.pi.128.4 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.mat.4)
  dens.pi.diag.128.4 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.diag.4)
  dens.pi.128.5 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.mat.5)
  dens.pi.diag.128.5 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.diag.5)
  dens.pi.128.10 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.mat.10)
  dens.pi.diag.128.10 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.diag.10)
  
  mean.Bw.mat.2b.pi <- mean.Bw.mat.2b.pi + Bw.pi.mat
  mean.Bw.diag.2b.pi <- mean.Bw.diag.2b.pi + Bw.pi.diag
  MISE.mat.2b.pi <- MISE.mat.2b.pi + F.computeISE(im.2b.128, dens.pi.128)
  MISE.diag.2b.pi <- MISE.diag.2b.pi + F.computeISE(im.2b.128, dens.pi.diag.128)
  mean.Bw.mat.2b.pi.1 <- mean.Bw.mat.2b.pi.1 + Bw.pi.mat.1
  mean.Bw.diag.2b.pi.1 <- mean.Bw.diag.2b.pi.1 + Bw.pi.diag.1
  MISE.mat.2b.pi.1 <- MISE.mat.2b.pi.1 + F.computeISE(im.2b.128, dens.pi.128.1)
  MISE.diag.2b.pi.1 <- MISE.diag.2b.pi.1 + F.computeISE(im.2b.128, dens.pi.diag.128.1)
  mean.Bw.mat.2b.pi.2 <- mean.Bw.mat.2b.pi.2 + Bw.pi.mat.2
  mean.Bw.diag.2b.pi.2 <- mean.Bw.diag.2b.pi.2 + Bw.pi.diag.2
  MISE.mat.2b.pi.2 <- MISE.mat.2b.pi.2 + F.computeISE(im.2b.128, dens.pi.128.2)
  MISE.diag.2b.pi.2 <- MISE.diag.2b.pi.2 + F.computeISE(im.2b.128, dens.pi.diag.128.2)
  mean.Bw.mat.2b.pi.3 <- mean.Bw.mat.2b.pi.3 + Bw.pi.mat.3
  mean.Bw.diag.2b.pi.3 <- mean.Bw.diag.2b.pi.3 + Bw.pi.diag.3
  MISE.mat.2b.pi.3 <- MISE.mat.2b.pi.3 + F.computeISE(im.2b.128, dens.pi.128.3)
  MISE.diag.2b.pi.3 <- MISE.diag.2b.pi.3 + F.computeISE(im.2b.128, dens.pi.diag.128.3)
  mean.Bw.mat.2b.pi.4 <- mean.Bw.mat.2b.pi.4 + Bw.pi.mat.4
  mean.Bw.diag.2b.pi.4 <- mean.Bw.diag.2b.pi.4 + Bw.pi.diag.4
  MISE.mat.2b.pi.4 <- MISE.mat.2b.pi.4 + F.computeISE(im.2b.128, dens.pi.128.4)
  MISE.diag.2b.pi.4 <- MISE.diag.2b.pi.4 + F.computeISE(im.2b.128, dens.pi.diag.128.4)
  mean.Bw.mat.2b.pi.5 <- mean.Bw.mat.2b.pi.5 + Bw.pi.mat.5
  mean.Bw.diag.2b.pi.5 <- mean.Bw.diag.2b.pi.5 + Bw.pi.diag.5
  MISE.mat.2b.pi.5 <- MISE.mat.2b.pi.5 + F.computeISE(im.2b.128, dens.pi.128.5)
  MISE.diag.2b.pi.5 <- MISE.diag.2b.pi.5 + F.computeISE(im.2b.128, dens.pi.diag.128.5)
  mean.Bw.mat.2b.pi.10 <- mean.Bw.mat.2b.pi.10 + Bw.pi.mat.10
  mean.Bw.diag.2b.pi.10 <- mean.Bw.diag.2b.pi.10 + Bw.pi.diag.10
  MISE.mat.2b.pi.10 <- MISE.mat.2b.pi.10 + F.computeISE(im.2b.128, dens.pi.128.10)
  MISE.diag.2b.pi.10 <- MISE.diag.2b.pi.10 + F.computeISE(im.2b.128, dens.pi.diag.128.10)
  time.2b.pi <- time.2b.pi + Sys.time() - st
  
  # ise = F.computeISE(im.2b.128, dens.128)
  # MISE.mat.2b <- MISE.mat.2b + ise
  # ise.diag = F.computeISE(im.2b.128, dens.diag.128)
  # MISE.diag.2b <- MISE.diag.2b + ise.diag
  # cat(ise, ise.diag, "\n")
}

time.2b.pi <- time.2b.pi/N
MISE.mat.2b.pi <- MISE.mat.2b.pi/N
MISE.diag.2b.pi <- MISE.diag.2b.pi/N
mean.Bw.mat.2b.pi <- mean.Bw.mat.2b.pi/N
mean.Bw.diag.2b.pi <- mean.Bw.diag.2b.pi/N
MISE.mat.2b.pi.1 <- MISE.mat.2b.pi.1/N
MISE.diag.2b.pi.1 <- MISE.diag.2b.pi.1/N
mean.Bw.mat.2b.pi.1 <- mean.Bw.mat.2b.pi.1/N
mean.Bw.diag.2b.pi.1 <- mean.Bw.diag.2b.pi.1/N
MISE.mat.2b.pi.2 <- MISE.mat.2b.pi.2/N
MISE.diag.2b.pi.2 <- MISE.diag.2b.pi.2/N
mean.Bw.mat.2b.pi.2 <- mean.Bw.mat.2b.pi.2/N
mean.Bw.diag.2b.pi.2 <- mean.Bw.diag.2b.pi.2/N
MISE.mat.2b.pi.3 <- MISE.mat.2b.pi.3/N
MISE.diag.2b.pi.3 <- MISE.diag.2b.pi.3/N
mean.Bw.mat.2b.pi.3 <- mean.Bw.mat.2b.pi.3/N
mean.Bw.diag.2b.pi.3 <- mean.Bw.diag.2b.pi.3/N
MISE.mat.2b.pi.4 <- MISE.mat.2b.pi.4/N
MISE.diag.2b.pi.4 <- MISE.diag.2b.pi.4/N
mean.Bw.mat.2b.pi.4 <- mean.Bw.mat.2b.pi.4/N
mean.Bw.diag.2b.pi.4 <- mean.Bw.diag.2b.pi.4/N
MISE.mat.2b.pi.5 <- MISE.mat.2b.pi.5/N
MISE.diag.2b.pi.5 <- MISE.diag.2b.pi.5/N
mean.Bw.mat.2b.pi.5 <- mean.Bw.mat.2b.pi.5/N
mean.Bw.diag.2b.pi.5 <- mean.Bw.diag.2b.pi.5/N
MISE.mat.2b.pi.10 <- MISE.mat.2b.pi.10/N
MISE.diag.2b.pi.10 <- MISE.diag.2b.pi.10/N
mean.Bw.mat.2b.pi.10 <- mean.Bw.mat.2b.pi.10/N
mean.Bw.diag.2b.pi.10 <- mean.Bw.diag.2b.pi.10/N



# 2b Thomas, scv

N <- 100

set.seed(887330)

MISE.mat.2b.thomas <- 0
MISE.diag.2b.thomas <- 0
MISE.mat.2b.thomas.1 <- 0
MISE.diag.2b.thomas.1 <- 0
MISE.mat.2b.thomas.2 <- 0
MISE.diag.2b.thomas.2 <- 0
MISE.mat.2b.thomas.3 <- 0
MISE.diag.2b.thomas.3 <- 0
MISE.mat.2b.thomas.4 <- 0
MISE.diag.2b.thomas.4 <- 0
MISE.mat.2b.thomas.5 <- 0
MISE.diag.2b.thomas.5 <- 0
MISE.mat.2b.thomas.10 <- 0
MISE.diag.2b.thomas.10 <- 0
mean.Bw.mat.2b.thomas <- diag(0,2,2)
mean.Bw.diag.2b.thomas <- diag(0,2,2)
mean.Bw.mat.2b.thomas.1 <- diag(0,2,2)
mean.Bw.diag.2b.thomas.1 <- diag(0,2,2)
mean.Bw.mat.2b.thomas.2 <- diag(0,2,2)
mean.Bw.diag.2b.thomas.2 <- diag(0,2,2)
mean.Bw.mat.2b.thomas.3 <- diag(0,2,2)
mean.Bw.diag.2b.thomas.3 <- diag(0,2,2)
mean.Bw.mat.2b.thomas.4 <- diag(0,2,2)
mean.Bw.diag.2b.thomas.4 <- diag(0,2,2)
mean.Bw.mat.2b.thomas.5 <- diag(0,2,2)
mean.Bw.diag.2b.thomas.5 <- diag(0,2,2)
mean.Bw.mat.2b.thomas.10 <- diag(0,2,2)
mean.Bw.diag.2b.thomas.10 <- diag(0,2,2)

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.2b.128 <- as.im(lambda2b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.2b <- rThomas(kappa=200, scale=0.1, mu=5*im.2b.128*128*128/sum(im.2b.128))
  pp.2b.df <- data.frame(x=pp.2b$x, y=pp.2b$y)
  Bw.mat <- Hscv(pp.2b.df)
  Bw.diag <- Hscv.diag(pp.2b.df)
  Bw.mat.1 <- Hscv.dir(pp.2b.df, alpha=1)
  Bw.diag.1 <- Hscv.diag.dir(pp.2b.df, alpha=1)
  Bw.mat.2 <- Hscv.dir(pp.2b.df, alpha=2)
  Bw.diag.2 <- Hscv.diag.dir(pp.2b.df, alpha=2)
  Bw.mat.3 <- Hscv.dir(pp.2b.df, alpha=3)
  Bw.diag.3 <- Hscv.diag.dir(pp.2b.df, alpha=3)
  Bw.mat.4 <- Hscv.dir(pp.2b.df, alpha=4)
  Bw.diag.4 <- Hscv.diag.dir(pp.2b.df, alpha=4)
  Bw.mat.5 <- Hscv.dir(pp.2b.df, alpha=5)
  Bw.diag.5 <- Hscv.diag.dir(pp.2b.df, alpha=5)
  Bw.mat.10 <- Hscv.dir(pp.2b.df, alpha=10)
  Bw.diag.10 <- Hscv.diag.dir(pp.2b.df, alpha=10)
  
  dens.128 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat)
  dens.diag.128 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag)
  dens.128.1 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat.1)
  dens.diag.128.1 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag.1)
  dens.128.2 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat.2)
  dens.diag.128.2 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag.2)
  dens.128.3 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat.3)
  dens.diag.128.3 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag.3)
  dens.128.4 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat.4)
  dens.diag.128.4 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag.4)
  dens.128.5 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat.5)
  dens.diag.128.5 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag.5)
  dens.128.10 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat.10)
  dens.diag.128.10 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag.10)
  
  mean.Bw.mat.2b.thomas <- mean.Bw.mat.2b.thomas + Bw.mat
  mean.Bw.diag.2b.thomas <- mean.Bw.diag.2b.thomas + Bw.diag
  MISE.mat.2b.thomas <- MISE.mat.2b.thomas + F.computeISE(im.2b.128, dens.128)
  MISE.diag.2b.thomas <- MISE.diag.2b.thomas + F.computeISE(im.2b.128, dens.diag.128)
  mean.Bw.mat.2b.thomas.1 <- mean.Bw.mat.2b.thomas.1 + Bw.mat.1
  mean.Bw.diag.2b.thomas.1 <- mean.Bw.diag.2b.thomas.1 + Bw.diag.1
  MISE.mat.2b.thomas.1 <- MISE.mat.2b.thomas.1 + F.computeISE(im.2b.128, dens.128.1)
  MISE.diag.2b.thomas.1 <- MISE.diag.2b.thomas.1 + F.computeISE(im.2b.128, dens.diag.128.1)
  mean.Bw.mat.2b.thomas.2 <- mean.Bw.mat.2b.thomas.2 + Bw.mat.2
  mean.Bw.diag.2b.thomas.2 <- mean.Bw.diag.2b.thomas.2 + Bw.diag.2
  MISE.mat.2b.thomas.2 <- MISE.mat.2b.thomas.2 + F.computeISE(im.2b.128, dens.128.2)
  MISE.diag.2b.thomas.2 <- MISE.diag.2b.thomas.2 + F.computeISE(im.2b.128, dens.diag.128.2)
  mean.Bw.mat.2b.thomas.3 <- mean.Bw.mat.2b.thomas.3 + Bw.mat.3
  mean.Bw.diag.2b.thomas.3 <- mean.Bw.diag.2b.thomas.3 + Bw.diag.3
  MISE.mat.2b.thomas.3 <- MISE.mat.2b.thomas.3 + F.computeISE(im.2b.128, dens.128.3)
  MISE.diag.2b.thomas.3 <- MISE.diag.2b.thomas.3 + F.computeISE(im.2b.128, dens.diag.128.3)
  mean.Bw.mat.2b.thomas.4 <- mean.Bw.mat.2b.thomas.4 + Bw.mat.4
  mean.Bw.diag.2b.thomas.4 <- mean.Bw.diag.2b.thomas.4 + Bw.diag.4
  MISE.mat.2b.thomas.4 <- MISE.mat.2b.thomas.4 + F.computeISE(im.2b.128, dens.128.4)
  MISE.diag.2b.thomas.4 <- MISE.diag.2b.thomas.4 + F.computeISE(im.2b.128, dens.diag.128.4)
  mean.Bw.mat.2b.thomas.5 <- mean.Bw.mat.2b.thomas.5 + Bw.mat.5
  mean.Bw.diag.2b.thomas.5 <- mean.Bw.diag.2b.thomas.5 + Bw.diag.5
  MISE.mat.2b.thomas.5 <- MISE.mat.2b.thomas.5 + F.computeISE(im.2b.128, dens.128.5)
  MISE.diag.2b.thomas.5 <- MISE.diag.2b.thomas.5 + F.computeISE(im.2b.128, dens.diag.128.5)
  mean.Bw.mat.2b.thomas.10 <- mean.Bw.mat.2b.thomas.10 + Bw.mat.10
  mean.Bw.diag.2b.thomas.10 <- mean.Bw.diag.2b.thomas.10 + Bw.diag.10
  MISE.mat.2b.thomas.10 <- MISE.mat.2b.thomas.10 + F.computeISE(im.2b.128, dens.128.10)
  MISE.diag.2b.thomas.10 <- MISE.diag.2b.thomas.10 + F.computeISE(im.2b.128, dens.diag.128.10)
  time.2b.thomas <- time.2b.thomas + Sys.time() - st
  # ise = F.computeISE(im.2b.128, dens.128)
  # MISE.mat.2b <- MISE.mat.2b + ise
  # ise.diag = F.computeISE(im.2b.128, dens.diag.128)
  # MISE.diag.2b <- MISE.diag.2b + ise.diag
  # cat(ise, ise.diag, "\n")
}

time.2b.thomas <- time.2b.thomas/N
MISE.mat.2b.thomas <- MISE.mat.2b.thomas/N
MISE.diag.2b.thomas <- MISE.diag.2b.thomas/N
mean.Bw.mat.2b.thomas <- mean.Bw.mat.2b.thomas/N
mean.Bw.diag.2b.thomas <- mean.Bw.diag.2b.thomas/N
MISE.mat.2b.thomas.1 <- MISE.mat.2b.thomas.1/N
MISE.diag.2b.thomas.1 <- MISE.diag.2b.thomas.1/N
mean.Bw.mat.2b.thomas.1 <- mean.Bw.mat.2b.thomas.1/N
mean.Bw.diag.2b.thomas.1 <- mean.Bw.diag.2b.thomas.1/N
MISE.mat.2b.thomas.2 <- MISE.mat.2b.thomas.2/N
MISE.diag.2b.thomas.2 <- MISE.diag.2b.thomas.2/N
mean.Bw.mat.2b.thomas.2 <- mean.Bw.mat.2b.thomas.2/N
mean.Bw.diag.2b.thomas.2 <- mean.Bw.diag.2b.thomas.2/N
MISE.mat.2b.thomas.3 <- MISE.mat.2b.thomas.3/N
MISE.diag.2b.thomas.3 <- MISE.diag.2b.thomas.3/N
mean.Bw.mat.2b.thomas.3 <- mean.Bw.mat.2b.thomas.3/N
mean.Bw.diag.2b.thomas.3 <- mean.Bw.diag.2b.thomas.3/N
MISE.mat.2b.thomas.4 <- MISE.mat.2b.thomas.4/N
MISE.diag.2b.thomas.4 <- MISE.diag.2b.thomas.4/N
mean.Bw.mat.2b.thomas.4 <- mean.Bw.mat.2b.thomas.4/N
mean.Bw.diag.2b.thomas.4 <- mean.Bw.diag.2b.thomas.4/N
MISE.mat.2b.thomas.5 <- MISE.mat.2b.thomas.5/N
MISE.diag.2b.thomas.5 <- MISE.diag.2b.thomas.5/N
mean.Bw.mat.2b.thomas.5 <- mean.Bw.mat.2b.thomas.5/N
mean.Bw.diag.2b.thomas.5 <- mean.Bw.diag.2b.thomas.5/N
MISE.mat.2b.thomas.10 <- MISE.mat.2b.thomas.10/N
MISE.diag.2b.thomas.10 <- MISE.diag.2b.thomas.10/N
mean.Bw.mat.2b.thomas.10 <- mean.Bw.mat.2b.thomas.10/N
mean.Bw.diag.2b.thomas.10 <- mean.Bw.diag.2b.thomas.10/N




# 2b Thomas, pi

N <- 100

set.seed(887330)

MISE.mat.2b.thomas.pi <- 0
MISE.diag.2b.thomas.pi <- 0
MISE.mat.2b.thomas.pi.1 <- 0
MISE.diag.2b.thomas.pi.1 <- 0
MISE.mat.2b.thomas.pi.2 <- 0
MISE.diag.2b.thomas.pi.2 <- 0
MISE.mat.2b.thomas.pi.3 <- 0
MISE.diag.2b.thomas.pi.3 <- 0
MISE.mat.2b.thomas.pi.4 <- 0
MISE.diag.2b.thomas.pi.4 <- 0
MISE.mat.2b.thomas.pi.5 <- 0
MISE.diag.2b.thomas.pi.5 <- 0
MISE.mat.2b.thomas.pi.10 <- 0
MISE.diag.2b.thomas.pi.10 <- 0
mean.Bw.mat.2b.thomas.pi <- diag(0,2,2)
mean.Bw.diag.2b.thomas.pi <- diag(0,2,2)
mean.Bw.mat.2b.thomas.pi.1 <- diag(0,2,2)
mean.Bw.diag.2b.thomas.pi.1 <- diag(0,2,2)
mean.Bw.mat.2b.thomas.pi.2 <- diag(0,2,2)
mean.Bw.diag.2b.thomas.pi.2 <- diag(0,2,2)
mean.Bw.mat.2b.thomas.pi.3 <- diag(0,2,2)
mean.Bw.diag.2b.thomas.pi.3 <- diag(0,2,2)
mean.Bw.mat.2b.thomas.pi.4 <- diag(0,2,2)
mean.Bw.diag.2b.thomas.pi.4 <- diag(0,2,2)
mean.Bw.mat.2b.thomas.pi.5 <- diag(0,2,2)
mean.Bw.diag.2b.thomas.pi.5 <- diag(0,2,2)
mean.Bw.mat.2b.thomas.pi.10 <- diag(0,2,2)
mean.Bw.diag.2b.thomas.pi.10 <- diag(0,2,2)
time.2b.thomas.pi <- 0

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.2b.128 <- as.im(lambda2b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.2b <- rThomas(kappa=200, scale=0.1, mu=5*im.2b.128*128*128/sum(im.2b.128))
  pp.2b.df <- data.frame(x=pp.2b$x, y=pp.2b$y)
  Bw.pi.mat <- Hpi(pp.2b.df)
  Bw.pi.diag <- Hpi.diag(pp.2b.df)
  Bw.pi.mat.1 <- Hpi.dir(pp.2b.df, alpha=1)
  Bw.pi.diag.1 <- Hpi.diag.dir(pp.2b.df, alpha=1)
  Bw.pi.mat.2 <- Hpi.dir(pp.2b.df, alpha=2)
  Bw.pi.diag.2 <- Hpi.diag.dir(pp.2b.df, alpha=2)
  Bw.pi.mat.3 <- Hpi.dir(pp.2b.df, alpha=3)
  Bw.pi.diag.3 <- Hpi.diag.dir(pp.2b.df, alpha=3)
  Bw.pi.mat.4 <- Hpi.dir(pp.2b.df, alpha=4)
  Bw.pi.diag.4 <- Hpi.diag.dir(pp.2b.df, alpha=4)
  Bw.pi.mat.5 <- Hpi.dir(pp.2b.df, alpha=5)
  Bw.pi.diag.5 <- Hpi.diag.dir(pp.2b.df, alpha=5)
  Bw.pi.mat.10 <- Hpi.dir(pp.2b.df, alpha=10)
  Bw.pi.diag.10 <- Hpi.diag.dir(pp.2b.df, alpha=10)
  
  dens.pi.128 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.mat)
  dens.pi.diag.128 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.diag)
  dens.pi.128.1 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.mat.1)
  dens.pi.diag.128.1 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.diag.1)
  dens.pi.128.2 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.mat.2)
  dens.pi.diag.128.2 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.diag.2)
  dens.pi.128.3 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.mat.3)
  dens.pi.diag.128.3 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.diag.3)
  dens.pi.128.4 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.mat.4)
  dens.pi.diag.128.4 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.diag.4)
  dens.pi.128.5 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.mat.5)
  dens.pi.diag.128.5 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.diag.5)
  dens.pi.128.10 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.mat.10)
  dens.pi.diag.128.10 <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.pi.diag.10)
  
  mean.Bw.mat.2b.thomas.pi <- mean.Bw.mat.2b.thomas.pi + Bw.pi.mat
  mean.Bw.diag.2b.thomas.pi <- mean.Bw.diag.2b.thomas.pi + Bw.pi.diag
  MISE.mat.2b.thomas.pi <- MISE.mat.2b.thomas.pi + F.computeISE(im.2b.128, dens.pi.128)
  MISE.diag.2b.thomas.pi <- MISE.diag.2b.thomas.pi + F.computeISE(im.2b.128, dens.pi.diag.128)
  mean.Bw.mat.2b.thomas.pi.1 <- mean.Bw.mat.2b.thomas.pi.1 + Bw.pi.mat.1
  mean.Bw.diag.2b.thomas.pi.1 <- mean.Bw.diag.2b.thomas.pi.1 + Bw.pi.diag.1
  MISE.mat.2b.thomas.pi.1 <- MISE.mat.2b.thomas.pi.1 + F.computeISE(im.2b.128, dens.pi.128.1)
  MISE.diag.2b.thomas.pi.1 <- MISE.diag.2b.thomas.pi.1 + F.computeISE(im.2b.128, dens.pi.diag.128.1)
  mean.Bw.mat.2b.thomas.pi.2 <- mean.Bw.mat.2b.thomas.pi.2 + Bw.pi.mat.2
  mean.Bw.diag.2b.thomas.pi.2 <- mean.Bw.diag.2b.thomas.pi.2 + Bw.pi.diag.2
  MISE.mat.2b.thomas.pi.2 <- MISE.mat.2b.thomas.pi.2 + F.computeISE(im.2b.128, dens.pi.128.2)
  MISE.diag.2b.thomas.pi.2 <- MISE.diag.2b.thomas.pi.2 + F.computeISE(im.2b.128, dens.pi.diag.128.2)
  mean.Bw.mat.2b.thomas.pi.3 <- mean.Bw.mat.2b.thomas.pi.3 + Bw.pi.mat.3
  mean.Bw.diag.2b.thomas.pi.3 <- mean.Bw.diag.2b.thomas.pi.3 + Bw.pi.diag.3
  MISE.mat.2b.thomas.pi.3 <- MISE.mat.2b.thomas.pi.3 + F.computeISE(im.2b.128, dens.pi.128.3)
  MISE.diag.2b.thomas.pi.3 <- MISE.diag.2b.thomas.pi.3 + F.computeISE(im.2b.128, dens.pi.diag.128.3)
  mean.Bw.mat.2b.thomas.pi.4 <- mean.Bw.mat.2b.thomas.pi.4 + Bw.pi.mat.4
  mean.Bw.diag.2b.thomas.pi.4 <- mean.Bw.diag.2b.thomas.pi.4 + Bw.pi.diag.4
  MISE.mat.2b.thomas.pi.4 <- MISE.mat.2b.thomas.pi.4 + F.computeISE(im.2b.128, dens.pi.128.4)
  MISE.diag.2b.thomas.pi.4 <- MISE.diag.2b.thomas.pi.4 + F.computeISE(im.2b.128, dens.pi.diag.128.4)
  mean.Bw.mat.2b.thomas.pi.5 <- mean.Bw.mat.2b.thomas.pi.5 + Bw.pi.mat.5
  mean.Bw.diag.2b.thomas.pi.5 <- mean.Bw.diag.2b.thomas.pi.5 + Bw.pi.diag.5
  MISE.mat.2b.thomas.pi.5 <- MISE.mat.2b.thomas.pi.5 + F.computeISE(im.2b.128, dens.pi.128.5)
  MISE.diag.2b.thomas.pi.5 <- MISE.diag.2b.thomas.pi.5 + F.computeISE(im.2b.128, dens.pi.diag.128.5)
  mean.Bw.mat.2b.thomas.pi.10 <- mean.Bw.mat.2b.thomas.pi.10 + Bw.pi.mat.10
  mean.Bw.diag.2b.thomas.pi.10 <- mean.Bw.diag.2b.thomas.pi.10 + Bw.pi.diag.10
  MISE.mat.2b.thomas.pi.10 <- MISE.mat.2b.thomas.pi.10 + F.computeISE(im.2b.128, dens.pi.128.10)
  MISE.diag.2b.thomas.pi.10 <- MISE.diag.2b.thomas.pi.10 + F.computeISE(im.2b.128, dens.pi.diag.128.10)
  time.2b.thomas.pi <- time.2b.thomas.pi + Sys.time() - st
  
  # ise = F.computeISE(im.2b.128, dens.128)
  # MISE.mat.2b <- MISE.mat.2b + ise
  # ise.diag = F.computeISE(im.2b.128, dens.diag.128)
  # MISE.diag.2b <- MISE.diag.2b + ise.diag
  # cat(ise, ise.diag, "\n")
}

time.2b.thomas.pi <- time.2b.thomas.pi/N
MISE.mat.2b.thomas.pi <- MISE.mat.2b.thomas.pi/N
MISE.diag.2b.thomas.pi <- MISE.diag.2b.thomas.pi/N
mean.Bw.mat.2b.thomas.pi <- mean.Bw.mat.2b.thomas.pi/N
mean.Bw.diag.2b.thomas.pi <- mean.Bw.diag.2b.thomas.pi/N
MISE.mat.2b.thomas.pi.1 <- MISE.mat.2b.thomas.pi.1/N
MISE.diag.2b.thomas.pi.1 <- MISE.diag.2b.thomas.pi.1/N
mean.Bw.mat.2b.thomas.pi.1 <- mean.Bw.mat.2b.thomas.pi.1/N
mean.Bw.diag.2b.thomas.pi.1 <- mean.Bw.diag.2b.thomas.pi.1/N
MISE.mat.2b.thomas.pi.2 <- MISE.mat.2b.thomas.pi.2/N
MISE.diag.2b.thomas.pi.2 <- MISE.diag.2b.thomas.pi.2/N
mean.Bw.mat.2b.thomas.pi.2 <- mean.Bw.mat.2b.thomas.pi.2/N
mean.Bw.diag.2b.thomas.pi.2 <- mean.Bw.diag.2b.thomas.pi.2/N
MISE.mat.2b.thomas.pi.3 <- MISE.mat.2b.thomas.pi.3/N
MISE.diag.2b.thomas.pi.3 <- MISE.diag.2b.thomas.pi.3/N
mean.Bw.mat.2b.thomas.pi.3 <- mean.Bw.mat.2b.thomas.pi.3/N
mean.Bw.diag.2b.thomas.pi.3 <- mean.Bw.diag.2b.thomas.pi.3/N
MISE.mat.2b.thomas.pi.4 <- MISE.mat.2b.thomas.pi.4/N
MISE.diag.2b.thomas.pi.4 <- MISE.diag.2b.thomas.pi.4/N
mean.Bw.mat.2b.thomas.pi.4 <- mean.Bw.mat.2b.thomas.pi.4/N
mean.Bw.diag.2b.thomas.pi.4 <- mean.Bw.diag.2b.thomas.pi.4/N
MISE.mat.2b.thomas.pi.5 <- MISE.mat.2b.thomas.pi.5/N
MISE.diag.2b.thomas.pi.5 <- MISE.diag.2b.thomas.pi.5/N
mean.Bw.mat.2b.thomas.pi.5 <- mean.Bw.mat.2b.thomas.pi.5/N
mean.Bw.diag.2b.thomas.pi.5 <- mean.Bw.diag.2b.thomas.pi.5/N
MISE.mat.2b.thomas.pi.10 <- MISE.mat.2b.thomas.pi.10/N
MISE.diag.2b.thomas.pi.10 <- MISE.diag.2b.thomas.pi.10/N
mean.Bw.mat.2b.thomas.pi.10 <- mean.Bw.mat.2b.thomas.pi.10/N
mean.Bw.diag.2b.thomas.pi.10 <- mean.Bw.diag.2b.thomas.pi.10/N



# 3b, scv

N <- 100

set.seed(887330)

MISE.mat.3b <- 0
MISE.diag.3b <- 0
MISE.mat.3b.1 <- 0
MISE.diag.3b.1 <- 0
MISE.mat.3b.2 <- 0
MISE.diag.3b.2 <- 0
MISE.mat.3b.3 <- 0
MISE.diag.3b.3 <- 0
MISE.mat.3b.4 <- 0
MISE.diag.3b.4 <- 0
MISE.mat.3b.5 <- 0
MISE.diag.3b.5 <- 0
MISE.mat.3b.10 <- 0
MISE.diag.3b.10 <- 0
mean.Bw.mat.3b <- diag(0,2,2)
mean.Bw.diag.3b <- diag(0,2,2)
mean.Bw.mat.3b.1 <- diag(0,2,2)
mean.Bw.diag.3b.1 <- diag(0,2,2)
mean.Bw.mat.3b.2 <- diag(0,2,2)
mean.Bw.diag.3b.2 <- diag(0,2,2)
mean.Bw.mat.3b.3 <- diag(0,2,2)
mean.Bw.diag.3b.3 <- diag(0,2,2)
mean.Bw.mat.3b.4 <- diag(0,2,2)
mean.Bw.diag.3b.4 <- diag(0,2,2)
mean.Bw.mat.3b.5 <- diag(0,2,2)
mean.Bw.diag.3b.5 <- diag(0,2,2)
mean.Bw.mat.3b.10 <- diag(0,2,2)
mean.Bw.diag.3b.10 <- diag(0,2,2)
time.3b <- 0

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.3b.128 <- as.im(lambda3b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.3b <- rpoispp(im.3b.128)
  pp.3b.df <- data.frame(x=pp.3b$x, y=pp.3b$y)
  Bw.mat <- Hscv(pp.3b.df)
  Bw.diag <- Hscv.diag(pp.3b.df)
  Bw.mat.1 <- Hscv.dir(pp.3b.df, alpha=1)
  Bw.diag.1 <- Hscv.diag.dir(pp.3b.df, alpha=1)
  Bw.mat.2 <- Hscv.dir(pp.3b.df, alpha=2)
  Bw.diag.2 <- Hscv.diag.dir(pp.3b.df, alpha=2)
  Bw.mat.3 <- Hscv.dir(pp.3b.df, alpha=3)
  Bw.diag.3 <- Hscv.diag.dir(pp.3b.df, alpha=3)
  Bw.mat.4 <- Hscv.dir(pp.3b.df, alpha=4)
  Bw.diag.4 <- Hscv.diag.dir(pp.3b.df, alpha=4)
  Bw.mat.5 <- Hscv.dir(pp.3b.df, alpha=5)
  Bw.diag.5 <- Hscv.diag.dir(pp.3b.df, alpha=5)
  Bw.mat.10 <- Hscv.dir(pp.3b.df, alpha=10)
  Bw.diag.10 <- Hscv.diag.dir(pp.3b.df, alpha=10)
  
  dens.128 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.mat)
  dens.diag.128 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.diag)
  dens.128.1 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.mat.1)
  dens.diag.128.1 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.diag.1)
  dens.128.2 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.mat.2)
  dens.diag.128.2 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.diag.2)
  dens.128.3 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.mat.3)
  dens.diag.128.3 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.diag.3)
  dens.128.4 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.mat.4)
  dens.diag.128.4 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.diag.4)
  dens.128.5 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.mat.5)
  dens.diag.128.5 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.diag.5)
  dens.128.10 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.mat.10)
  dens.diag.128.10 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.diag.10)
  
  mean.Bw.mat.3b <- mean.Bw.mat.3b + Bw.mat
  mean.Bw.diag.3b <- mean.Bw.diag.3b + Bw.diag
  MISE.mat.3b <- MISE.mat.3b + F.computeISE(im.3b.128, dens.128)
  MISE.diag.3b <- MISE.diag.3b + F.computeISE(im.3b.128, dens.diag.128)
  mean.Bw.mat.3b.1 <- mean.Bw.mat.3b.1 + Bw.mat.1
  mean.Bw.diag.3b.1 <- mean.Bw.diag.3b.1 + Bw.diag.1
  MISE.mat.3b.1 <- MISE.mat.3b.1 + F.computeISE(im.3b.128, dens.128.1)
  MISE.diag.3b.1 <- MISE.diag.3b.1 + F.computeISE(im.3b.128, dens.diag.128.1)
  mean.Bw.mat.3b.2 <- mean.Bw.mat.3b.2 + Bw.mat.2
  mean.Bw.diag.3b.2 <- mean.Bw.diag.3b.2 + Bw.diag.2
  MISE.mat.3b.2 <- MISE.mat.3b.2 + F.computeISE(im.3b.128, dens.128.2)
  MISE.diag.3b.2 <- MISE.diag.3b.2 + F.computeISE(im.3b.128, dens.diag.128.2)
  mean.Bw.mat.3b.3 <- mean.Bw.mat.3b.3 + Bw.mat.3
  mean.Bw.diag.3b.3 <- mean.Bw.diag.3b.3 + Bw.diag.3
  MISE.mat.3b.3 <- MISE.mat.3b.3 + F.computeISE(im.3b.128, dens.128.3)
  MISE.diag.3b.3 <- MISE.diag.3b.3 + F.computeISE(im.3b.128, dens.diag.128.3)
  mean.Bw.mat.3b.4 <- mean.Bw.mat.3b.4 + Bw.mat.4
  mean.Bw.diag.3b.4 <- mean.Bw.diag.3b.4 + Bw.diag.4
  MISE.mat.3b.4 <- MISE.mat.3b.4 + F.computeISE(im.3b.128, dens.128.4)
  MISE.diag.3b.4 <- MISE.diag.3b.4 + F.computeISE(im.3b.128, dens.diag.128.4)
  mean.Bw.mat.3b.5 <- mean.Bw.mat.3b.5 + Bw.mat.5
  mean.Bw.diag.3b.5 <- mean.Bw.diag.3b.5 + Bw.diag.5
  MISE.mat.3b.5 <- MISE.mat.3b.5 + F.computeISE(im.3b.128, dens.128.5)
  MISE.diag.3b.5 <- MISE.diag.3b.5 + F.computeISE(im.3b.128, dens.diag.128.5)
  mean.Bw.mat.3b.10 <- mean.Bw.mat.3b.10 + Bw.mat.10
  mean.Bw.diag.3b.10 <- mean.Bw.diag.3b.10 + Bw.diag.10
  MISE.mat.3b.10 <- MISE.mat.3b.10 + F.computeISE(im.3b.128, dens.128.10)
  MISE.diag.3b.10 <- MISE.diag.3b.10 + F.computeISE(im.3b.128, dens.diag.128.10)
  time.3b <- time.3b + Sys.time() - st
}

time.3b <- time.3b/N
MISE.mat.3b <- MISE.mat.3b/N
MISE.diag.3b <- MISE.diag.3b/N
mean.Bw.mat.3b <- mean.Bw.mat.3b/N
mean.Bw.diag.3b <- mean.Bw.diag.3b/N
MISE.mat.3b.1 <- MISE.mat.3b.1/N
MISE.diag.3b.1 <- MISE.diag.3b.1/N
mean.Bw.mat.3b.1 <- mean.Bw.mat.3b.1/N
mean.Bw.diag.3b.1 <- mean.Bw.diag.3b.1/N
MISE.mat.3b.2 <- MISE.mat.3b.2/N
MISE.diag.3b.2 <- MISE.diag.3b.2/N
mean.Bw.mat.3b.2 <- mean.Bw.mat.3b.2/N
mean.Bw.diag.3b.2 <- mean.Bw.diag.3b.2/N
MISE.mat.3b.3 <- MISE.mat.3b.3/N
MISE.diag.3b.3 <- MISE.diag.3b.3/N
mean.Bw.mat.3b.3 <- mean.Bw.mat.3b.3/N
mean.Bw.diag.3b.3 <- mean.Bw.diag.3b.3/N
MISE.mat.3b.4 <- MISE.mat.3b.4/N
MISE.diag.3b.4 <- MISE.diag.3b.4/N
mean.Bw.mat.3b.4 <- mean.Bw.mat.3b.4/N
mean.Bw.diag.3b.4 <- mean.Bw.diag.3b.4/N
MISE.mat.3b.5 <- MISE.mat.3b.5/N
MISE.diag.3b.5 <- MISE.diag.3b.5/N
mean.Bw.mat.3b.5 <- mean.Bw.mat.3b.5/N
mean.Bw.diag.3b.5 <- mean.Bw.diag.3b.5/N
MISE.mat.3b.10 <- MISE.mat.3b.10/N
MISE.diag.3b.10 <- MISE.diag.3b.10/N
mean.Bw.mat.3b.10 <- mean.Bw.mat.3b.10/N
mean.Bw.diag.3b.10 <- mean.Bw.diag.3b.10/N


# 3b, pi

N <- 100

set.seed(887330)

MISE.mat.3b.pi <- 0
MISE.diag.3b.pi <- 0
MISE.mat.3b.pi.1 <- 0
MISE.diag.3b.pi.1 <- 0
MISE.mat.3b.pi.2 <- 0
MISE.diag.3b.pi.2 <- 0
MISE.mat.3b.pi.3 <- 0
MISE.diag.3b.pi.3 <- 0
MISE.mat.3b.pi.4 <- 0
MISE.diag.3b.pi.4 <- 0
MISE.mat.3b.pi.5 <- 0
MISE.diag.3b.pi.5 <- 0
MISE.mat.3b.pi.10 <- 0
MISE.diag.3b.pi.10 <- 0
mean.Bw.mat.3b.pi <- diag(0,2,2)
mean.Bw.diag.3b.pi <- diag(0,2,2)
mean.Bw.mat.3b.pi.1 <- diag(0,2,2)
mean.Bw.diag.3b.pi.1 <- diag(0,2,2)
mean.Bw.mat.3b.pi.2 <- diag(0,2,2)
mean.Bw.diag.3b.pi.2 <- diag(0,2,2)
mean.Bw.mat.3b.pi.3 <- diag(0,2,2)
mean.Bw.diag.3b.pi.3 <- diag(0,2,2)
mean.Bw.mat.3b.pi.4 <- diag(0,2,2)
mean.Bw.diag.3b.pi.4 <- diag(0,2,2)
mean.Bw.mat.3b.pi.5 <- diag(0,2,2)
mean.Bw.diag.3b.pi.5 <- diag(0,2,2)
mean.Bw.mat.3b.pi.10 <- diag(0,2,2)
mean.Bw.diag.3b.pi.10 <- diag(0,2,2)
time.3b.pi <- 0

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.3b.128 <- as.im(lambda3b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.3b <- rpoispp(im.3b.128)
  pp.3b.df <- data.frame(x=pp.3b$x, y=pp.3b$y)
  Bw.mat <- Hpi(pp.3b.df)
  Bw.diag <- Hpi.diag(pp.3b.df)
  Bw.pi.mat.1 <- Hpi.dir(pp.3b.df, alpha=1)
  Bw.pi.diag.1 <- Hpi.diag.dir(pp.3b.df, alpha=1)
  Bw.pi.mat.2 <- Hpi.dir(pp.3b.df, alpha=2)
  Bw.pi.diag.2 <- Hpi.diag.dir(pp.3b.df, alpha=2)
  Bw.pi.mat.3 <- Hpi.dir(pp.3b.df, alpha=3)
  Bw.pi.diag.3 <- Hpi.diag.dir(pp.3b.df, alpha=3)
  Bw.pi.mat.4 <- Hpi.dir(pp.3b.df, alpha=4)
  Bw.pi.diag.4 <- Hpi.diag.dir(pp.3b.df, alpha=4)
  Bw.pi.mat.5 <- Hpi.dir(pp.3b.df, alpha=5)
  Bw.pi.diag.5 <- Hpi.diag.dir(pp.3b.df, alpha=5)
  Bw.pi.mat.10 <- Hpi.dir(pp.3b.df, alpha=10)
  Bw.pi.diag.10 <- Hpi.diag.dir(pp.3b.df, alpha=10)
  
  dens.128 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.mat)
  dens.diag.128 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.diag)
  dens.pi.128.1 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.pi.mat.1)
  dens.pi.diag.128.1 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.pi.diag.1)
  dens.pi.128.2 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.pi.mat.2)
  dens.pi.diag.128.2 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.pi.diag.2)
  dens.pi.128.3 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.pi.mat.3)
  dens.pi.diag.128.3 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.pi.diag.3)
  dens.pi.128.4 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.pi.mat.4)
  dens.pi.diag.128.4 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.pi.diag.4)
  dens.pi.128.5 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.pi.mat.5)
  dens.pi.diag.128.5 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.pi.diag.5)
  dens.pi.128.10 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.pi.mat.10)
  dens.pi.diag.128.10 <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.pi.diag.10)
  
  
  mean.Bw.mat.3b.pi <- mean.Bw.mat.3b.pi + Bw.mat
  mean.Bw.diag.3b.pi <- mean.Bw.diag.3b.pi + Bw.diag
  MISE.mat.3b.pi <- MISE.mat.3b.pi + F.computeISE(im.3b.128, dens.128)
  MISE.diag.3b.pi <- MISE.diag.3b.pi + F.computeISE(im.3b.128, dens.diag.128)
  mean.Bw.mat.3b.pi.1 <- mean.Bw.mat.3b.pi.1 + Bw.pi.mat.1
  mean.Bw.diag.3b.pi.1 <- mean.Bw.diag.3b.pi.1 + Bw.pi.diag.1
  MISE.mat.3b.pi.1 <- MISE.mat.3b.pi.1 + F.computeISE(im.3b.128, dens.pi.128.1)
  MISE.diag.3b.pi.1 <- MISE.diag.3b.pi.1 + F.computeISE(im.3b.128, dens.pi.diag.128.1)
  mean.Bw.mat.3b.pi.2 <- mean.Bw.mat.3b.pi.2 + Bw.pi.mat.2
  mean.Bw.diag.3b.pi.2 <- mean.Bw.diag.3b.pi.2 + Bw.pi.diag.2
  MISE.mat.3b.pi.2 <- MISE.mat.3b.pi.2 + F.computeISE(im.3b.128, dens.pi.128.2)
  MISE.diag.3b.pi.2 <- MISE.diag.3b.pi.2 + F.computeISE(im.3b.128, dens.pi.diag.128.2)
  mean.Bw.mat.3b.pi.3 <- mean.Bw.mat.3b.pi.3 + Bw.pi.mat.3
  mean.Bw.diag.3b.pi.3 <- mean.Bw.diag.3b.pi.3 + Bw.pi.diag.3
  MISE.mat.3b.pi.3 <- MISE.mat.3b.pi.3 + F.computeISE(im.3b.128, dens.pi.128.3)
  MISE.diag.3b.pi.3 <- MISE.diag.3b.pi.3 + F.computeISE(im.3b.128, dens.pi.diag.128.3)
  mean.Bw.mat.3b.pi.4 <- mean.Bw.mat.3b.pi.4 + Bw.pi.mat.4
  mean.Bw.diag.3b.pi.4 <- mean.Bw.diag.3b.pi.4 + Bw.pi.diag.4
  MISE.mat.3b.pi.4 <- MISE.mat.3b.pi.4 + F.computeISE(im.3b.128, dens.pi.128.4)
  MISE.diag.3b.pi.4 <- MISE.diag.3b.pi.4 + F.computeISE(im.3b.128, dens.pi.diag.128.4)
  mean.Bw.mat.3b.pi.5 <- mean.Bw.mat.3b.pi.5 + Bw.pi.mat.5
  mean.Bw.diag.3b.pi.5 <- mean.Bw.diag.3b.pi.5 + Bw.pi.diag.5
  MISE.mat.3b.pi.5 <- MISE.mat.3b.pi.5 + F.computeISE(im.3b.128, dens.pi.128.5)
  MISE.diag.3b.pi.5 <- MISE.diag.3b.pi.5 + F.computeISE(im.3b.128, dens.pi.diag.128.5)
  mean.Bw.mat.3b.pi.10 <- mean.Bw.mat.3b.pi.10 + Bw.pi.mat.10
  mean.Bw.diag.3b.pi.10 <- mean.Bw.diag.3b.pi.10 + Bw.pi.diag.10
  MISE.mat.3b.pi.10 <- MISE.mat.3b.pi.10 + F.computeISE(im.3b.128, dens.pi.128.10)
  MISE.diag.3b.pi.10 <- MISE.diag.3b.pi.10 + F.computeISE(im.3b.128, dens.pi.diag.128.10)
  time.3b.pi <- time.3b.pi + Sys.time() - st
  
  # ise = F.computeISE(im.3b.128, dens.128)
  # MISE.mat.3b <- MISE.mat.3b + ise
  # ise.diag = F.computeISE(im.3b.128, dens.diag.128)
  # MISE.diag.3b <- MISE.diag.3b + ise.diag
  # cat(ise, ise.diag, "\n")
}

time.3b.pi <- time.3b.pi/N
MISE.mat.3b.pi <- MISE.mat.3b.pi/N
MISE.diag.3b.pi <- MISE.diag.3b.pi/N
mean.Bw.mat.3b.pi <- mean.Bw.mat.3b.pi/N
mean.Bw.diag.3b.pi <- mean.Bw.diag.3b.pi/N
MISE.mat.3b.pi.1 <- MISE.mat.3b.pi.1/N
MISE.diag.3b.pi.1 <- MISE.diag.3b.pi.1/N
mean.Bw.mat.3b.pi.1 <- mean.Bw.mat.3b.pi.1/N
mean.Bw.diag.3b.pi.1 <- mean.Bw.diag.3b.pi.1/N
MISE.mat.3b.pi.2 <- MISE.mat.3b.pi.2/N
MISE.diag.3b.pi.2 <- MISE.diag.3b.pi.2/N
mean.Bw.mat.3b.pi.2 <- mean.Bw.mat.3b.pi.2/N
mean.Bw.diag.3b.pi.2 <- mean.Bw.diag.3b.pi.2/N
MISE.mat.3b.pi.3 <- MISE.mat.3b.pi.3/N
MISE.diag.3b.pi.3 <- MISE.diag.3b.pi.3/N
mean.Bw.mat.3b.pi.3 <- mean.Bw.mat.3b.pi.3/N
mean.Bw.diag.3b.pi.3 <- mean.Bw.diag.3b.pi.3/N
MISE.mat.3b.pi.4 <- MISE.mat.3b.pi.4/N
MISE.diag.3b.pi.4 <- MISE.diag.3b.pi.4/N
mean.Bw.mat.3b.pi.4 <- mean.Bw.mat.3b.pi.4/N
mean.Bw.diag.3b.pi.4 <- mean.Bw.diag.3b.pi.4/N
MISE.mat.3b.pi.5 <- MISE.mat.3b.pi.5/N
MISE.diag.3b.pi.5 <- MISE.diag.3b.pi.5/N
mean.Bw.mat.3b.pi.5 <- mean.Bw.mat.3b.pi.5/N
mean.Bw.diag.3b.pi.5 <- mean.Bw.diag.3b.pi.5/N
MISE.mat.3b.pi.10 <- MISE.mat.3b.pi.10/N
MISE.diag.3b.pi.10 <- MISE.diag.3b.pi.10/N
mean.Bw.mat.3b.pi.10 <- mean.Bw.mat.3b.pi.10/N
mean.Bw.diag.3b.pi.10 <- mean.Bw.diag.3b.pi.10/N



# 4b, scv

N <- 100

set.seed(887330)

MISE.mat.4b <- 0
MISE.diag.4b <- 0
MISE.mat.4b.1 <- 0
MISE.diag.4b.1 <- 0
MISE.mat.4b.2 <- 0
MISE.diag.4b.2 <- 0
MISE.mat.4b.3 <- 0
MISE.diag.4b.3 <- 0
MISE.mat.4b.4 <- 0
MISE.diag.4b.4 <- 0
MISE.mat.4b.5 <- 0
MISE.diag.4b.5 <- 0
MISE.mat.4b.10 <- 0
MISE.diag.4b.10 <- 0
mean.Bw.mat.4b <- diag(0,2,2)
mean.Bw.diag.4b <- diag(0,2,2)
mean.Bw.mat.4b.1 <- diag(0,2,2)
mean.Bw.diag.4b.1 <- diag(0,2,2)
mean.Bw.mat.4b.2 <- diag(0,2,2)
mean.Bw.diag.4b.2 <- diag(0,2,2)
mean.Bw.mat.4b.3 <- diag(0,2,2)
mean.Bw.diag.4b.3 <- diag(0,2,2)
mean.Bw.mat.4b.4 <- diag(0,2,2)
mean.Bw.diag.4b.4 <- diag(0,2,2)
mean.Bw.mat.4b.5 <- diag(0,2,2)
mean.Bw.diag.4b.5 <- diag(0,2,2)
mean.Bw.mat.4b.10 <- diag(0,2,2)
mean.Bw.diag.4b.10 <- diag(0,2,2)

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.4b.128 <- as.im(lambda4b(my.sim.RF), W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.4b <- rpoispp(im.4b.128)
  pp.4b.df <- data.frame(x=pp.4b$x, y=pp.4b$y)
  Bw.mat <- Hscv(pp.4b.df)
  Bw.diag <- Hscv.diag(pp.4b.df)
  Bw.mat.1 <- Hscv.dir(pp.4b.df, alpha=1)
  Bw.diag.1 <- Hscv.diag.dir(pp.4b.df, alpha=1)
  Bw.mat.2 <- Hscv.dir(pp.4b.df, alpha=2)
  Bw.diag.2 <- Hscv.diag.dir(pp.4b.df, alpha=2)
  Bw.mat.3 <- Hscv.dir(pp.4b.df, alpha=3)
  Bw.diag.3 <- Hscv.diag.dir(pp.4b.df, alpha=3)
  Bw.mat.4 <- Hscv.dir(pp.4b.df, alpha=4)
  Bw.diag.4 <- Hscv.diag.dir(pp.4b.df, alpha=4)
  Bw.mat.5 <- Hscv.dir(pp.4b.df, alpha=5)
  Bw.diag.5 <- Hscv.diag.dir(pp.4b.df, alpha=5)
  Bw.mat.10 <- Hscv.dir(pp.4b.df, alpha=10)
  Bw.diag.10 <- Hscv.diag.dir(pp.4b.df, alpha=10)
  
  dens.128 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.mat)
  dens.diag.128 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.diag)
  dens.128.1 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.mat.1)
  dens.diag.128.1 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.diag.1)
  dens.128.2 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.mat.2)
  dens.diag.128.2 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.diag.2)
  dens.128.3 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.mat.3)
  dens.diag.128.3 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.diag.3)
  dens.128.4 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.mat.4)
  dens.diag.128.4 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.diag.4)
  dens.128.5 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.mat.5)
  dens.diag.128.5 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.diag.5)
  dens.128.10 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.mat.10)
  dens.diag.128.10 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.diag.10)
  
  mean.Bw.mat.4b <- mean.Bw.mat.4b + Bw.mat
  mean.Bw.diag.4b <- mean.Bw.diag.4b + Bw.diag
  MISE.mat.4b <- MISE.mat.4b + F.computeISE(im.4b.128, dens.128)
  MISE.diag.4b <- MISE.diag.4b + F.computeISE(im.4b.128, dens.diag.128)
  mean.Bw.mat.4b.1 <- mean.Bw.mat.4b.1 + Bw.mat.1
  mean.Bw.diag.4b.1 <- mean.Bw.diag.4b.1 + Bw.diag.1
  MISE.mat.4b.1 <- MISE.mat.4b.1 + F.computeISE(im.4b.128, dens.128.1)
  MISE.diag.4b.1 <- MISE.diag.4b.1 + F.computeISE(im.4b.128, dens.diag.128.1)
  mean.Bw.mat.4b.2 <- mean.Bw.mat.4b.2 + Bw.mat.2
  mean.Bw.diag.4b.2 <- mean.Bw.diag.4b.2 + Bw.diag.2
  MISE.mat.4b.2 <- MISE.mat.4b.2 + F.computeISE(im.4b.128, dens.128.2)
  MISE.diag.4b.2 <- MISE.diag.4b.2 + F.computeISE(im.4b.128, dens.diag.128.2)
  mean.Bw.mat.4b.3 <- mean.Bw.mat.4b.3 + Bw.mat.3
  mean.Bw.diag.4b.3 <- mean.Bw.diag.4b.3 + Bw.diag.3
  MISE.mat.4b.3 <- MISE.mat.4b.3 + F.computeISE(im.4b.128, dens.128.3)
  MISE.diag.4b.3 <- MISE.diag.4b.3 + F.computeISE(im.4b.128, dens.diag.128.3)
  mean.Bw.mat.4b.4 <- mean.Bw.mat.4b.4 + Bw.mat.4
  mean.Bw.diag.4b.4 <- mean.Bw.diag.4b.4 + Bw.diag.4
  MISE.mat.4b.4 <- MISE.mat.4b.4 + F.computeISE(im.4b.128, dens.128.4)
  MISE.diag.4b.4 <- MISE.diag.4b.4 + F.computeISE(im.4b.128, dens.diag.128.4)
  mean.Bw.mat.4b.5 <- mean.Bw.mat.4b.5 + Bw.mat.5
  mean.Bw.diag.4b.5 <- mean.Bw.diag.4b.5 + Bw.diag.5
  MISE.mat.4b.5 <- MISE.mat.4b.5 + F.computeISE(im.4b.128, dens.128.5)
  MISE.diag.4b.5 <- MISE.diag.4b.5 + F.computeISE(im.4b.128, dens.diag.128.5)
  mean.Bw.mat.4b.10 <- mean.Bw.mat.4b.10 + Bw.mat.10
  mean.Bw.diag.4b.10 <- mean.Bw.diag.4b.10 + Bw.diag.10
  MISE.mat.4b.10 <- MISE.mat.4b.10 + F.computeISE(im.4b.128, dens.128.10)
  MISE.diag.4b.10 <- MISE.diag.4b.10 + F.computeISE(im.4b.128, dens.diag.128.10)
  time.4b <- time.4b + Sys.time() - st
  # ise = F.computeISE(im.4b.128, dens.128)
  # MISE.mat.4b <- MISE.mat.4b + ise
  # ise.diag = F.computeISE(im.4b.128, dens.diag.128)
  # MISE.diag.4b <- MISE.diag.4b + ise.diag
  # cat(ise, ise.diag, "\n")
}

time.4b <- time.4b/N
MISE.mat.4b <- MISE.mat.4b/N
MISE.diag.4b <- MISE.diag.4b/N
mean.Bw.mat.4b <- mean.Bw.mat.4b/N
mean.Bw.diag.4b <- mean.Bw.diag.4b/N
MISE.mat.4b.1 <- MISE.mat.4b.1/N
MISE.diag.4b.1 <- MISE.diag.4b.1/N
mean.Bw.mat.4b.1 <- mean.Bw.mat.4b.1/N
mean.Bw.diag.4b.1 <- mean.Bw.diag.4b.1/N
MISE.mat.4b.2 <- MISE.mat.4b.2/N
MISE.diag.4b.2 <- MISE.diag.4b.2/N
mean.Bw.mat.4b.2 <- mean.Bw.mat.4b.2/N
mean.Bw.diag.4b.2 <- mean.Bw.diag.4b.2/N
MISE.mat.4b.3 <- MISE.mat.4b.3/N
MISE.diag.4b.3 <- MISE.diag.4b.3/N
mean.Bw.mat.4b.3 <- mean.Bw.mat.4b.3/N
mean.Bw.diag.4b.3 <- mean.Bw.diag.4b.3/N
MISE.mat.4b.4 <- MISE.mat.4b.4/N
MISE.diag.4b.4 <- MISE.diag.4b.4/N
mean.Bw.mat.4b.4 <- mean.Bw.mat.4b.4/N
mean.Bw.diag.4b.4 <- mean.Bw.diag.4b.4/N
MISE.mat.4b.5 <- MISE.mat.4b.5/N
MISE.diag.4b.5 <- MISE.diag.4b.5/N
mean.Bw.mat.4b.5 <- mean.Bw.mat.4b.5/N
mean.Bw.diag.4b.5 <- mean.Bw.diag.4b.5/N
MISE.mat.4b.10 <- MISE.mat.4b.10/N
MISE.diag.4b.10 <- MISE.diag.4b.10/N
mean.Bw.mat.4b.10 <- mean.Bw.mat.4b.10/N
mean.Bw.diag.4b.10 <- mean.Bw.diag.4b.10/N




# 4b, pi

N <- 100

set.seed(887330)

MISE.mat.4b.pi <- 0
MISE.diag.4b.pi <- 0
MISE.mat.4b.pi.1 <- 0
MISE.diag.4b.pi.1 <- 0
MISE.mat.4b.pi.2 <- 0
MISE.diag.4b.pi.2 <- 0
MISE.mat.4b.pi.3 <- 0
MISE.diag.4b.pi.3 <- 0
MISE.mat.4b.pi.4 <- 0
MISE.diag.4b.pi.4 <- 0
MISE.mat.4b.pi.5 <- 0
MISE.diag.4b.pi.5 <- 0
MISE.mat.4b.pi.10 <- 0
MISE.diag.4b.pi.10 <- 0
mean.Bw.mat.4b.pi <- diag(0,2,2)
mean.Bw.diag.4b.pi <- diag(0,2,2)
mean.Bw.mat.4b.pi.1 <- diag(0,2,2)
mean.Bw.diag.4b.pi.1 <- diag(0,2,2)
mean.Bw.mat.4b.pi.2 <- diag(0,2,2)
mean.Bw.diag.4b.pi.2 <- diag(0,2,2)
mean.Bw.mat.4b.pi.3 <- diag(0,2,2)
mean.Bw.diag.4b.pi.3 <- diag(0,2,2)
mean.Bw.mat.4b.pi.4 <- diag(0,2,2)
mean.Bw.diag.4b.pi.4 <- diag(0,2,2)
mean.Bw.mat.4b.pi.5 <- diag(0,2,2)
mean.Bw.diag.4b.pi.5 <- diag(0,2,2)
mean.Bw.mat.4b.pi.10 <- diag(0,2,2)
mean.Bw.diag.4b.pi.10 <- diag(0,2,2)
time.4b.pi <- 0

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.4b.128 <- as.im(lambda4b(my.sim.RF), W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.4b <- rpoispp(im.4b.128)
  pp.4b.df <- data.frame(x=pp.4b$x, y=pp.4b$y)
  Bw.pi.mat <- Hpi(pp.4b.df)
  Bw.pi.diag <- Hpi.diag(pp.4b.df)
  Bw.pi.mat.1 <- Hpi.dir(pp.4b.df, alpha=1)
  Bw.pi.diag.1 <- Hpi.diag.dir(pp.4b.df, alpha=1)
  Bw.pi.mat.2 <- Hpi.dir(pp.4b.df, alpha=2)
  Bw.pi.diag.2 <- Hpi.diag.dir(pp.4b.df, alpha=2)
  Bw.pi.mat.3 <- Hpi.dir(pp.4b.df, alpha=3)
  Bw.pi.diag.3 <- Hpi.diag.dir(pp.4b.df, alpha=3)
  Bw.pi.mat.4 <- Hpi.dir(pp.4b.df, alpha=4)
  Bw.pi.diag.4 <- Hpi.diag.dir(pp.4b.df, alpha=4)
  Bw.pi.mat.5 <- Hpi.dir(pp.4b.df, alpha=5)
  Bw.pi.diag.5 <- Hpi.diag.dir(pp.4b.df, alpha=5)
  Bw.pi.mat.10 <- Hpi.dir(pp.4b.df, alpha=10)
  Bw.pi.diag.10 <- Hpi.diag.dir(pp.4b.df, alpha=10)
  
  dens.pi.128 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.pi.mat)
  dens.pi.diag.128 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.pi.diag)
  dens.pi.128.1 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.pi.mat.1)
  dens.pi.diag.128.1 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.pi.diag.1)
  dens.pi.128.2 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.pi.mat.2)
  dens.pi.diag.128.2 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.pi.diag.2)
  dens.pi.128.3 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.pi.mat.3)
  dens.pi.diag.128.3 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.pi.diag.3)
  dens.pi.128.4 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.pi.mat.4)
  dens.pi.diag.128.4 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.pi.diag.4)
  dens.pi.128.5 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.pi.mat.5)
  dens.pi.diag.128.5 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.pi.diag.5)
  dens.pi.128.10 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.pi.mat.10)
  dens.pi.diag.128.10 <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.pi.diag.10)
  
  mean.Bw.mat.4b.pi <- mean.Bw.mat.4b.pi + Bw.mat
  mean.Bw.diag.4b.pi <- mean.Bw.diag.4b.pi + Bw.diag
  MISE.mat.4b.pi <- MISE.mat.4b.pi + F.computeISE(im.4b.128, dens.pi.128)
  MISE.diag.4b.pi <- MISE.diag.4b.pi + F.computeISE(im.4b.128, dens.pi.diag.128)
  mean.Bw.mat.4b.pi.1 <- mean.Bw.mat.4b.pi.1 + Bw.mat.1
  mean.Bw.diag.4b.pi.1 <- mean.Bw.diag.4b.pi.1 + Bw.diag.1
  MISE.mat.4b.pi.1 <- MISE.mat.4b.pi.1 + F.computeISE(im.4b.128, dens.pi.128.1)
  MISE.diag.4b.pi.1 <- MISE.diag.4b.pi.1 + F.computeISE(im.4b.128, dens.pi.diag.128.1)
  mean.Bw.mat.4b.pi.2 <- mean.Bw.mat.4b.pi.2 + Bw.mat.2
  mean.Bw.diag.4b.pi.2 <- mean.Bw.diag.4b.pi.2 + Bw.diag.2
  MISE.mat.4b.pi.2 <- MISE.mat.4b.pi.2 + F.computeISE(im.4b.128, dens.pi.128.2)
  MISE.diag.4b.pi.2 <- MISE.diag.4b.pi.2 + F.computeISE(im.4b.128, dens.pi.diag.128.2)
  mean.Bw.mat.4b.pi.3 <- mean.Bw.mat.4b.pi.3 + Bw.mat.3
  mean.Bw.diag.4b.pi.3 <- mean.Bw.diag.4b.pi.3 + Bw.diag.3
  MISE.mat.4b.pi.3 <- MISE.mat.4b.pi.3 + F.computeISE(im.4b.128, dens.pi.128.3)
  MISE.diag.4b.pi.3 <- MISE.diag.4b.pi.3 + F.computeISE(im.4b.128, dens.pi.diag.128.3)
  mean.Bw.mat.4b.pi.4 <- mean.Bw.mat.4b.pi.4 + Bw.mat.4
  mean.Bw.diag.4b.pi.4 <- mean.Bw.diag.4b.pi.4 + Bw.diag.4
  MISE.mat.4b.pi.4 <- MISE.mat.4b.pi.4 + F.computeISE(im.4b.128, dens.pi.128.4)
  MISE.diag.4b.pi.4 <- MISE.diag.4b.pi.4 + F.computeISE(im.4b.128, dens.pi.diag.128.4)
  mean.Bw.mat.4b.pi.5 <- mean.Bw.mat.4b.pi.5 + Bw.mat.5
  mean.Bw.diag.4b.pi.5 <- mean.Bw.diag.4b.pi.5 + Bw.diag.5
  MISE.mat.4b.pi.5 <- MISE.mat.4b.pi.5 + F.computeISE(im.4b.128, dens.pi.128.5)
  MISE.diag.4b.pi.5 <- MISE.diag.4b.pi.5 + F.computeISE(im.4b.128, dens.pi.diag.128.5)
  mean.Bw.mat.4b.pi.10 <- mean.Bw.mat.4b.pi.10 + Bw.mat.10
  mean.Bw.diag.4b.pi.10 <- mean.Bw.diag.4b.pi.10 + Bw.diag.10
  MISE.mat.4b.pi.10 <- MISE.mat.4b.pi.10 + F.computeISE(im.4b.128, dens.pi.128.10)
  MISE.diag.4b.pi.10 <- MISE.diag.4b.pi.10 + F.computeISE(im.4b.128, dens.pi.diag.128.10)
  time.4b.pi <- time.4b.pi + Sys.time() - st
  # ise = F.computeISE(im.4b.128, dens.128)
  # MISE.mat.4b <- MISE.mat.4b + ise
  # ise.diag = F.computeISE(im.4b.128, dens.diag.128)
  # MISE.diag.4b <- MISE.diag.4b + ise.diag
  # cat(ise, ise.diag, "\n")
}

time.4b.pi <- time.4b.pi/N
MISE.mat.4b.pi <- MISE.mat.4b.pi/N
MISE.diag.4b.pi <- MISE.diag.4b.pi/N
mean.Bw.mat.4b.pi <- mean.Bw.mat.4b.pi/N
mean.Bw.diag.4b.pi <- mean.Bw.diag.4b.pi/N
MISE.mat.4b.pi.1 <- MISE.mat.4b.pi.1/N
MISE.diag.4b.pi.1 <- MISE.diag.4b.pi.1/N
mean.Bw.mat.4b.pi.1 <- mean.Bw.mat.4b.pi.1/N
mean.Bw.diag.4b.pi.1 <- mean.Bw.diag.4b.pi.1/N
MISE.mat.4b.pi.2 <- MISE.mat.4b.pi.2/N
MISE.diag.4b.pi.2 <- MISE.diag.4b.pi.2/N
mean.Bw.mat.4b.pi.2 <- mean.Bw.mat.4b.pi.2/N
mean.Bw.diag.4b.pi.2 <- mean.Bw.diag.4b.pi.2/N
MISE.mat.4b.pi.3 <- MISE.mat.4b.pi.3/N
MISE.diag.4b.pi.3 <- MISE.diag.4b.pi.3/N
mean.Bw.mat.4b.pi.3 <- mean.Bw.mat.4b.pi.3/N
mean.Bw.diag.4b.pi.3 <- mean.Bw.diag.4b.pi.3/N
MISE.mat.4b.pi.4 <- MISE.mat.4b.pi.4/N
MISE.diag.4b.pi.4 <- MISE.diag.4b.pi.4/N
mean.Bw.mat.4b.pi.4 <- mean.Bw.mat.4b.pi.4/N
mean.Bw.diag.4b.pi.4 <- mean.Bw.diag.4b.pi.4/N
MISE.mat.4b.pi.5 <- MISE.mat.4b.pi.5/N
MISE.diag.4b.pi.5 <- MISE.diag.4b.pi.5/N
mean.Bw.mat.4b.pi.5 <- mean.Bw.mat.4b.pi.5/N
mean.Bw.diag.4b.pi.5 <- mean.Bw.diag.4b.pi.5/N
MISE.mat.4b.pi.10 <- MISE.mat.4b.pi.10/N
MISE.diag.4b.pi.10 <- MISE.diag.4b.pi.10/N
mean.Bw.mat.4b.pi.10 <- mean.Bw.mat.4b.pi.10/N
mean.Bw.diag.4b.pi.10 <- mean.Bw.diag.4b.pi.10/N









# 5b, scv

N <- 100

set.seed(887330)

MISE.mat.5b <- 0
MISE.diag.5b <- 0
MISE.mat.5b.1 <- 0
MISE.diag.5b.1 <- 0
MISE.mat.5b.2 <- 0
MISE.diag.5b.2 <- 0
MISE.mat.5b.3 <- 0
MISE.diag.5b.3 <- 0
MISE.mat.5b.4 <- 0
MISE.diag.5b.4 <- 0
MISE.mat.5b.5 <- 0
MISE.diag.5b.5 <- 0
MISE.mat.5b.10 <- 0
MISE.diag.5b.10 <- 0
mean.Bw.mat.5b <- diag(0,2,2)
mean.Bw.diag.5b <- diag(0,2,2)
mean.Bw.mat.5b.1 <- diag(0,2,2)
mean.Bw.diag.5b.1 <- diag(0,2,2)
mean.Bw.mat.5b.2 <- diag(0,2,2)
mean.Bw.diag.5b.2 <- diag(0,2,2)
mean.Bw.mat.5b.3 <- diag(0,2,2)
mean.Bw.diag.5b.3 <- diag(0,2,2)
mean.Bw.mat.5b.4 <- diag(0,2,2)
mean.Bw.diag.5b.4 <- diag(0,2,2)
mean.Bw.mat.5b.5 <- diag(0,2,2)
mean.Bw.diag.5b.5 <- diag(0,2,2)
mean.Bw.mat.5b.10 <- diag(0,2,2)
mean.Bw.diag.5b.10 <- diag(0,2,2)
time.5b <- 0

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.5b.128 <- as.im(lambda5b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.5b <- rpoispp(im.5b.128)
  pp.5b.df <- data.frame(x=pp.5b$x, y=pp.5b$y)
  Bw.mat <- Hscv(pp.5b.df)
  Bw.diag <- Hscv.diag(pp.5b.df)
  Bw.mat.1 <- Hscv.dir(pp.5b.df, alpha=1)
  Bw.diag.1 <- Hscv.diag.dir(pp.5b.df, alpha=1)
  Bw.mat.2 <- Hscv.dir(pp.5b.df, alpha=2)
  Bw.diag.2 <- Hscv.diag.dir(pp.5b.df, alpha=2)
  Bw.mat.3 <- Hscv.dir(pp.5b.df, alpha=3)
  Bw.diag.3 <- Hscv.diag.dir(pp.5b.df, alpha=3)
  Bw.mat.4 <- Hscv.dir(pp.5b.df, alpha=4)
  Bw.diag.4 <- Hscv.diag.dir(pp.5b.df, alpha=4)
  Bw.mat.5 <- Hscv.dir(pp.5b.df, alpha=5)
  Bw.diag.5 <- Hscv.diag.dir(pp.5b.df, alpha=5)
  Bw.mat.10 <- Hscv.dir(pp.5b.df, alpha=10)
  Bw.diag.10 <- Hscv.diag.dir(pp.5b.df, alpha=10)
  
  dens.128 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.mat)
  dens.diag.128 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.diag)
  dens.128.1 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.mat.1)
  dens.diag.128.1 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.diag.1)
  dens.128.2 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.mat.2)
  dens.diag.128.2 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.diag.2)
  dens.128.3 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.mat.3)
  dens.diag.128.3 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.diag.3)
  dens.128.4 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.mat.4)
  dens.diag.128.4 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.diag.4)
  dens.128.5 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.mat.5)
  dens.diag.128.5 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.diag.5)
  dens.128.10 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.mat.10)
  dens.diag.128.10 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.diag.10)
  
  mean.Bw.mat.5b <- mean.Bw.mat.5b + Bw.mat
  mean.Bw.diag.5b <- mean.Bw.diag.5b + Bw.diag
  MISE.mat.5b <- MISE.mat.5b + F.computeISE(im.5b.128, dens.128)
  MISE.diag.5b <- MISE.diag.5b + F.computeISE(im.5b.128, dens.diag.128)
  mean.Bw.mat.5b.1 <- mean.Bw.mat.5b.1 + Bw.mat.1
  mean.Bw.diag.5b.1 <- mean.Bw.diag.5b.1 + Bw.diag.1
  MISE.mat.5b.1 <- MISE.mat.5b.1 + F.computeISE(im.5b.128, dens.128.1)
  MISE.diag.5b.1 <- MISE.diag.5b.1 + F.computeISE(im.5b.128, dens.diag.128.1)
  mean.Bw.mat.5b.2 <- mean.Bw.mat.5b.2 + Bw.mat.2
  mean.Bw.diag.5b.2 <- mean.Bw.diag.5b.2 + Bw.diag.2
  MISE.mat.5b.2 <- MISE.mat.5b.2 + F.computeISE(im.5b.128, dens.128.2)
  MISE.diag.5b.2 <- MISE.diag.5b.2 + F.computeISE(im.5b.128, dens.diag.128.2)
  mean.Bw.mat.5b.3 <- mean.Bw.mat.5b.3 + Bw.mat.3
  mean.Bw.diag.5b.3 <- mean.Bw.diag.5b.3 + Bw.diag.3
  MISE.mat.5b.3 <- MISE.mat.5b.3 + F.computeISE(im.5b.128, dens.128.3)
  MISE.diag.5b.3 <- MISE.diag.5b.3 + F.computeISE(im.5b.128, dens.diag.128.3)
  mean.Bw.mat.5b.4 <- mean.Bw.mat.5b.4 + Bw.mat.4
  mean.Bw.diag.5b.4 <- mean.Bw.diag.5b.4 + Bw.diag.4
  MISE.mat.5b.4 <- MISE.mat.5b.4 + F.computeISE(im.5b.128, dens.128.4)
  MISE.diag.5b.4 <- MISE.diag.5b.4 + F.computeISE(im.5b.128, dens.diag.128.4)
  mean.Bw.mat.5b.5 <- mean.Bw.mat.5b.5 + Bw.mat.5
  mean.Bw.diag.5b.5 <- mean.Bw.diag.5b.5 + Bw.diag.5
  MISE.mat.5b.5 <- MISE.mat.5b.5 + F.computeISE(im.5b.128, dens.128.5)
  MISE.diag.5b.5 <- MISE.diag.5b.5 + F.computeISE(im.5b.128, dens.diag.128.5)
  mean.Bw.mat.5b.10 <- mean.Bw.mat.5b.10 + Bw.mat.10
  mean.Bw.diag.5b.10 <- mean.Bw.diag.5b.10 + Bw.diag.10
  MISE.mat.5b.10 <- MISE.mat.5b.10 + F.computeISE(im.5b.128, dens.128.10)
  MISE.diag.5b.10 <- MISE.diag.5b.10 + F.computeISE(im.5b.128, dens.diag.128.10)
  time.5b <- time.5b + Sys.time() -st
  # ise = F.computeISE(im.5b.128, dens.128)
  # MISE.mat.5b <- MISE.mat.5b + ise
  # ise.diag = F.computeISE(im.5b.128, dens.diag.128)
  # MISE.diag.5b <- MISE.diag.5b + ise.diag
  # cat(ise, ise.diag, "\n")
}

time.5b <- time.5b/N
MISE.mat.5b <- MISE.mat.5b/N
MISE.diag.5b <- MISE.diag.5b/N
mean.Bw.mat.5b <- mean.Bw.mat.5b/N
mean.Bw.diag.5b <- mean.Bw.diag.5b/N
MISE.mat.5b.1 <- MISE.mat.5b.1/N
MISE.diag.5b.1 <- MISE.diag.5b.1/N
mean.Bw.mat.5b.1 <- mean.Bw.mat.5b.1/N
mean.Bw.diag.5b.1 <- mean.Bw.diag.5b.1/N
MISE.mat.5b.2 <- MISE.mat.5b.2/N
MISE.diag.5b.2 <- MISE.diag.5b.2/N
mean.Bw.mat.5b.2 <- mean.Bw.mat.5b.2/N
mean.Bw.diag.5b.2 <- mean.Bw.diag.5b.2/N
MISE.mat.5b.3 <- MISE.mat.5b.3/N
MISE.diag.5b.3 <- MISE.diag.5b.3/N
mean.Bw.mat.5b.3 <- mean.Bw.mat.5b.3/N
mean.Bw.diag.5b.3 <- mean.Bw.diag.5b.3/N
MISE.mat.5b.4 <- MISE.mat.5b.4/N
MISE.diag.5b.4 <- MISE.diag.5b.4/N
mean.Bw.mat.5b.4 <- mean.Bw.mat.5b.4/N
mean.Bw.diag.5b.4 <- mean.Bw.diag.5b.4/N
MISE.mat.5b.5 <- MISE.mat.5b.5/N
MISE.diag.5b.5 <- MISE.diag.5b.5/N
mean.Bw.mat.5b.5 <- mean.Bw.mat.5b.5/N
mean.Bw.diag.5b.5 <- mean.Bw.diag.5b.5/N
MISE.mat.5b.10 <- MISE.mat.5b.10/N
MISE.diag.5b.10 <- MISE.diag.5b.10/N
mean.Bw.mat.5b.10 <- mean.Bw.mat.5b.10/N
mean.Bw.diag.5b.10 <- mean.Bw.diag.5b.10/N



# 5b, pi

N <- 100

set.seed(887330)

MISE.mat.5b.pi <- 0
MISE.diag.5b.pi <- 0
MISE.mat.5b.pi.1 <- 0
MISE.diag.5b.pi.1 <- 0
MISE.mat.5b.pi.2 <- 0
MISE.diag.5b.pi.2 <- 0
MISE.mat.5b.pi.3 <- 0
MISE.diag.5b.pi.3 <- 0
MISE.mat.5b.pi.4 <- 0
MISE.diag.5b.pi.4 <- 0
MISE.mat.5b.pi.5 <- 0
MISE.diag.5b.pi.5 <- 0
MISE.mat.5b.pi.10 <- 0
MISE.diag.5b.pi.10 <- 0
mean.Bw.mat.5b.pi <- diag(0,2,2)
mean.Bw.diag.5b.pi <- diag(0,2,2)
mean.Bw.mat.5b.pi.1 <- diag(0,2,2)
mean.Bw.diag.5b.pi.1 <- diag(0,2,2)
mean.Bw.mat.5b.pi.2 <- diag(0,2,2)
mean.Bw.diag.5b.pi.2 <- diag(0,2,2)
mean.Bw.mat.5b.pi.3 <- diag(0,2,2)
mean.Bw.diag.5b.pi.3 <- diag(0,2,2)
mean.Bw.mat.5b.pi.4 <- diag(0,2,2)
mean.Bw.diag.5b.pi.4 <- diag(0,2,2)
mean.Bw.mat.5b.pi.5 <- diag(0,2,2)
mean.Bw.diag.5b.pi.5 <- diag(0,2,2)
mean.Bw.mat.5b.pi.10 <- diag(0,2,2)
mean.Bw.diag.5b.pi.10 <- diag(0,2,2)
time.5b.pi <- 0

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.5b.128 <- as.im(lambda5b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.5b <- rpoispp(im.5b.128)
  pp.5b.df <- data.frame(x=pp.5b$x, y=pp.5b$y)
  Bw.pi.mat <- Hpi(pp.5b.df)
  Bw.pi.diag <- Hpi.diag(pp.5b.df)
  Bw.pi.mat.1 <- Hpi.dir(pp.5b.df, alpha=1)
  Bw.pi.diag.1 <- Hpi.diag.dir(pp.5b.df, alpha=1)
  Bw.pi.mat.2 <- Hpi.dir(pp.5b.df, alpha=2)
  Bw.pi.diag.2 <- Hpi.diag.dir(pp.5b.df, alpha=2)
  Bw.pi.mat.3 <- Hpi.dir(pp.5b.df, alpha=3)
  Bw.pi.diag.3 <- Hpi.diag.dir(pp.5b.df, alpha=3)
  Bw.pi.mat.4 <- Hpi.dir(pp.5b.df, alpha=4)
  Bw.pi.diag.4 <- Hpi.diag.dir(pp.5b.df, alpha=4)
  Bw.pi.mat.5 <- Hpi.dir(pp.5b.df, alpha=5)
  Bw.pi.diag.5 <- Hpi.diag.dir(pp.5b.df, alpha=5)
  Bw.pi.mat.10 <- Hpi.dir(pp.5b.df, alpha=10)
  Bw.pi.diag.10 <- Hpi.diag.dir(pp.5b.df, alpha=10)
  
  dens.pi.128 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.pi.mat)
  dens.pi.diag.128 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.pi.diag)
  dens.pi.128.1 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.pi.mat.1)
  dens.pi.diag.128.1 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.pi.diag.1)
  dens.pi.128.2 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.pi.mat.2)
  dens.pi.diag.128.2 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.pi.diag.2)
  dens.pi.128.3 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.pi.mat.3)
  dens.pi.diag.128.3 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.pi.diag.3)
  dens.pi.128.4 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.pi.mat.4)
  dens.pi.diag.128.4 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.pi.diag.4)
  dens.pi.128.5 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.pi.mat.5)
  dens.pi.diag.128.5 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.pi.diag.5)
  dens.pi.128.10 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.pi.mat.10)
  dens.pi.diag.128.10 <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.pi.diag.10)
  
  mean.Bw.mat.5b.pi <- mean.Bw.mat.5b.pi + Bw.mat
  mean.Bw.diag.5b.pi <- mean.Bw.diag.5b.pi + Bw.diag
  MISE.mat.5b.pi <- MISE.mat.5b.pi + F.computeISE(im.5b.128, dens.pi.128)
  MISE.diag.5b.pi <- MISE.diag.5b.pi + F.computeISE(im.5b.128, dens.pi.diag.128)
  mean.Bw.mat.5b.pi.1 <- mean.Bw.mat.5b.pi.1 + Bw.mat.1
  mean.Bw.diag.5b.pi.1 <- mean.Bw.diag.5b.pi.1 + Bw.diag.1
  MISE.mat.5b.pi.1 <- MISE.mat.5b.pi.1 + F.computeISE(im.5b.128, dens.pi.128.1)
  MISE.diag.5b.pi.1 <- MISE.diag.5b.pi.1 + F.computeISE(im.5b.128, dens.pi.diag.128.1)
  mean.Bw.mat.5b.pi.2 <- mean.Bw.mat.5b.pi.2 + Bw.mat.2
  mean.Bw.diag.5b.pi.2 <- mean.Bw.diag.5b.pi.2 + Bw.diag.2
  MISE.mat.5b.pi.2 <- MISE.mat.5b.pi.2 + F.computeISE(im.5b.128, dens.pi.128.2)
  MISE.diag.5b.pi.2 <- MISE.diag.5b.pi.2 + F.computeISE(im.5b.128, dens.pi.diag.128.2)
  mean.Bw.mat.5b.pi.3 <- mean.Bw.mat.5b.pi.3 + Bw.mat.3
  mean.Bw.diag.5b.pi.3 <- mean.Bw.diag.5b.pi.3 + Bw.diag.3
  MISE.mat.5b.pi.3 <- MISE.mat.5b.pi.3 + F.computeISE(im.5b.128, dens.pi.128.3)
  MISE.diag.5b.pi.3 <- MISE.diag.5b.pi.3 + F.computeISE(im.5b.128, dens.pi.diag.128.3)
  mean.Bw.mat.5b.pi.4 <- mean.Bw.mat.5b.pi.4 + Bw.mat.4
  mean.Bw.diag.5b.pi.4 <- mean.Bw.diag.5b.pi.4 + Bw.diag.4
  MISE.mat.5b.pi.4 <- MISE.mat.5b.pi.4 + F.computeISE(im.5b.128, dens.pi.128.4)
  MISE.diag.5b.pi.4 <- MISE.diag.5b.pi.4 + F.computeISE(im.5b.128, dens.pi.diag.128.4)
  mean.Bw.mat.5b.pi.5 <- mean.Bw.mat.5b.pi.5 + Bw.mat.5
  mean.Bw.diag.5b.pi.5 <- mean.Bw.diag.5b.pi.5 + Bw.diag.5
  MISE.mat.5b.pi.5 <- MISE.mat.5b.pi.5 + F.computeISE(im.5b.128, dens.pi.128.5)
  MISE.diag.5b.pi.5 <- MISE.diag.5b.pi.5 + F.computeISE(im.5b.128, dens.pi.diag.128.5)
  mean.Bw.mat.5b.pi.10 <- mean.Bw.mat.5b.pi.10 + Bw.mat.10
  mean.Bw.diag.5b.pi.10 <- mean.Bw.diag.5b.pi.10 + Bw.diag.10
  MISE.mat.5b.pi.10 <- MISE.mat.5b.pi.10 + F.computeISE(im.5b.128, dens.pi.128.10)
  MISE.diag.5b.pi.10 <- MISE.diag.5b.pi.10 + F.computeISE(im.5b.128, dens.pi.diag.128.10)
  time.5b.pi <- time.5b.pi + Sys.time() - st
  # ise = F.computeISE(im.5b.128, dens.128)
  # MISE.mat.5b <- MISE.mat.5b + ise
  # ise.diag = F.computeISE(im.5b.128, dens.diag.128)
  # MISE.diag.5b <- MISE.diag.5b + ise.diag
  # cat(ise, ise.diag, "\n")
}

time.5b.pi <- time.5b.pi/N
MISE.mat.5b.pi <- MISE.mat.5b.pi/N
MISE.diag.5b.pi <- MISE.diag.5b.pi/N
mean.Bw.mat.5b.pi <- mean.Bw.mat.5b.pi/N
mean.Bw.diag.5b.pi <- mean.Bw.diag.5b.pi/N
MISE.mat.5b.pi.1 <- MISE.mat.5b.pi.1/N
MISE.diag.5b.pi.1 <- MISE.diag.5b.pi.1/N
mean.Bw.mat.5b.pi.1 <- mean.Bw.mat.5b.pi.1/N
mean.Bw.diag.5b.pi.1 <- mean.Bw.diag.5b.pi.1/N
MISE.mat.5b.pi.2 <- MISE.mat.5b.pi.2/N
MISE.diag.5b.pi.2 <- MISE.diag.5b.pi.2/N
mean.Bw.mat.5b.pi.2 <- mean.Bw.mat.5b.pi.2/N
mean.Bw.diag.5b.pi.2 <- mean.Bw.diag.5b.pi.2/N
MISE.mat.5b.pi.3 <- MISE.mat.5b.pi.3/N
MISE.diag.5b.pi.3 <- MISE.diag.5b.pi.3/N
mean.Bw.mat.5b.pi.3 <- mean.Bw.mat.5b.pi.3/N
mean.Bw.diag.5b.pi.3 <- mean.Bw.diag.5b.pi.3/N
MISE.mat.5b.pi.4 <- MISE.mat.5b.pi.4/N
MISE.diag.5b.pi.4 <- MISE.diag.5b.pi.4/N
mean.Bw.mat.5b.pi.4 <- mean.Bw.mat.5b.pi.4/N
mean.Bw.diag.5b.pi.4 <- mean.Bw.diag.5b.pi.4/N
MISE.mat.5b.pi.5 <- MISE.mat.5b.pi.5/N
MISE.diag.5b.pi.5 <- MISE.diag.5b.pi.5/N
mean.Bw.mat.5b.pi.5 <- mean.Bw.mat.5b.pi.5/N
mean.Bw.diag.5b.pi.5 <- mean.Bw.diag.5b.pi.5/N
MISE.mat.5b.pi.10 <- MISE.mat.5b.pi.10/N
MISE.diag.5b.pi.10 <- MISE.diag.5b.pi.10/N
mean.Bw.mat.5b.pi.10 <- mean.Bw.mat.5b.pi.10/N
mean.Bw.diag.5b.pi.10 <- mean.Bw.diag.5b.pi.10/N



# 6b, scv


N <- 100

set.seed(887330)

MISE.mat.6b <- 0
MISE.diag.6b <- 0
MISE.mat.6b.1 <- 0
MISE.diag.6b.1 <- 0
MISE.mat.6b.2 <- 0
MISE.diag.6b.2 <- 0
MISE.mat.6b.3 <- 0
MISE.diag.6b.3 <- 0
MISE.mat.6b.4 <- 0
MISE.diag.6b.4 <- 0
MISE.mat.6b.5 <- 0
MISE.diag.6b.5 <- 0
MISE.mat.6b.10 <- 0
MISE.diag.6b.10 <- 0
mean.Bw.mat.6b <- diag(0,2,2)
mean.Bw.diag.6b <- diag(0,2,2)
mean.Bw.mat.6b.1 <- diag(0,2,2)
mean.Bw.diag.6b.1 <- diag(0,2,2)
mean.Bw.mat.6b.2 <- diag(0,2,2)
mean.Bw.diag.6b.2 <- diag(0,2,2)
mean.Bw.mat.6b.3 <- diag(0,2,2)
mean.Bw.diag.6b.3 <- diag(0,2,2)
mean.Bw.mat.6b.4 <- diag(0,2,2)
mean.Bw.diag.6b.4 <- diag(0,2,2)
mean.Bw.mat.6b.5 <- diag(0,2,2)
mean.Bw.diag.6b.5 <- diag(0,2,2)
mean.Bw.mat.6b.10 <- diag(0,2,2)
mean.Bw.diag.6b.10 <- diag(0,2,2)
time.6b <- 0

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.6b.128 <- as.im(lambda6b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.6b <- rpoispp(im.6b.128)
  pp.6b.df <- data.frame(x=pp.6b$x, y=pp.6b$y)
  Bw.mat <- Hscv(pp.6b.df)
  Bw.diag <- Hscv.diag(pp.6b.df)
  Bw.mat.1 <- Hscv.dir(pp.6b.df, alpha=1)
  Bw.diag.1 <- Hscv.diag.dir(pp.6b.df, alpha=1)
  Bw.mat.2 <- Hscv.dir(pp.6b.df, alpha=2)
  Bw.diag.2 <- Hscv.diag.dir(pp.6b.df, alpha=2)
  Bw.mat.3 <- Hscv.dir(pp.6b.df, alpha=3)
  Bw.diag.3 <- Hscv.diag.dir(pp.6b.df, alpha=3)
  Bw.mat.4 <- Hscv.dir(pp.6b.df, alpha=4)
  Bw.diag.4 <- Hscv.diag.dir(pp.6b.df, alpha=4)
  Bw.mat.5 <- Hscv.dir(pp.6b.df, alpha=5)
  Bw.diag.5 <- Hscv.diag.dir(pp.6b.df, alpha=5)
  Bw.mat.10 <- Hscv.dir(pp.6b.df, alpha=10)
  Bw.diag.10 <- Hscv.diag.dir(pp.6b.df, alpha=10)
  
  dens.128 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.mat)
  dens.diag.128 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.diag)
  dens.128.1 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.mat.1)
  dens.diag.128.1 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.diag.1)
  dens.128.2 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.mat.2)
  dens.diag.128.2 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.diag.2)
  dens.128.3 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.mat.3)
  dens.diag.128.3 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.diag.3)
  dens.128.4 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.mat.4)
  dens.diag.128.4 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.diag.4)
  dens.128.5 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.mat.5)
  dens.diag.128.5 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.diag.5)
  dens.128.10 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.mat.10)
  dens.diag.128.10 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.diag.10)
  
  mean.Bw.mat.6b <- mean.Bw.mat.6b + Bw.mat
  mean.Bw.diag.6b <- mean.Bw.diag.6b + Bw.diag
  MISE.mat.6b <- MISE.mat.6b + F.computeISE(im.6b.128, dens.128)
  MISE.diag.6b <- MISE.diag.6b + F.computeISE(im.6b.128, dens.diag.128)
  mean.Bw.mat.6b.1 <- mean.Bw.mat.6b.1 + Bw.mat.1
  mean.Bw.diag.6b.1 <- mean.Bw.diag.6b.1 + Bw.diag.1
  MISE.mat.6b.1 <- MISE.mat.6b.1 + F.computeISE(im.6b.128, dens.128.1)
  MISE.diag.6b.1 <- MISE.diag.6b.1 + F.computeISE(im.6b.128, dens.diag.128.1)
  mean.Bw.mat.6b.2 <- mean.Bw.mat.6b.2 + Bw.mat.2
  mean.Bw.diag.6b.2 <- mean.Bw.diag.6b.2 + Bw.diag.2
  MISE.mat.6b.2 <- MISE.mat.6b.2 + F.computeISE(im.6b.128, dens.128.2)
  MISE.diag.6b.2 <- MISE.diag.6b.2 + F.computeISE(im.6b.128, dens.diag.128.2)
  mean.Bw.mat.6b.3 <- mean.Bw.mat.6b.3 + Bw.mat.3
  mean.Bw.diag.6b.3 <- mean.Bw.diag.6b.3 + Bw.diag.3
  MISE.mat.6b.3 <- MISE.mat.6b.3 + F.computeISE(im.6b.128, dens.128.3)
  MISE.diag.6b.3 <- MISE.diag.6b.3 + F.computeISE(im.6b.128, dens.diag.128.3)
  mean.Bw.mat.6b.4 <- mean.Bw.mat.6b.4 + Bw.mat.4
  mean.Bw.diag.6b.4 <- mean.Bw.diag.6b.4 + Bw.diag.4
  MISE.mat.6b.4 <- MISE.mat.6b.4 + F.computeISE(im.6b.128, dens.128.4)
  MISE.diag.6b.4 <- MISE.diag.6b.4 + F.computeISE(im.6b.128, dens.diag.128.4)
  mean.Bw.mat.6b.5 <- mean.Bw.mat.6b.5 + Bw.mat.5
  mean.Bw.diag.6b.5 <- mean.Bw.diag.6b.5 + Bw.diag.5
  MISE.mat.6b.5 <- MISE.mat.6b.5 + F.computeISE(im.6b.128, dens.128.5)
  MISE.diag.6b.5 <- MISE.diag.6b.5 + F.computeISE(im.6b.128, dens.diag.128.5)
  mean.Bw.mat.6b.10 <- mean.Bw.mat.6b.10 + Bw.mat.10
  mean.Bw.diag.6b.10 <- mean.Bw.diag.6b.10 + Bw.diag.10
  MISE.mat.6b.10 <- MISE.mat.6b.10 + F.computeISE(im.6b.128, dens.128.10)
  MISE.diag.6b.10 <- MISE.diag.6b.10 + F.computeISE(im.6b.128, dens.diag.128.10)
  time.6b <- time.6b + Sys.time() - st
  
  # ise = F.computeISE(im.6b.128, dens.128)
  # MISE.mat.6b <- MISE.mat.6b + ise
  # ise.diag = F.computeISE(im.6b.128, dens.diag.128)
  # MISE.diag.6b <- MISE.diag.6b + ise.diag
  # cat(ise, ise.diag, "\n")
}

time.6b <- time.6b/N
MISE.mat.6b <- MISE.mat.6b/N
MISE.diag.6b <- MISE.diag.6b/N
mean.Bw.mat.6b <- mean.Bw.mat.6b/N
mean.Bw.diag.6b <- mean.Bw.diag.6b/N
MISE.mat.6b.1 <- MISE.mat.6b.1/N
MISE.diag.6b.1 <- MISE.diag.6b.1/N
mean.Bw.mat.6b.1 <- mean.Bw.mat.6b.1/N
mean.Bw.diag.6b.1 <- mean.Bw.diag.6b.1/N
MISE.mat.6b.2 <- MISE.mat.6b.2/N
MISE.diag.6b.2 <- MISE.diag.6b.2/N
mean.Bw.mat.6b.2 <- mean.Bw.mat.6b.2/N
mean.Bw.diag.6b.2 <- mean.Bw.diag.6b.2/N
MISE.mat.6b.3 <- MISE.mat.6b.3/N
MISE.diag.6b.3 <- MISE.diag.6b.3/N
mean.Bw.mat.6b.3 <- mean.Bw.mat.6b.3/N
mean.Bw.diag.6b.3 <- mean.Bw.diag.6b.3/N
MISE.mat.6b.4 <- MISE.mat.6b.4/N
MISE.diag.6b.4 <- MISE.diag.6b.4/N
mean.Bw.mat.6b.4 <- mean.Bw.mat.6b.4/N
mean.Bw.diag.6b.4 <- mean.Bw.diag.6b.4/N
MISE.mat.6b.5 <- MISE.mat.6b.5/N
MISE.diag.6b.5 <- MISE.diag.6b.5/N
mean.Bw.mat.6b.5 <- mean.Bw.mat.6b.5/N
mean.Bw.diag.6b.5 <- mean.Bw.diag.6b.5/N
MISE.mat.6b.10 <- MISE.mat.6b.10/N
MISE.diag.6b.10 <- MISE.diag.6b.10/N
mean.Bw.mat.6b.10 <- mean.Bw.mat.6b.10/N
mean.Bw.diag.6b.10 <- mean.Bw.diag.6b.10/N


# 6b, pi

N <- 100

set.seed(887330)

MISE.mat.6b.pi <- 0
MISE.diag.6b.pi <- 0
MISE.mat.6b.pi.1 <- 0
MISE.diag.6b.pi.1 <- 0
MISE.mat.6b.pi.2 <- 0
MISE.diag.6b.pi.2 <- 0
MISE.mat.6b.pi.3 <- 0
MISE.diag.6b.pi.3 <- 0
MISE.mat.6b.pi.4 <- 0
MISE.diag.6b.pi.4 <- 0
MISE.mat.6b.pi.5 <- 0
MISE.diag.6b.pi.5 <- 0
MISE.mat.6b.pi.10 <- 0
MISE.diag.6b.pi.10 <- 0
mean.Bw.mat.6b.pi <- diag(0,2,2)
mean.Bw.diag.6b.pi <- diag(0,2,2)
mean.Bw.mat.6b.pi.1 <- diag(0,2,2)
mean.Bw.diag.6b.pi.1 <- diag(0,2,2)
mean.Bw.mat.6b.pi.2 <- diag(0,2,2)
mean.Bw.diag.6b.pi.2 <- diag(0,2,2)
mean.Bw.mat.6b.pi.3 <- diag(0,2,2)
mean.Bw.diag.6b.pi.3 <- diag(0,2,2)
mean.Bw.mat.6b.pi.4 <- diag(0,2,2)
mean.Bw.diag.6b.pi.4 <- diag(0,2,2)
mean.Bw.mat.6b.pi.5 <- diag(0,2,2)
mean.Bw.diag.6b.pi.5 <- diag(0,2,2)
mean.Bw.mat.6b.pi.10 <- diag(0,2,2)
mean.Bw.diag.6b.pi.10 <- diag(0,2,2)
time.6b.pi <- 0

for (n in 1:N){
  cat(n, "\n")
  st = Sys.time()
  im.6b.128 <- as.im(lambda6b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.6b <- rpoispp(im.6b.128)
  pp.6b.df <- data.frame(x=pp.6b$x, y=pp.6b$y)
  Bw.pi.mat <- Hpi(pp.6b.df)
  Bw.pi.diag <- Hpi.diag(pp.6b.df)
  Bw.pi.mat.1 <- Hpi.dir(pp.6b.df, alpha=1)
  Bw.pi.diag.1 <- Hpi.diag.dir(pp.6b.df, alpha=1)
  Bw.pi.mat.2 <- Hpi.dir(pp.6b.df, alpha=2)
  Bw.pi.diag.2 <- Hpi.diag.dir(pp.6b.df, alpha=2)
  Bw.pi.mat.3 <- Hpi.dir(pp.6b.df, alpha=3)
  Bw.pi.diag.3 <- Hpi.diag.dir(pp.6b.df, alpha=3)
  Bw.pi.mat.4 <- Hpi.dir(pp.6b.df, alpha=4)
  Bw.pi.diag.4 <- Hpi.diag.dir(pp.6b.df, alpha=4)
  Bw.pi.mat.5 <- Hpi.dir(pp.6b.df, alpha=5)
  Bw.pi.diag.5 <- Hpi.diag.dir(pp.6b.df, alpha=5)
  Bw.pi.mat.10 <- Hpi.dir(pp.6b.df, alpha=10)
  Bw.pi.diag.10 <- Hpi.diag.dir(pp.6b.df, alpha=10)
  
  dens.pi.128 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.pi.mat)
  dens.pi.diag.128 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.pi.diag)
  dens.pi.128.1 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.pi.mat.1)
  dens.pi.diag.128.1 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.pi.diag.1)
  dens.pi.128.2 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.pi.mat.2)
  dens.pi.diag.128.2 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.pi.diag.2)
  dens.pi.128.3 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.pi.mat.3)
  dens.pi.diag.128.3 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.pi.diag.3)
  dens.pi.128.4 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.pi.mat.4)
  dens.pi.diag.128.4 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.pi.diag.4)
  dens.pi.128.5 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.pi.mat.5)
  dens.pi.diag.128.5 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.pi.diag.5)
  dens.pi.128.10 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.pi.mat.10)
  dens.pi.diag.128.10 <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.pi.diag.10)
  
  mean.Bw.mat.6b.pi <- mean.Bw.mat.6b.pi + Bw.pi.mat
  mean.Bw.diag.6b.pi <- mean.Bw.diag.6b.pi + Bw.pi.diag
  MISE.mat.6b.pi <- MISE.mat.6b.pi + F.computeISE(im.6b.128, dens.pi.128)
  MISE.diag.6b.pi <- MISE.diag.6b.pi + F.computeISE(im.6b.128, dens.pi.diag.128)
  mean.Bw.mat.6b.pi.1 <- mean.Bw.mat.6b.pi.1 + Bw.mat.1
  mean.Bw.diag.6b.pi.1 <- mean.Bw.diag.6b.pi.1 + Bw.diag.1
  MISE.mat.6b.pi.1 <- MISE.mat.6b.pi.1 + F.computeISE(im.6b.128, dens.pi.128.1)
  MISE.diag.6b.pi.1 <- MISE.diag.6b.pi.1 + F.computeISE(im.6b.128, dens.pi.diag.128.1)
  mean.Bw.mat.6b.pi.2 <- mean.Bw.mat.6b.pi.2 + Bw.mat.2
  mean.Bw.diag.6b.pi.2 <- mean.Bw.diag.6b.pi.2 + Bw.diag.2
  MISE.mat.6b.pi.2 <- MISE.mat.6b.pi.2 + F.computeISE(im.6b.128, dens.pi.128.2)
  MISE.diag.6b.pi.2 <- MISE.diag.6b.pi.2 + F.computeISE(im.6b.128, dens.pi.diag.128.2)
  mean.Bw.mat.6b.pi.3 <- mean.Bw.mat.6b.pi.3 + Bw.mat.3
  mean.Bw.diag.6b.pi.3 <- mean.Bw.diag.6b.pi.3 + Bw.diag.3
  MISE.mat.6b.pi.3 <- MISE.mat.6b.pi.3 + F.computeISE(im.6b.128, dens.pi.128.3)
  MISE.diag.6b.pi.3 <- MISE.diag.6b.pi.3 + F.computeISE(im.6b.128, dens.pi.diag.128.3)
  mean.Bw.mat.6b.pi.4 <- mean.Bw.mat.6b.pi.4 + Bw.mat.4
  mean.Bw.diag.6b.pi.4 <- mean.Bw.diag.6b.pi.4 + Bw.diag.4
  MISE.mat.6b.pi.4 <- MISE.mat.6b.pi.4 + F.computeISE(im.6b.128, dens.pi.128.4)
  MISE.diag.6b.pi.4 <- MISE.diag.6b.pi.4 + F.computeISE(im.6b.128, dens.pi.diag.128.4)
  mean.Bw.mat.6b.pi.5 <- mean.Bw.mat.6b.pi.5 + Bw.mat.5
  mean.Bw.diag.6b.pi.5 <- mean.Bw.diag.6b.pi.5 + Bw.diag.5
  MISE.mat.6b.pi.5 <- MISE.mat.6b.pi.5 + F.computeISE(im.6b.128, dens.pi.128.5)
  MISE.diag.6b.pi.5 <- MISE.diag.6b.pi.5 + F.computeISE(im.6b.128, dens.pi.diag.128.5)
  mean.Bw.mat.6b.pi.10 <- mean.Bw.mat.6b.pi.10 + Bw.mat.10
  mean.Bw.diag.6b.pi.10 <- mean.Bw.diag.6b.pi.10 + Bw.diag.10
  MISE.mat.6b.pi.10 <- MISE.mat.6b.pi.10 + F.computeISE(im.6b.128, dens.pi.128.10)
  MISE.diag.6b.pi.10 <- MISE.diag.6b.pi.10 + F.computeISE(im.6b.128, dens.pi.diag.128.10)
  time.6b.pi <- time.6b.pi + Sys.time() - st
  # ise = F.computeISE(im.6b.128, dens.128)
  # MISE.mat.6b <- MISE.mat.6b + ise
  # ise.diag = F.computeISE(im.6b.128, dens.diag.128)
  # MISE.diag.6b <- MISE.diag.6b + ise.diag
  # cat(ise, ise.diag, "\n")
}

time.6b.pi <- time.6b.pi/N
MISE.mat.6b.pi <- MISE.mat.6b.pi/N
MISE.diag.6b.pi <- MISE.diag.6b.pi/N
mean.Bw.mat.6b.pi <- mean.Bw.mat.6b.pi/N
mean.Bw.diag.6b.pi <- mean.Bw.diag.6b.pi/N
MISE.mat.6b.pi.1 <- MISE.mat.6b.pi.1/N
MISE.diag.6b.pi.1 <- MISE.diag.6b.pi.1/N
mean.Bw.mat.6b.pi.1 <- mean.Bw.mat.6b.pi.1/N
mean.Bw.diag.6b.pi.1 <- mean.Bw.diag.6b.pi.1/N
MISE.mat.6b.pi.2 <- MISE.mat.6b.pi.2/N
MISE.diag.6b.pi.2 <- MISE.diag.6b.pi.2/N
mean.Bw.mat.6b.pi.2 <- mean.Bw.mat.6b.pi.2/N
mean.Bw.diag.6b.pi.2 <- mean.Bw.diag.6b.pi.2/N
MISE.mat.6b.pi.3 <- MISE.mat.6b.pi.3/N
MISE.diag.6b.pi.3 <- MISE.diag.6b.pi.3/N
mean.Bw.mat.6b.pi.3 <- mean.Bw.mat.6b.pi.3/N
mean.Bw.diag.6b.pi.3 <- mean.Bw.diag.6b.pi.3/N
MISE.mat.6b.pi.4 <- MISE.mat.6b.pi.4/N
MISE.diag.6b.pi.4 <- MISE.diag.6b.pi.4/N
mean.Bw.mat.6b.pi.4 <- mean.Bw.mat.6b.pi.4/N
mean.Bw.diag.6b.pi.4 <- mean.Bw.diag.6b.pi.4/N
MISE.mat.6b.pi.5 <- MISE.mat.6b.pi.5/N
MISE.diag.6b.pi.5 <- MISE.diag.6b.pi.5/N
mean.Bw.mat.6b.pi.5 <- mean.Bw.mat.6b.pi.5/N
mean.Bw.diag.6b.pi.5 <- mean.Bw.diag.6b.pi.5/N
MISE.mat.6b.pi.10 <- MISE.mat.6b.pi.10/N
MISE.diag.6b.pi.10 <- MISE.diag.6b.pi.10/N
mean.Bw.mat.6b.pi.10 <- mean.Bw.mat.6b.pi.10/N
mean.Bw.diag.6b.pi.10 <- mean.Bw.diag.6b.pi.10/N









# 6a


N <- 100

set.seed(887330)

MISE.mat.6a.1 <- 0
MISE.diag.6a.1 <- 0
MISE.mat.6a.2 <- 0
MISE.diag.6a.2 <- 0
MISE.mat.6a.3 <- 0
MISE.diag.6a.3 <- 0
MISE.mat.6a.4 <- 0
MISE.diag.6a.4 <- 0
MISE.mat.6a.5 <- 0
MISE.diag.6a.5 <- 0
MISE.mat.6a.10 <- 0
MISE.diag.6a.10 <- 0
mean.Bw.mat.6a.1 <- diag(0,2,2)
mean.Bw.diag.6a.1 <- diag(0,2,2)
mean.Bw.mat.6a.2 <- diag(0,2,2)
mean.Bw.diag.6a.2 <- diag(0,2,2)
mean.Bw.mat.6a.3 <- diag(0,2,2)
mean.Bw.diag.6a.3 <- diag(0,2,2)
mean.Bw.mat.6a.4 <- diag(0,2,2)
mean.Bw.diag.6a.4 <- diag(0,2,2)
mean.Bw.mat.6a.5 <- diag(0,2,2)
mean.Bw.diag.6a.5 <- diag(0,2,2)
mean.Bw.mat.6a.10 <- diag(0,2,2)
mean.Bw.diag.6a.10 <- diag(0,2,2)

for (n in 1:N){
  cat(n, "\n")
  im.6a.128 <- as.im(lambda6a, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.6a <- rpoispp(im.6a.128)
  pp.6a.df <- data.frame(x=pp.6a$x, y=pp.6a$y)
  Bw.mat.1 <- Hscv.dir(pp.6a.df, alpha=1)
  Bw.diag.1 <- Hscv.diag.dir(pp.6a.df, alpha=1)
  Bw.mat.2 <- Hscv.dir(pp.6a.df, alpha=2)
  Bw.diag.2 <- Hscv.diag.dir(pp.6a.df, alpha=2)
  Bw.mat.3 <- Hscv.dir(pp.6a.df, alpha=3)
  Bw.diag.3 <- Hscv.diag.dir(pp.6a.df, alpha=3)
  Bw.mat.4 <- Hscv.dir(pp.6a.df, alpha=4)
  Bw.diag.4 <- Hscv.diag.dir(pp.6a.df, alpha=4)
  Bw.mat.5 <- Hscv.dir(pp.6a.df, alpha=5)
  Bw.diag.5 <- Hscv.diag.dir(pp.6a.df, alpha=5)
  Bw.mat.10 <- Hscv.dir(pp.6a.df, alpha=10)
  Bw.diag.10 <- Hscv.diag.dir(pp.6a.df, alpha=10)
  
  dens.128.1 <- density.ppp(pp.6a, dimyx=c(128,128), varcov=Bw.mat.1)
  dens.diag.128.1 <- density.ppp(pp.6a, dimyx=c(128,128), varcov=Bw.diag.1)
  dens.128.2 <- density.ppp(pp.6a, dimyx=c(128,128), varcov=Bw.mat.2)
  dens.diag.128.2 <- density.ppp(pp.6a, dimyx=c(128,128), varcov=Bw.diag.2)
  dens.128.3 <- density.ppp(pp.6a, dimyx=c(128,128), varcov=Bw.mat.3)
  dens.diag.128.3 <- density.ppp(pp.6a, dimyx=c(128,128), varcov=Bw.diag.3)
  dens.128.4 <- density.ppp(pp.6a, dimyx=c(128,128), varcov=Bw.mat.4)
  dens.diag.128.4 <- density.ppp(pp.6a, dimyx=c(128,128), varcov=Bw.diag.4)
  dens.128.5 <- density.ppp(pp.6a, dimyx=c(128,128), varcov=Bw.mat.5)
  dens.diag.128.5 <- density.ppp(pp.6a, dimyx=c(128,128), varcov=Bw.diag.5)
  dens.128.10 <- density.ppp(pp.6a, dimyx=c(128,128), varcov=Bw.mat.10)
  dens.diag.128.10 <- density.ppp(pp.6a, dimyx=c(128,128), varcov=Bw.diag.10)
  
  mean.Bw.mat.6a.1 <- mean.Bw.mat.6a.1 + Bw.mat.1
  mean.Bw.diag.6a.1 <- mean.Bw.diag.6a.1 + Bw.diag.1
  MISE.mat.6a.1 <- MISE.mat.6a.1 + F.computeISE(im.6a.128, dens.128.1)
  MISE.diag.6a.1 <- MISE.diag.6a.1 + F.computeISE(im.6a.128, dens.diag.128.1)
  mean.Bw.mat.6a.2 <- mean.Bw.mat.6a.2 + Bw.mat.2
  mean.Bw.diag.6a.2 <- mean.Bw.diag.6a.2 + Bw.diag.2
  MISE.mat.6a.2 <- MISE.mat.6a.2 + F.computeISE(im.6a.128, dens.128.2)
  MISE.diag.6a.2 <- MISE.diag.6a.2 + F.computeISE(im.6a.128, dens.diag.128.2)
  mean.Bw.mat.6a.3 <- mean.Bw.mat.6a.3 + Bw.mat.3
  mean.Bw.diag.6a.3 <- mean.Bw.diag.6a.3 + Bw.diag.3
  MISE.mat.6a.3 <- MISE.mat.6a.3 + F.computeISE(im.6a.128, dens.128.3)
  MISE.diag.6a.3 <- MISE.diag.6a.3 + F.computeISE(im.6a.128, dens.diag.128.3)
  mean.Bw.mat.6a.4 <- mean.Bw.mat.6a.4 + Bw.mat.4
  mean.Bw.diag.6a.4 <- mean.Bw.diag.6a.4 + Bw.diag.4
  MISE.mat.6a.4 <- MISE.mat.6a.4 + F.computeISE(im.6a.128, dens.128.4)
  MISE.diag.6a.4 <- MISE.diag.6a.4 + F.computeISE(im.6a.128, dens.diag.128.4)
  mean.Bw.mat.6a.5 <- mean.Bw.mat.6a.5 + Bw.mat.5
  mean.Bw.diag.6a.5 <- mean.Bw.diag.6a.5 + Bw.diag.5
  MISE.mat.6a.5 <- MISE.mat.6a.5 + F.computeISE(im.6a.128, dens.128.5)
  MISE.diag.6a.5 <- MISE.diag.6a.5 + F.computeISE(im.6a.128, dens.diag.128.5)
  mean.Bw.mat.6a.10 <- mean.Bw.mat.6a.10 + Bw.mat.10
  mean.Bw.diag.6a.10 <- mean.Bw.diag.6a.10 + Bw.diag.10
  MISE.mat.6a.10 <- MISE.mat.6a.10 + F.computeISE(im.6a.128, dens.128.10)
  MISE.diag.6a.10 <- MISE.diag.6a.10 + F.computeISE(im.6a.128, dens.diag.128.10)
  
  # ise = F.computeISE(im.6a.128, dens.128)
  # MISE.mat.6a <- MISE.mat.6a + ise
  # ise.diag = F.computeISE(im.6a.128, dens.diag.128)
  # MISE.diag.6a <- MISE.diag.6a + ise.diag
  # cat(ise, ise.diag, "\n")
}

MISE.mat.6a.1 <- MISE.mat.6a.1/N
MISE.diag.6a.1 <- MISE.diag.6a.1/N
mean.Bw.mat.6a.1 <- mean.Bw.mat.6a.1/N
mean.Bw.diag.6a.1 <- mean.Bw.diag.6a.1/N
MISE.mat.6a.2 <- MISE.mat.6a.2/N
MISE.diag.6a.2 <- MISE.diag.6a.2/N
mean.Bw.mat.6a.2 <- mean.Bw.mat.6a.2/N
mean.Bw.diag.6a.2 <- mean.Bw.diag.6a.2/N
MISE.mat.6a.3 <- MISE.mat.6a.3/N
MISE.diag.6a.3 <- MISE.diag.6a.3/N
mean.Bw.mat.6a.3 <- mean.Bw.mat.6a.3/N
mean.Bw.diag.6a.3 <- mean.Bw.diag.6a.3/N
MISE.mat.6a.4 <- MISE.mat.6a.4/N
MISE.diag.6a.4 <- MISE.diag.6a.4/N
mean.Bw.mat.6a.4 <- mean.Bw.mat.6a.4/N
mean.Bw.diag.6a.4 <- mean.Bw.diag.6a.4/N
MISE.mat.6a.5 <- MISE.mat.6a.5/N
MISE.diag.6a.5 <- MISE.diag.6a.5/N
mean.Bw.mat.6a.5 <- mean.Bw.mat.6a.5/N
mean.Bw.diag.6a.5 <- mean.Bw.diag.6a.5/N
MISE.mat.6a.10 <- MISE.mat.6a.10/N
MISE.diag.6a.10 <- MISE.diag.6a.10/N
mean.Bw.mat.6a.10 <- mean.Bw.mat.6a.10/N
mean.Bw.diag.6a.10 <- mean.Bw.diag.6a.10/N
