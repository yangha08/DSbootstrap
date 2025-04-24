require(ks)
require(mclust)
require(spatstat)
require(RandomFields) # 4b


# 1b, others

N <- 100

set.seed(887330)

ISE.1b.diggle <- rep(NA, N)
ISE.1b.mat.lscv <- rep(NA, N)
ISE.1b.diag.lscv <- rep(NA, N)
ISE.1b.scott <- rep(NA, N)
ISE.1b.lcv <- rep(NA, N)

for (n in 1:N){
  cat(n,"\t")
  im.1b.128 <- as.im(lambda1b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.1b <- rpoispp(im.1b.128)
  pp.1b.df <- data.frame(x=pp.1b$x, y=pp.1b$y)
  
  Bw.diggle <- bw.diggle(pp.1b)
  Bw.mat.lscv <- Hlscv(pp.1b.df)
  Bw.diag.lscv <- Hlscv.diag(pp.1b.df)
  Bw.sco <- bw.scott(pp.1b)
  Bw.lcv <- bw.ppl(pp.1b)
  
  dens.diggle <- density.ppp(pp.1b, dimyx=c(128,128), varcov=diag(Bw.diggle^2,2))
  dens.mat.lscv <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat.lscv)
  dens.diag.lscv <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag.lscv)
  dens.scott <- density.ppp(pp.1b, dimyx=c(128,128), varcov=diag(Bw.sco^2))
  dens.lcv <- density.ppp(pp.1b, dimyx=c(128,128), varcov=diag(Bw.lcv^2,2))
  
  ISE.1b.diggle[n] <- F.computeISE(im.1b.128, dens.diggle)
  ISE.1b.mat.lscv[n] <- F.computeISE(im.1b.128, dens.mat.lscv)
  ISE.1b.diag.lscv[n] <- F.computeISE(im.1b.128, dens.diag.lscv)
  ISE.1b.scott[n] <- F.computeISE(im.1b.128, dens.scott)
  ISE.1b.lcv[n] <- F.computeISE(im.1b.128, dens.lcv)
}

save(ISE.1b.diggle,
     ISE.1b.mat.lscv,
     ISE.1b.diag.lscv,
     ISE.1b.scott,
     ISE.1b.lcv,
     file = "ISE1other.Rdata")



# 1b thomas, others

N <- 100

set.seed(887330)

ISE.1b.thomas.diggle <- rep(NA, N)
ISE.1b.thomas.mat.lscv <- rep(NA, N)
ISE.1b.thomas.diag.lscv <- rep(NA, N)
ISE.1b.thomas.scott <- rep(NA, N)
ISE.1b.thomas.lcv <- rep(NA, N)

for (n in 1:N){
  cat(n,"\t")
  im.1b.128 <- as.im(lambda1b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.1b <- rThomas(kappa=200, scale=0.1, mu=5*im.1b.128*128*128/sum(im.1b.128))
  pp.1b.df <- data.frame(x=pp.1b$x, y=pp.1b$y)
  
  Bw.diggle <- bw.diggle(pp.1b)
  Bw.mat.lscv <- Hlscv(pp.1b.df)
  Bw.diag.lscv <- Hlscv.diag(pp.1b.df)
  Bw.sco <- bw.scott(pp.1b)
  Bw.lcv <- bw.ppl(pp.1b)
  
  dens.diggle <- density.ppp(pp.1b, dimyx=c(128,128), varcov=diag(Bw.diggle^2,2))
  dens.mat.lscv <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.mat.lscv)
  dens.diag.lscv <- density.ppp(pp.1b, dimyx=c(128,128), varcov=Bw.diag.lscv)
  dens.scott <- density.ppp(pp.1b, dimyx=c(128,128), varcov=diag(Bw.sco^2))
  dens.lcv <- density.ppp(pp.1b, dimyx=c(128,128), varcov=diag(Bw.lcv^2,2))
  
  ISE.1b.thomas.diggle[n] <- F.computeISE(im.1b.128, dens.diggle)
  ISE.1b.thomas.mat.lscv[n] <- F.computeISE(im.1b.128, dens.mat.lscv)
  ISE.1b.thomas.diag.lscv[n] <- F.computeISE(im.1b.128, dens.diag.lscv)
  ISE.1b.thomas.scott[n] <- F.computeISE(im.1b.128, dens.scott)
  ISE.1b.thomas.lcv[n] <- F.computeISE(im.1b.128, dens.lcv)
}

save(ISE.1b.thomas.diggle,
     ISE.1b.thomas.mat.lscv,
     ISE.1b.thomas.diag.lscv,
     ISE.1b.thomas.scott,
     ISE.1b.thomas.lcv,
     file = "ISE1thother.Rdata")


# 2b, others

N <- 100

set.seed(887330)

ISE.2b.diggle <- rep(NA, N)
ISE.2b.mat.lscv <- rep(NA, N)
ISE.2b.diag.lscv <- rep(NA, N)
ISE.2b.scott <- rep(NA, N)
ISE.2b.lcv <- rep(NA, N)

for (n in 1:N){
  cat(n,"\t")
  im.2b.128 <- as.im(lambda2b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.2b <- rpoispp(im.2b.128)
  pp.2b.df <- data.frame(x=pp.2b$x, y=pp.2b$y)
  
  Bw.diggle <- bw.diggle(pp.2b)
  Bw.mat.lscv <- Hlscv(pp.2b.df)
  Bw.diag.lscv <- Hlscv.diag(pp.2b.df)
  Bw.sco <- bw.scott(pp.2b)
  Bw.lcv <- bw.ppl(pp.2b)
  
  dens.diggle <- density.ppp(pp.2b, dimyx=c(128,128), varcov=diag(Bw.diggle^2,2))
  dens.mat.lscv <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat.lscv)
  dens.diag.lscv <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag.lscv)
  dens.scott <- density.ppp(pp.2b, dimyx=c(128,128), varcov=diag(Bw.sco^2))
  dens.lcv <- density.ppp(pp.2b, dimyx=c(128,128), varcov=diag(Bw.lcv^2,2))
  
  ISE.2b.diggle[n] <- F.computeISE(im.2b.128, dens.diggle)
  ISE.2b.mat.lscv[n] <- F.computeISE(im.2b.128, dens.mat.lscv)
  ISE.2b.diag.lscv[n] <- F.computeISE(im.2b.128, dens.diag.lscv)
  ISE.2b.scott[n] <- F.computeISE(im.2b.128, dens.scott)
  ISE.2b.lcv[n] <- F.computeISE(im.2b.128, dens.lcv)
}

save(ISE.2b.diggle,
     ISE.2b.mat.lscv,
     ISE.2b.diag.lscv,
     ISE.2b.scott,
     ISE.2b.lcv,
     file = "ISE2other.Rdata")


# 2b thomas, others

N <- 100

set.seed(887330)

ISE.2b.thomas.diggle <- rep(NA, N)
ISE.2b.thomas.mat.lscv <- rep(NA, N)
ISE.2b.thomas.diag.lscv <- rep(NA, N)
ISE.2b.thomas.scott <- rep(NA, N)
ISE.2b.thomas.lcv <- rep(NA, N)

for (n in 1:N){
  cat(n,"\t")
  im.2b.128 <- as.im(lambda2b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.2b <- rThomas(kappa=200, scale=0.1, mu=5*im.2b.128*128*128/sum(im.2b.128))
  pp.2b.df <- data.frame(x=pp.2b$x, y=pp.2b$y)
  
  Bw.diggle <- bw.diggle(pp.2b)
  Bw.mat.lscv <- Hlscv(pp.2b.df)
  Bw.diag.lscv <- Hlscv.diag(pp.2b.df)
  Bw.sco <- bw.scott(pp.2b)
  Bw.lcv <- bw.ppl(pp.2b)
  
  dens.diggle <- density.ppp(pp.2b, dimyx=c(128,128), varcov=diag(Bw.diggle^2,2))
  dens.mat.lscv <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.mat.lscv)
  dens.diag.lscv <- density.ppp(pp.2b, dimyx=c(128,128), varcov=Bw.diag.lscv)
  dens.scott <- density.ppp(pp.2b, dimyx=c(128,128), varcov=diag(Bw.sco^2))
  dens.lcv <- density.ppp(pp.2b, dimyx=c(128,128), varcov=diag(Bw.lcv^2,2))
  
  ISE.2b.thomas.diggle[n] <- F.computeISE(im.2b.128, dens.diggle)
  ISE.2b.thomas.mat.lscv[n] <- F.computeISE(im.2b.128, dens.mat.lscv)
  ISE.2b.thomas.diag.lscv[n] <- F.computeISE(im.2b.128, dens.diag.lscv)
  ISE.2b.thomas.scott[n] <- F.computeISE(im.2b.128, dens.scott)
  ISE.2b.thomas.lcv[n] <- F.computeISE(im.2b.128, dens.lcv)
}

save(ISE.2b.thomas.diggle,
     ISE.2b.thomas.mat.lscv,
     ISE.2b.thomas.diag.lscv,
     ISE.2b.thomas.scott,
     ISE.2b.thomas.lcv,
     file = "ISE2thother.Rdata")


# 3b, others

N <- 100

set.seed(887330)

ISE.3b.diggle <- rep(NA, N)
ISE.3b.mat.lscv <- rep(NA, N)
ISE.3b.diag.lscv <- rep(NA, N)
ISE.3b.scott <- rep(NA, N)
ISE.3b.lcv <- rep(NA, N)

for (n in 1:N){
  cat(n,"\t")
  im.3b.128 <- as.im(lambda3b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.3b <- rpoispp(im.3b.128)
  pp.3b.df <- data.frame(x=pp.3b$x, y=pp.3b$y)
  
  Bw.diggle <- bw.diggle(pp.3b)
  Bw.mat.lscv <- Hlscv(pp.3b.df)
  Bw.diag.lscv <- Hlscv.diag(pp.3b.df)
  Bw.sco <- bw.scott(pp.3b)
  Bw.lcv <- bw.ppl(pp.3b)
  
  dens.diggle <- density.ppp(pp.3b, dimyx=c(128,128), varcov=diag(Bw.diggle^2,2))
  dens.mat.lscv <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.mat.lscv)
  dens.diag.lscv <- density.ppp(pp.3b, dimyx=c(128,128), varcov=Bw.diag.lscv)
  dens.scott <- density.ppp(pp.3b, dimyx=c(128,128), varcov=diag(Bw.sco^2))
  dens.lcv <- density.ppp(pp.3b, dimyx=c(128,128), varcov=diag(Bw.lcv^2,2))
  
  ISE.3b.diggle[n] <- F.computeISE(im.3b.128, dens.diggle)
  ISE.3b.mat.lscv[n] <- F.computeISE(im.3b.128, dens.mat.lscv)
  ISE.3b.diag.lscv[n] <- F.computeISE(im.3b.128, dens.diag.lscv)
  ISE.3b.scott[n] <- F.computeISE(im.3b.128, dens.scott)
  ISE.3b.lcv[n] <- F.computeISE(im.3b.128, dens.lcv)
}

save(ISE.3b.diggle,
     ISE.3b.mat.lscv,
     ISE.3b.diag.lscv,
     ISE.3b.scott,
     ISE.3b.lcv,
     file = "ISE3other.Rdata")



# 4b, others

N <- 100

set.seed(887330)

ISE.4b.diggle <- rep(NA, N)
ISE.4b.mat.lscv <- rep(NA, N)
ISE.4b.diag.lscv <- rep(NA, N)
ISE.4b.scott <- rep(NA, N)
ISE.4b.lcv <- rep(NA, N)

for (n in 1:N){
  cat(n,"\t")
  im.4b.128 <- as.im(lambda4b(my.sim.RF), W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.4b <- rpoispp(im.4b.128)
  pp.4b.df <- data.frame(x=pp.4b$x, y=pp.4b$y)
  
  Bw.diggle <- bw.diggle(pp.4b)
  Bw.mat.lscv <- Hlscv(pp.4b.df)
  Bw.diag.lscv <- Hlscv.diag(pp.4b.df)
  Bw.sco <- bw.scott(pp.4b)
  Bw.lcv <- bw.ppl(pp.4b)
  
  dens.diggle <- density.ppp(pp.4b, dimyx=c(128,128), varcov=diag(Bw.diggle^2,2))
  dens.mat.lscv <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.mat.lscv)
  dens.diag.lscv <- density.ppp(pp.4b, dimyx=c(128,128), varcov=Bw.diag.lscv)
  dens.scott <- density.ppp(pp.4b, dimyx=c(128,128), varcov=diag(Bw.sco^2))
  dens.lcv <- density.ppp(pp.4b, dimyx=c(128,128), varcov=diag(Bw.lcv^2,2))
  
  ISE.4b.diggle[n] <- F.computeISE(im.4b.128, dens.diggle)
  ISE.4b.mat.lscv[n] <- F.computeISE(im.4b.128, dens.mat.lscv)
  ISE.4b.diag.lscv[n] <- F.computeISE(im.4b.128, dens.diag.lscv)
  ISE.4b.scott[n] <- F.computeISE(im.4b.128, dens.scott)
  ISE.4b.lcv[n] <- F.computeISE(im.4b.128, dens.lcv)
}

save(ISE.4b.diggle,
     ISE.4b.mat.lscv,
     ISE.4b.diag.lscv,
     ISE.4b.scott,
     ISE.4b.lcv,
     file = "ISE4other.Rdata")



# 5b, others

N <- 100

set.seed(887330)

ISE.5b.diggle <- rep(NA, N)
ISE.5b.mat.lscv <- rep(NA, N)
ISE.5b.diag.lscv <- rep(NA, N)
ISE.5b.scott <- rep(NA, N)
ISE.5b.lcv <- rep(NA, N)

for (n in 1:N){
  cat(n,"\t")
  im.5b.128 <- as.im(lambda5b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.5b <- rpoispp(im.5b.128)
  pp.5b.df <- data.frame(x=pp.5b$x, y=pp.5b$y)
  
  Bw.diggle <- bw.diggle(pp.5b)
  Bw.mat.lscv <- Hlscv(pp.5b.df)
  Bw.diag.lscv <- Hlscv.diag(pp.5b.df)
  Bw.sco <- bw.scott(pp.5b)
  Bw.lcv <- bw.ppl(pp.5b)
  
  dens.diggle <- density.ppp(pp.5b, dimyx=c(128,128), varcov=diag(Bw.diggle^2,2))
  dens.mat.lscv <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.mat.lscv)
  dens.diag.lscv <- density.ppp(pp.5b, dimyx=c(128,128), varcov=Bw.diag.lscv)
  dens.scott <- density.ppp(pp.5b, dimyx=c(128,128), varcov=diag(Bw.sco^2))
  dens.lcv <- density.ppp(pp.5b, dimyx=c(128,128), varcov=diag(Bw.lcv^2,2))
  
  ISE.5b.diggle[n] <- F.computeISE(im.5b.128, dens.diggle)
  ISE.5b.mat.lscv[n] <- F.computeISE(im.5b.128, dens.mat.lscv)
  ISE.5b.diag.lscv[n] <- F.computeISE(im.5b.128, dens.diag.lscv)
  ISE.5b.scott[n] <- F.computeISE(im.5b.128, dens.scott)
  ISE.5b.lcv[n] <- F.computeISE(im.5b.128, dens.lcv)
}

save(ISE.5b.diggle,
     ISE.5b.mat.lscv,
     ISE.5b.diag.lscv,
     ISE.5b.scott,
     ISE.5b.lcv,
     file = "ISE5other.Rdata")



# 6b, others

N <- 100

set.seed(887330)

ISE.6b.diggle <- rep(NA, N)
ISE.6b.mat.lscv <- rep(NA, N)
ISE.6b.diag.lscv <- rep(NA, N)
ISE.6b.scott <- rep(NA, N)
ISE.6b.lcv <- rep(NA, N)

for (n in 1:N){
  cat(n,"\t")
  im.6b.128 <- as.im(lambda6b, W=owin(c(0,1), c(0,1)), dimyx=c(128,128))
  pp.6b <- rpoispp(im.6b.128)
  pp.6b.df <- data.frame(x=pp.6b$x, y=pp.6b$y)
  
  Bw.diggle <- bw.diggle(pp.6b)
  Bw.mat.lscv <- Hlscv(pp.6b.df)
  Bw.diag.lscv <- Hlscv.diag(pp.6b.df)
  Bw.sco <- bw.scott(pp.6b)
  Bw.lcv <- bw.ppl(pp.6b)
  
  dens.diggle <- density.ppp(pp.6b, dimyx=c(128,128), varcov=diag(Bw.diggle^2,2))
  dens.mat.lscv <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.mat.lscv)
  dens.diag.lscv <- density.ppp(pp.6b, dimyx=c(128,128), varcov=Bw.diag.lscv)
  dens.scott <- density.ppp(pp.6b, dimyx=c(128,128), varcov=diag(Bw.sco^2))
  dens.lcv <- density.ppp(pp.6b, dimyx=c(128,128), varcov=diag(Bw.lcv^2,2))
  
  ISE.6b.diggle[n] <- F.computeISE(im.6b.128, dens.diggle)
  ISE.6b.mat.lscv[n] <- F.computeISE(im.6b.128, dens.mat.lscv)
  ISE.6b.diag.lscv[n] <- F.computeISE(im.6b.128, dens.diag.lscv)
  ISE.6b.scott[n] <- F.computeISE(im.6b.128, dens.scott)
  ISE.6b.lcv[n] <- F.computeISE(im.6b.128, dens.lcv)
}

save(ISE.6b.diggle,
     ISE.6b.mat.lscv,
     ISE.6b.diag.lscv,
     ISE.6b.scott,
     ISE.6b.lcv,
     file = "ISE6other.Rdata")
