

#### July 26, 2021
#### Defining the intensity functions to reproduce the point patterns in
#### Fuentes-Santos' SJS 2015 paper

library(spatstat)
library(RandomFields)

lambda1a <- function(x, y){
    ### only the x coordinate
    intensity <- 1700 * exp(-3 * x)
    return(intensity)

}

lambda1b <- function(x, y){
    ### only the x coordinate
    intensity <- 3300 * exp(-3 * x)
    return(intensity)

}

lambda2a <- function(x, y){
    intensity <- 500 * dnorm(x, mean=(0.3-.2*y), sd=0.1) + 25
    return(intensity)
}


lambda2b <- function(x, y){
    intensity <- 1000 * dnorm(x, mean=(0.3-.2*y), sd=0.1) + 25
    return(intensity)
}


lambda3a <- function(x, y){
    intensity <- 500 * dnorm(x, mean=(0.3-.2*y), sd=0.02) + 25
    return(intensity)
}


lambda3b <- function(x, y){
    intensity <- 1000 * dnorm(x, mean=(0.3-.2*y), sd=0.02) + 25
    return(intensity)
}

set.seed(19388)
mymodel <- RMexp(var=0.01, scale=0.1)
x.seq <- seq(0,1, length=128)
y.seq <- seq(0,1, length=128)
my.sim.RF <- RFsimulate(mymodel, x=x.seq, y=y.seq)
my.sim.RF  <- matrix(my.sim.RF@data[,1], nrow=128, byrow=T)


lambda4a <- function(myRF){
    intensity <- exp(6+4 * myRF)
    return(intensity)
}


lambda4b <- function(myRF){
    intensity <- 2*exp(6+4 * myRF)
    return(intensity)
}


lambda5a <- function(x, y){
    intensity <- 275 * (dnorm(x, mean=y, sd=0.1) + dnorm(x, mean=1-y, sd=0.1) -
                        0.2 * dnorm(x, mean=0.5, sd=0.1)*dnorm(y, mean=0.5, sd=0.1)) + 30
    
    return(intensity)
}

lambda5b <- function(x, y){
    intensity <- 600 * (dnorm(x, mean=y, sd=0.1) + dnorm(x, mean=1-y, sd=0.1) -
                        0.2 * dnorm(x, mean=0.5, sd=0.1)*dnorm(y, mean=0.5, sd=0.1)) + 30

    return(intensity)
}


### note: in the paper, their image of the intensity for lambda6 is flipped
lambda6a<- function(x, y){
    intensity <- 275 * (dnorm(x, mean=y, sd=0.1) + dnorm(x, mean=1-y, sd=0.2) -
                        0.05 * dnorm(x, mean=0.5, sd=0.1)*dnorm(y, mean=0.5, sd=0.1)) + 50

    return(intensity)
}


lambda6b <- function(x, y){
    intensity <- 600 * (dnorm(x, mean=y, sd=0.1) + dnorm(x, mean=1-y, sd=0.2) -
                        0.05 * dnorm(x, mean=0.5, sd=0.1)*dnorm(y, mean=0.5, sd=0.1)) + 50

    return(intensity)
}

save.image("IntensityFunctions.RData")
