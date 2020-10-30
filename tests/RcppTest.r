## Test script for RCPP material

## wipe
rm(list=ls())

library(devtools)

setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex")
load_all(".")

###############################
#### Test code for rcpp stuff

library(Rcpp)

#setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex\\src")
#sourceCpp('test.cpp')
#getRe(rep(1i,3))

#setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex\\R")
#sourceCpp('test.cpp')
#getRe(rep(1i,3))

setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\SigExNew")
#setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex\\cpp")
sourceCpp('mvar_midcast.cpp')

# testing of mvar_midcast.cpp
x.acf <- array(0,c(2,2,3))
x.acf[,,1] <- 6*diag(2)
x.acf[,,2] <- -2*diag(2)
x.acf[,,3] <- 1*diag(2)
z <- array(0,c(2,3))
z[,1] <- c(2,2+1i)
z[,2] <- c(1,-1)
z[,3] <- c(3+1i,4+1i)
delta <- c(1,-1)
debug <- TRUE
mvar_midcast(x.acf,z,delta,debug)

x.acf <- array(0,c(2,2,3))
x.acf[,,1] <- 6*diag(2)
x.acf[,,2] <- -2*diag(2)
x.acf[,,3] <- 1*diag(2)
z <- array(0,c(2,3))
z[,1] <- c(2,5)
z[,2] <- c(1,-1)
z[,3] <- c(3,4)
delta <- c(1,-2,1)
debug <- TRUE
mvar_midcast(x.acf,z,delta,debug)




# VARMAauto tests

setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\SigExNew")
sourceCpp('VARMAauto.cpp')

N <- 7
p <- 100
q <- 0
phi <- array(rnorm(N*N*p),c(N,N,p))
theta <- NULL
sigma <- diag(N)
param <- cbind(matrix(phi,nrow=N),sigma)

tic()
temp <- VARMA_auto(param,p,0,20)
toc()

tic()
out <- VARMAauto(phi,theta,sigma,20)
toc()



