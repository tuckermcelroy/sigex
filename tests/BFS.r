#########################################
###  Script for Weekly BFS Data
#########################################

## wipe
rm(list=ls())

library(devtools)

setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex")
load_all(".")

######################
### Part I: load data

load("C:\\Users\\neide\\OneDrive\\Documents\\Research\\bfs.RData")


#############################################################
### Part II: Metadata Specifications and Exploratory Analysis

begin <- c(2006,1) 
end <- c(2020,18)
period <- 52
#first.day <- 1

## create ts object and plot
dataALL.ts <- sigex.load(bfs,begin,period,"bfs",TRUE)

#############################
## select span and transforms 

transform <- "log"
aggregate <- FALSE
subseries <- 1
range <- NULL
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)


#######################
## spectral exploratory

## levels
par(mfrow=c(1,1))
for(i in 1:length(subseries))
{
  sigex.specar(data.ts,FALSE,i,period)
}
dev.off()


## growth rates
par(mfrow=c(1,1))
for(i in 1:length(subseries))
{
  sigex.specar(data.ts,TRUE,i,period)
}
dev.off()

###############################
### Part III: Model Declaration

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]

############################## 
## Generate holiday regressors

#easter.dates <- read.table("data\\easter500.txt")
#easter.reg1 <- gethol(easter.dates,0,0,start.date,end.date)
#easter.reg2 <- gethol(easter.dates,8,-1,start.date,end.date)


##############
## Basic Model

# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"sarma",c(1,1,1,1,52),list(1,1,1,1),"process",1)
mdl <- sigex.meaninit(mdl,data.ts,0)	

##################################
### PART IV: Model Fitting

constraint <- NULL
par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:  
fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## manage output
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint) 
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]

## input parameter from previous fit (MLE on entire span)
#  divergence:    -2084.366 lik
#psi.mle <- c(-3.86830013490459, 6.08187233213583, 3.6806305674002, 3.98416618713543, 
#             1.65616209571557, 10.9517270937398)
#par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

##  model checking 
resid.mle <- sigex.resid(psi.mle,mdl,data.ts)[[1]]
resid.mle <- sigex.load(t(Re(resid.mle)),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*53,plot=FALSE)$acf

#pdf(file="retResidAcf.pdf",height=10,width=10)
par(mfrow=c(N,N),mar=c(3,2,2,0)+0.1,cex.lab=.8,cex.axis=.5,bty="n")
for(j in 1:N)
{
  for(k in 1:N)
  {
    plot.ts(resid.acf[,j,k],ylab="",xlab="Lag",ylim=c(-1,1),cex=.5)
    abline(h=1.96/sqrt(T),lty=3)
    abline(h=-1.96/sqrt(T),lty=3)
  }
}
#dev.off()

# bundle 
analysis.mle <- sigex.bundle(data.ts,transform,mdl,psi.mle)


##########################################
### Part V: Signal Extraction  

setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\Casting\\Figures")

## load up the fitted model for signal extraction
data.ts <- analysis.mle[[1]]
mdl <- analysis.mle[[3]]
psi <- analysis.mle[[4]]
param <- sigex.psi2par(psi,mdl,data.ts)

## define SA weekly filter
sa.filter <- array(c(1,rep(2,7),1)/(2*7),c(1,1,9))
len <- 4
shift <- len

td.low <- sigex.adhocextract(psi,mdl,data.ts,sa.filter,shift,0,TRUE)
td.low.daily <- list()
td.low.daily[[1]] <- sigex.weekly2daily(ts(td.low[[1]],start=start(data.ts),frequency=frequency(data.ts)),first.day)
td.low.daily[[2]] <- sigex.weekly2daily(ts(td.low[[2]],start=start(data.ts),frequency=frequency(data.ts)),first.day)
td.low.daily[[3]] <- sigex.weekly2daily(ts(td.low[[3]],start=start(data.ts),frequency=frequency(data.ts)),first.day)

