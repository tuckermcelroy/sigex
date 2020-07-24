########################
#### Script for NDC Data
########################

## wipe
rm(list=ls())

library(devtools)

setwd("C:\\Users\\Tucker\\Documents\\GitHub\\sigex")
load_all(".")

######################
### Part I: load data

# automatic

#############################################################
### Part II: Metadata Specifications and Exploratory Analysis

start.date <- c(1992,1)
end.date <- c(2020,5)
period <- 12

## create ts object and plot
dataALL.ts <- sigex.load(ndc,start.date,period,c("Shipments","NewOrders"),TRUE)


#############################
## select span and transforms

## all data for NE-MW with log transform
transform <- "log"
aggregate <- FALSE
subseries <- c(1,2)
begin.date <- start(dataALL.ts)
end.date <- end(dataALL.ts)
range <- NULL
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)

#######################
## spectral exploratory

## levels
par(mfrow=c(2,1))
for(i in subseries) { sigex.specar(data.ts,FALSE,i,period) }
dev.off()

## growth rates
par(mfrow=c(2,1))
for(i in subseries) {	sigex.specar(data.ts,TRUE,i,period) }
dev.off()


###############################
### Part III: Model Declaration

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]

#######################################
## Default Model: Trend-Cycle-Irregular

## model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"trend",c(1,-1))
mdl <- sigex.add(mdl,seq(1,N),"cycleBAL",1,0,"cycle",1)
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"irregular",1)
# regressors:
mdl <- sigex.meaninit(mdl,data.ts,0)


#############################
### Part IV: Model Estimation

constraint <- NULL
par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## manage output
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]
resid.mle <- sigex.resid(psi.mle,mdl,data.ts)[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
acf(resid.mle,lag.max=40)

## examine condition numbers
log(sigex.conditions(data.ts,psi.mle,mdl))

## model checking
sigex.portmanteau(resid.mle,4*period,length(psi.mle))
sigex.gausscheck(resid.mle)

## bundle
analysis.mle <- sigex.bundle(data.ts,transform,mdl.mle,psi.mle)



###########################################
### Part V: Signal Extraction based on MLE

## load up the MLE fit for signal extraction
data.ts <- analysis.mle[[1]]
mdl <- analysis.mle[[3]]
psi <- analysis.mle[[4]]
param <- sigex.psi2par(psi,mdl,data.ts)


####################################################################
############################# METHOD 2: FORECASTING and WK SIGEX

#grid <- 70000	# high accuracy, close to method 1
#grid <- 700	# low accuracy, but pretty fast
grid <- 7000	# need grid > filter length
window <- 100
horizon <- 0
target <- array(diag(N),c(N,N,1))

extract.trend <- sigex.wkextract(psi,mdl,data.ts,1,target,grid,window,horizon,TRUE)
extract.cycle <- sigex.wkextract(psi,mdl,data.ts,2,target,grid,window,horizon,TRUE)
extract.irr <- sigex.wkextract(psi,mdl,data.ts,3,target,grid,window,horizon,TRUE)
extract.noncyc <- sigex.wkextract(psi,mdl,data.ts,c(1,3),target,grid,window,horizon,TRUE)
extract.hp <- sigex.wkextract(psi,mdl,data.ts,c(2,3),target,grid,window,horizon,TRUE)
extract.lp <- sigex.wkextract(psi,mdl,data.ts,c(1,2),target,grid,window,horizon,TRUE)


#########################################
### get fixed effects

reg.trend <- NULL
for(i in 1:N) { reg.trend <- cbind(reg.trend,sigex.fixed(data,mdl,i,param,"Trend")) }

## plotting
trendcol <- "tomato"
cyccol <- "orchid"
seascol <- "seagreen"
sacol <- "navyblue"
fade <- 60
#pdf(file="NdcSignals.pdf")
par(mfrow=c(2,2))
for(i in 1:N)
{
  plot(data.ts[,i],xlab="Year",ylab="",ylim=c(min(data.ts[,i])-25,max(data.ts[,i])),lwd=1)
  sigex.graph(extract.trend,reg.trend,begin.date,period,i,0,sacol,fade)
  sigex.graph(extract.lp,reg.trend,begin.date,period,i,0,trendcol,fade)
  sigex.graph(extract.cycle,NULL,begin.date,period,i,min(data.ts[,i])-10,seascol,fade)
}
#dev.off()


## spectral diagnostics: trend
par(mfrow=c(2,2))
for(i in 1:N)
{
  sigex.specar(ts(extract.trend[[1]],frequency=period,names=colnames(data.ts)),FALSE,i,period)
}
#dev.off()

## spectral diagnostics: sa
par(mfrow=c(2,2))
for(i in 1:N)
{
  sigex.specar(ts(extract.sa[[1]],frequency=period,names=colnames(data.ts)),FALSE,i,period)
}
#dev.off()

## transfer function analysis
grid <- 200
#pdf(file="StartsTrendfrf.pdf")
frf.trend <- sigex.getfrf(data.ts,param,mdl,1,TRUE,grid)
#dev.off()
#pdf(file="StartsSeasfrf.pdf")
frf.seas <- sigex.getfrf(data.ts,param,mdl,seq(2,7),TRUE,grid)
#dev.off()
#pdf(file="StartsSAfrf.pdf")
frf.sa <- sigex.getfrf(data.ts,param,mdl,c(1,8),TRUE,grid)
#dev.off()

## filter analysis
len <- 50
target <- array(diag(N),c(N,N,1))
#pdf(file="StartsTrendwk.pdf")
wk.trend <- sigex.wk(data.ts,param,mdl,1,target,TRUE,grid,len)
#dev.off()
#pdf(file="StartsSeaswk.pdf")
wk.seas <- sigex.wk(data.ts,param,mdl,seq(2,7),target,TRUE,grid,len)
#dev.off()
#pdf(file="StartsSAwk.pdf")
wk.sa <- sigex.wk(data.ts,param,mdl,c(1,8),target,TRUE,grid,len)
#dev.off()


