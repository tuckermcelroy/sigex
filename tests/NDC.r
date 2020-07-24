########################
#### Script for NDC Data
########################

## wipe
rm(list=ls())

library(devtools)

#setwd("C:\\Users\\Tucker\\Documents\\GitHub\\sigex")
setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex")
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

## all data with log transform
transform <- "log"
aggregate <- FALSE
subseries <- c(1,2)
begin.date <- start(dataALL.ts)
end.date <- end(dataALL.ts)
range <- NULL
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)


###############################
### Part III: Model Declaration

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]

#######################################
## Default Model: Trend-Cycle-Irregular

## model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"arma",c(1,1),0,"trend",c(1,-1))
mdl <- sigex.add(mdl,seq(1,N),"bal",1,c(0,1,0,1),"cycle",1)
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"irregular",1)
# regressors:
mdl <- sigex.meaninit(mdl,data.ts,0)

# HERE  debug VARMA ???

## model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"varma",c(2,1),list(c(1,1),1),"process",c(1,-1))
# regressors:
mdl <- sigex.meaninit(mdl,data.ts,0)


#############################
### Part IV: Model Estimation

constraint <- NULL
par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## MLE fitting results
#  divergence:    -3734.032
#psi.mle <- c(1.49058304904929, -8.32792347054058, -7.90797764406646, 0.0499570161046461,
#    -10.0917624787578, -9.56466665775856, 0.324064725662732, -9.25526788166362,
#    -5.8707556709411, -6.16202241709588, -3.3411108685506, 0.00154937934174211,
#    0.00107578875176897)
#par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

## manage output
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]

## residual analysis
resid.mle <- sigex.resid(psi.mle,mdl,data.ts)[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*period,plot=TRUE)$acf

## examine condition numbers
log(sigex.conditions(data.ts,psi.mle,mdl))

## model checking
sigex.portmanteau(resid.mle,4*period,length(psi.mle))
sigex.gausscheck(resid.mle)

## check on standard errors and get t statistics
print(eigen(hess)$values)
tstats <- sigex.tstats(mdl,psi.mle,hess,constraint)
print(tstats)

## bundle
analysis.mle <- sigex.bundle(data.ts,transform,mdl,psi.mle)



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

