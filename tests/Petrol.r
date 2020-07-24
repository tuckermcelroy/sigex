###########################
#### Script for Petrol Data
###########################

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

start.date <- c(1973,1)
period <- 12

## create ts object and plot
dataALL.ts <- sigex.load(petrol[,c(1,2)],start.date,period,c("Consumption","Imports"),TRUE)

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


################################
## Default Model: Related Trends

## model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"trend",c(1,-1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"irregular",1)
# regressors:
mdl <- sigex.meaninit(mdl,data.ts,0)

#################################
## Alternate Model: Common Trends

## model construction
mdl2 <- NULL
mdl2 <- sigex.add(mdl2,1,"arma",c(0,0),0,"trend",c(1,-1))
mdl2 <- sigex.add(mdl2,seq(1,N),"arma",c(0,0),0,"irregular",1)
# regressors:
mdl2 <- sigex.meaninit(mdl2,data.ts,0)



###################################################
### Part IV: MLE Estimation of Related Trends Model

## parameter initialization
constraint <- NULL
par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## MLE fitting results
#  divergence:    -5043.938
#psi.mle <- c(1.1893661498523, -8.51181967077802, -5.78741462564482, 0.211097637688758,
#             -6.79777756698125, -6.77733857261275, -0.000182396168600933,
#             0.000440082914531049)
#par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

## manage output
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]

##  model checking

# residual analysis
resid.mle <- sigex.resid(psi.mle,mdl,data.ts)[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*period,plot=FALSE)$acf
sigex.portmanteau(resid.mle,4*period,length(psi.mle))
sigex.gausscheck(resid.mle)
acf(resid.mle,lag.max=4*period)

## check on standard errors and condition numbers
print(eigen(hess)$values)
tstats <- sigex.tstats(mdl,psi.mle,hess,constraint)
print(tstats)
taus <- log(sigex.conditions(data.ts,psi.mle,mdl))
print(taus)

## bundle
analysis.mle <- sigex.bundle(data.ts,transform,mdl,psi.mle)



#################################################
### Part V: MLE Estimation of Common Trends Model

## parameter initialization
constraint <- NULL
par.mle2 <- sigex.default(mdl2,data.ts,constraint)
psi.mle2 <- sigex.par2psi(par.mle2,mdl2)

## run fitting:
fit.mle2 <- sigex.mlefit(data.ts,par.mle2,constraint,mdl2,"bfgs",debug=TRUE)

## MLE fitting results
#  divergence:    -5043.938
#psi.mle <- c(1.1893661498523, -8.51181967077802, -5.78741462564482, 0.211097637688758,
#             -6.79777756698125, -6.77733857261275, -0.000182396168600933,
#             0.000440082914531049)
#par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

## manage output
psi.mle2 <- sigex.eta2psi(fit.mle2[[1]]$par,constraint)
hess2 <- fit.mle2[[1]]$hessian
par.mle2 <- fit.mle2[[2]]

##  model checking

# residual analysis
resid.mle2 <- sigex.resid(psi.mle2,mdl2,data.ts)[[1]]
resid.mle2 <- sigex.load(t(resid.mle2),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
resid.acf2 <- acf(resid.mle2,lag.max=4*period,plot=FALSE)$acf
sigex.portmanteau(resid.mle2,4*period,length(psi.mle2))
sigex.gausscheck(resid.mle2)
acf(resid.mle2,lag.max=4*period)

## check on standard errors and condition numbers
print(eigen(hess2)$values)
tstats2 <- sigex.tstats(mdl2,psi.mle2,hess2,constraint)
print(tstats2)
taus2 <- log(sigex.conditions(data.ts,psi.mle2,mdl2))
print(taus2)

## bundle
analysis.mle2 <- sigex.bundle(data.ts,transform,mdl2,psi.mle2)


## model comparison

test.glr <- sigex.glr(data.ts,psi.mle2,psi.mle,mdl2,mdl)
1-pchisq(test.glr[1],df=test.glr[2])

#HERE





###########################################
### Part VI: Signal Extraction based on MLE

## load up the MLE fit for signal extraction
data.ts <- analysis.mle[[1]]
mdl <- analysis.mle[[3]]
psi <- analysis.mle[[4]]
param <- sigex.psi2par(psi,mdl,data.ts)

## get signal filters
signal.trend <- sigex.signal(data.ts,param,mdl,1)
signal.irr <- sigex.signal(data.ts,param,mdl,2)

## get extractions
extract.trend <- sigex.extract(data.ts,signal.trend,mdl,param)
extract.irr <- sigex.extract(data.ts,signal.irr,mdl,param)

## get fixed effects
reg.trend <- NULL
for(i in 1:N) {
reg.trend <- cbind(reg.trend,sigex.fixed(data.ts,mdl,i,param,"Trend")) }

## plotting
trendcol <- "tomato"
cyccol <- "navyblue"
fade <- 40
par(mfrow=c(2,1),mar=c(5,4,4,5)+.1)
for(i in 1:N)
{
	plot(data.ts[,i],xlab="Year",ylab="Trend",ylim=c(min(data.ts[,i]),max(data.ts[,i])),
		lwd=2,yaxt="n",xaxt="n")
	sigex.graph(extract.trend,reg.trend,begin.date,period,i,0,trendcol,fade)
	axis(1,cex.axis=1)
	axis(2,cex.axis=1)
}
dev.off()

## transfer function analysis
grid <- 200
frf.trend <- sigex.getfrf(data.ts,param,mdl,1,TRUE,grid)

## filter analysis
len <- 50
wk.trend <- sigex.wk(data.ts,param,mdl,1,TRUE,grid,len)


