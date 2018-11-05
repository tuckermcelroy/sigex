###########################
#### Script for Petrol Data
###########################

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
for(i in subseries)
{
	sigex.specar(data.ts,FALSE,i,period)
}
dev.off()

## growth rates
par(mfrow=c(2,1))
for(i in subseries)
{
	sigex.specar(data.ts,TRUE,i,period)
}
dev.off()


###############################
### Part III: Model Declaration

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]


################################
## Default Model: Related Trends
  
## model construction
mdl <- NULL
# trend:
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,c(1,-1))		
# irregular:
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,1)	
# regressors:
mdl <- sigex.meaninit(mdl,data.ts,0)				
 
## parameter initialization and checks
par.default <- sigex.default(mdl,data.ts)[[1]]
flag.default <- sigex.default(mdl,data.ts)[[2]]
psi.default <- sigex.par2psi(par.default,flag.default,mdl)
resid.init <- sigex.resid(psi.default,mdl,data.ts)
resid.init <- sigex.load(t(resid.init),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
acf(resid.init,lag.max=40)


###################################################
### Part IV: MLE Estimation of Related Trends Model

## setup, and fix parameters as desired
mdl.mle <- mdl
psi.mle <- psi.default
flag.mle <- Im(psi.mle)
par.mle <- sigex.psi2par(psi.mle,mdl.mle,data.ts)

## run fitting
fit.mle <- sigex.mlefit(data.ts,par.mle,flag.mle,mdl.mle,"bfgs")

## manage output
psi.mle[flag.mle==1] <- fit.mle[[1]]$par 
psi.mle <- psi.mle + 1i*flag.mle
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]

##  model checking 
resid.mle <- sigex.resid(psi.mle,mdl.mle,data.ts)
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
sigex.portmanteau(resid.mle,4*period,length(psi.mle))
sigex.gausscheck(resid.mle)
acf(resid.mle,lag.max=40)
      
## check on standard errors and condition numbers
print(eigen(hess)$values)
taus <- log(sigex.conditions(data.ts,psi.mle,mdl.mle))
print(taus)
tstats <- sigex.tstats(mdl.mle,psi.mle,hess)
stderrs <- sigex.psi2par(tstats,mdl,data.ts)
print(tstats)

## bundle  
analysis.mle <- sigex.bundle(data.ts,transform,mdl.mle,psi.mle)



#################################################
### Part V: MLE Estimation of Common Trends Model

mdl.mle2 <- mdl.mle
mdl.mle2[[1]][[1]] <- 1

par.mle2 <- sigex.default(mdl.mle2,data.ts)[[1]]
flag.ml2 <- sigex.default(mdl.mle2,data.ts)[[2]]
psi.mle2 <- sigex.par2psi(par.mle2,flag.mle2,mdl.mle2)

# run fitting
fit.mle2 <- sigex.mlefit(data.ts,par.mle2,flag.mle2,mdl.mle2,"bfgs")

## manage output
psi.mle2[flag.mle2==1] <- fit.mle2[[1]]$par 
psi.mle2 <- psi.mle2 + 1i*flag.mle2
hess2 <- fit.mle2[[1]]$hessian
par.mle2 <- fit.mle2[[2]]

##  model checking 
resid.mle2 <- sigex.resid(psi.mle2,mdl.mle2,data.ts)
resid.mle2 <- sigex.load(t(resid.mle2),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
sigex.portmanteau(resid.mle2,4*period,length(psi.mle2))
sigex.gausscheck(resid.mle2)
acf(resid.mle2,lag.max=40)
 
## check on standard errors and condition numbers
print(eigen(hess2)$values)
taus2 <- log(sigex.conditions(data.ts,psi.mle2,mdl.mle2))
print(taus2)
tstats <- sigex.tstats(mdl.mle2,psi.mle2,hess2)
stderrs <- sigex.psi2par(tstats,mdl.mle2,data.ts)
print(tstats)

# model comparison
test.glr <- sigex.glr(data.ts,psi.mle2,psi.mle,mdl.mle2,mdl.mle)
1-pchisq(test.glr[1],df=test.glr[2])

# bundle 
analysis.mle2 <- sigex.bundle(data.ts,transform,mdl.mle2,psi.mle2)
 


##########################################
### Part V: Signal Extraction based on MLE

## load up the MLE fit for signal extraction
data.ts <- analysis.mle[[1]]
mdl <- analysis.mle[[3]]
psi <- analysis.mle[[4]]
param <- sigex.psi2par(psi,mdl,data.ts)

## get signal filters
signal.trend <- sigex.signal(data.ts,param,mdl,1)
signal.irr <- sigex.signal(data.ts,param,mdl,1)

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
	plot(data.ts[,i],xlab="Year",ylab="Trend",ylim=c(min(data.ts[,i])-2,max(data.ts[,i])),
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

 
