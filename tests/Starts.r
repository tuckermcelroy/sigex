#######################################
###  Script for Housing Starts Data
#########################################

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

start.date = c(1964,1)
period <- 12

## create ts object and plot
dataALL.ts <- sigex.load(starts,start.date,period,c("South","West","NE","MW"),TRUE)

#############################
## select span and transforms

## all data for NE-MW with log transform
transform <- "log"
aggregate <- FALSE
subseries <- c(3,4)
begin.date <- start(dataALL.ts)
end.date <- end(dataALL.ts)
range <- NULL
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)

## recent span with no transform
transform <- "none"
aggregate <- FALSE
subseries <- c(1,2)
begin.date <- c(2004,1)
end.date <- end(dataALL.ts)
range <- list(begin.date,end.date)
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)

#######################
## spectral exploratory

## levels
par(mfrow=c(2,2))
for(i in subseries)
{
	sigex.specar(data.ts,FALSE,i,period)
}
dev.off()

## growth rates
par(mfrow=c(2,2))
for(i in subseries)
{
	sigex.specar(data.ts,TRUE,i,period)
}
dev.off()


###############################
### Part III: Model Declaration

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]


################
## Default Model

def <- c(0,1,0,1)

## model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"trend",c(1,-2,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"first seasonal",c(1,-sqrt(3),1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"second seasonal",c(1,-1,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"third seasonal",c(1,0,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"fourth seasonal",c(1,1,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"fifth seasonal",c(1,sqrt(3),1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"sixth seasonal",c(1,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"irregular",1)
# regressors:
mdl <- sigex.meaninit(mdl,data.ts,0)

## parameter initialization and checks
par.default <- sigex.default(mdl,data.ts)[[1]]
flag.default <- sigex.default(mdl,data.ts)[[2]]
psi.default <- sigex.par2psi(par.default,flag.default,mdl)
resid.init <- sigex.resid(psi.default,mdl,data.ts)[[1]]
resid.init <- sigex.load(t(resid.init),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
acf(resid.init,lag.max=40)



###########################
### Part IV: MOM Estimation

## fitting
mdl.mom <- mdl
par.mom <- sigex.momfit(data.ts,par.default,mdl.mom)
psi.mom <- sigex.par2psi(par.mom,flag.default,mdl.mom)
resid.mom <- sigex.resid(psi.mom,mdl.mom,data.ts)[[1]]
resid.mom <- sigex.load(t(resid.mom),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)

## compute a reduced rank model
thresh <- -6.22
reduced.mom <- sigex.reduce(data.ts,par.mom,flag.default,mdl.mom,thresh,FALSE)
mdl.mom <- reduced.mom[[1]]
par.mom <- reduced.mom[[2]]
flag.mom <- sigex.default(mdl.mom,data.ts)[[2]]
psi.mom <- sigex.par2psi(par.mom,flag.mom,mdl.mom)
resid.mom <- sigex.resid(psi.mom,mdl.mom,data.ts)[[1]]
resid.mom <- sigex.load(t(resid.mom),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
acf(resid.mom,lag.max=40)

## examine condition numbers
log(sigex.conditions(data.ts,psi.mom,mdl.mom))

## model checking
sigex.portmanteau(resid.mom,4*period,length(psi.mom))
sigex.gausscheck(resid.mom)

## bundle
analysis.mom <- sigex.bundle(data.ts,transform,mdl.mom,psi.mom)


##########################################
### Part V: Signal Extraction based on MOM

## load up the MOM fit for signal extraction
data.ts <- analysis.mom[[1]]
mdl <- analysis.mom[[3]]
psi <- analysis.mom[[4]]
param <- sigex.psi2par(psi,mdl,data.ts)

## get signal filters
signal.trend <- sigex.signal(data.ts,param,mdl,1)
signal.seas <- sigex.signal(data.ts,param,mdl,seq(2,7))
signal.sa <- sigex.signal(data.ts,param,mdl,c(1,8))

## get extractions
extract.trend <- sigex.extract(data.ts,signal.trend,mdl,param)
extract.seas <- sigex.extract(data.ts,signal.seas,mdl,param)
extract.sa <- sigex.extract(data.ts,signal.sa,mdl,param)

## get fixed effects
reg.trend <- NULL
for(i in 1:N) {
reg.trend <- cbind(reg.trend,param[[4]][i]*mdl[[4]][[i]]) }

## plotting
trendcol <- "tomato"
cyccol <- "orchid"
seascol <- "seagreen"
sacol <- "navyblue"
fade <- 60
par(mfrow=c(1,2))
for(i in 1:N)
{
	plot(data.ts[,i],xlab="Year",ylab="",ylim=c(min(data.ts[,i])-20,max(data.ts[,i])),lwd=1)
	sigex.graph(extract.sa,reg.trend,begin.date,period,i,0,sacol,fade)
	sigex.graph(extract.trend,reg.trend,begin.date,period,i,0,trendcol,fade)
	sigex.graph(extract.seas,NULL,begin.date,period,i,-5,seascol,fade)
}

## spectral diagnostics: trend
par(mfrow=c(1,2))
for(i in 1:N)
{
	sigex.specar(ts(extract.trend[[1]],frequency=period,names=colnames(data.ts)),FALSE,i,period)
}

## spectral diagnostics: sa
par(mfrow=c(1,2))
for(i in 1:N)
{
	sigex.specar(ts(extract.sa[[1]],frequency=period,names=colnames(data.ts)),FALSE,i,period)
}

## transfer function analysis
grid <- 200
frf.trend <- sigex.getfrf(data.ts,param,mdl,1,TRUE,grid)
frf.seas <- sigex.getfrf(data.ts,param,mdl,seq(2,7),TRUE,grid)
frf.sa <- sigex.getfrf(data.ts,param,mdl,c(1,8),TRUE,grid)

## filter analysis
len <- 50
wk.trend <- sigex.wk(data.ts,param,mdl,1,TRUE,grid,len)
wk.seas <- sigex.wk(data.ts,param,mdl,seq(2,7),TRUE,grid,len)
wk.sa <- sigex.wk(data.ts,param,mdl,c(1,8),TRUE,grid,len)


