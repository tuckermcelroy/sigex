#######################################
###  Script for Daily Immigration Data
#########################################

## wipe
rm(list=ls())

library(devtools)

setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex")
load_all(".")
 
library(Rcpp)
setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex\\R")
sourceCpp('mvar_midcast.cpp')

# testing of rcpp
x.acf <- array(0,c(2,3,2))
z <- matrix(1i,2,3)
delta <- c(1,-2,1)
debug <- TRUE
mvar_midcast(x.acf,z,delta,debug)

######################
### Part I: load data
 
# automatic: raw data 

# processing
 
n.months <- dim(imm)[1]/32
imm <- imm[-seq(1,n.months)*32,]	# strip out every 32nd row (totals)
imm <- matrix(imm[imm != 0],ncol=6) # strip out 31st False days



#############################################################
### Part II: Metadata Specifications and Exploratory Analysis

start.date <- c(9,1,1997)
end.date <- day2date(dim(imm)[1]-1,start.date)
#end.date <- c(7,31,2012)
period <- 365

# calendar calculations
start.day <- date2day(start.date[1],start.date[2],start.date[3])
end.day <- date2day(end.date[1],end.date[2],end.date[3])
begin <- c(start.date[3],start.day) 
end <- c(end.date[3],end.day)

## create ts object and plot
dataALL.ts <- sigex.load(imm,begin,period,
	c("NZArr","NZDep","VisArr","VisDep","PLTArr","PLTDep"),TRUE)
 
   
#############################
## select span and transforms 

## all data with no transform
transform <- "log"
aggregate <- FALSE
subseries <- 5
range <- NULL
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)

#######################
## spectral exploratory

## levels
par(mfrow=c(1,1))
for(i in 1:length(subseries))
{
	sigex.specar(data.ts,FALSE,i,7)
}
dev.off()

## growth rates
par(mfrow=c(1,1))
for(i in 1:length(subseries))
{
	sigex.specar(data.ts,TRUE,i,7)
}
dev.off()



#########################################
########  Modeling

transform = "log"
agg <- FALSE	# set TRUE to aggregate
series <- 1:6
#series <- 5

range <- seq(1,dim(imm)[1])
#range <- 1:500
data <- ts(as.matrix(sigex.transform(imm[range,series],transform,agg)),
	start=begin,frequency=period,names=colnames(imm)[series])
plot(data,xlab="Year")
x <- t(data)
N <- dim(x)[1]
T <- dim(x)[2]

## setup holiday regressors

# Easter Day
easter.reg1 <- gethol(easter.dates,0,0,start.date,end.date)
# pre-Easter
easter.reg2 <- gethol(easter.dates,8,-1,start.date,end.date)
# Cyber Monday
#cyber.reg <- gethol(cyber.dates,0,0,start.date,end.date)
# Black Friday
#black.reg <- gethol(black.dates,0,0,start.date,end.date)
# Superbowl Sunday
#super.reg <- gethol(super.dates,0,0,start.date,end.date)
# Labor Day
#labor.reg <- gethol(labor.dates,0,0,start.date,end.date)
# Chinese New Year
#cny.reg <- gethol(cny.dates,0,0,start.date,end.date)


#######################
# Default Model

# stochastic effects
delta.trend <- c(1,-1)
#delta.ann <- c(1,-2*cos(2*pi/365),1)

mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"trend-cycle",delta.trend)		
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"first weekly seasonal",c(1,-2*cos(2*pi/7),1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"second weekly seasonal",c(1,-2*cos(4*pi/7),1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"third weekly seasonal",c(1,-2*cos(6*pi/7),1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"irregular",1)

# fixed effects

mdl <- sigex.meaninit(mdl,data,0)
for(i in 1:N) {
mdl <- sigex.reg(mdl,i,ts(as.matrix(easter.reg1[range]),
	start=begin,frequency=period,names="Easter Day"))
mdl <- sigex.reg(mdl,i,ts(as.matrix(easter.reg2[range]),
	start=begin,frequency=period,names="pre-Easter"))
#mdl <- sigex.reg(mdl,i,ts(as.matrix(cyber.reg[range]),
#	start=begin,frequency=period,names="Cyber"))
#mdl <- sigex.reg(mdl,i,ts(as.matrix(black.reg[range]),
#	start=begin,frequency=period,names="Black"))
#mdl <- sigex.reg(mdl,i,ts(as.matrix(super.reg[range]),
#	start=begin,frequency=period,names="Super"))
#mdl <- sigex.reg(mdl,i,ts(as.matrix(labor.reg[range]),
#	start=begin,frequency=period,names="Labor"))
#mdl <- sigex.reg(mdl,i,ts(as.matrix(cny.reg[range]),
#	start=begin,frequency=period,names="CNY"))
#mdl <- sigex.reg(mdl,i,ts(as.matrix(ao.reg[,1]),
#	start=begin,frequency=period,names="AO1"))
#mdl <- sigex.reg(mdl,i,ts(as.matrix(ao.reg[,2]),
#	start=begin,frequency=period,names="AO2"))
#mdl <- sigex.reg(mdl,i,ts(as.matrix(ao.reg[,3]),
#	start=begin,frequency=period,names="AO3"))
}


par.default <- sigex.default(mdl,data)[[1]]
flag.default <- sigex.default(mdl,data)[[2]]
psi.default <- sigex.par2psi(par.default,flag.default,mdl)
#sigex.psi2par(psi.default,mdl,data)	# check
#sigex.lik(psi.default,mdl,data)
#resid.init <- sigex.resid(psi.default,mdl,data)[[1]]
#acf(t(resid.init),lag.max=40)



######################################
# MOM estimation and reduced specification

# fitting

mdl.mom <- mdl
par.mom <- sigex.momfit(data,par.default,mdl.mom)
psi.mom <- sigex.par2psi(par.mom,flag.default,mdl.mom)
resid.mom <- sigex.resid(psi.mom,mdl.mom,data)[[1]]

#thresh <- -6.22
#thresh <- -3.92
thresh <- -1.66

if(N > 1) {
reduced.mom <- sigex.reduce(data,par.mom,flag.default,mdl.mom,thresh,FALSE)
mdl.mom <- reduced.mom[[1]]
par.mom <- reduced.mom[[2]]
flag.mom <- sigex.default(mdl.mom,data)[[2]]
psi.mom <- sigex.par2psi(par.mom,flag.mom,mdl.mom)
#resid.mom <- sigex.resid(psi.mom,mdl.mom,data)[[1]]
}

log(sigex.conditions(data,psi.mom,mdl.mom))


# model checking

sigex.portmanteau(t(resid.mom),48,length(psi.mom))
sigex.gausscheck(t(resid.mom[,1:4000]))
#acf(t(resid.mom),lag.max=40)

# bundle for default span
analysis.mom <- sigex.bundle(data,transform,mdl.mom,psi.mom)

##########################################################################
################ MLE Routines (not recommended for these series) ##############
#########################################################

# MLE estimation of full model

# setup, and fix parameters as desired
mdl.mle <- mdl
psi.mle <- psi.default
flag.mle <- Im(psi.mle)
par.mle <- sigex.psi2par(psi.mle,mdl.mle,data)

# run fitting
fit.mle <- sigex.mlefit(data,par.mle,flag.mle,mdl.mle,"bfgs")

psi.mle[flag.mle==1] <- fit.mle[[1]]$par
psi.mle <- psi.mle + 1i*flag.mle
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]
resid.mle <- sigex.resid(psi.mle,mdl.mle,data)[[1]]
acf(t(resid.mle),lag.max=40)

print(eigen(hess)$values)
taus <- log(sigex.conditions(data,psi.mle,mdl.mle))
print(taus)

tstats <- sigex.tstats(mdl.mle,psi.mle,hess)
stderrs <- sigex.psi2par(tstats,mdl,data)
print(tstats)

# bundle for default span
analysis.mle <- sigex.bundle(data,transform,mdl.mle,psi.mle)


##############################
# MLE estimation of reduced model

thresh <- -6.22
#thresh <- -3.92
mdl.mle2 <- mdl.mle
par.mle2 <- par.mle
psi.mle2 <- psi.mle

if((N > 1) && (min(taus) < thresh)) {

mdl.mle2 <- sigex.reduce(data,par.mle,flag.mle,mdl.mle,thresh)[[1]]
par.mle2 <- sigex.reduce(data,par.mle,flag.mle,mdl.mle,thresh)[[2]]
flag.mle2 <- sigex.default(mdl.mle2,data)[[2]]
psi.mle2 <- sigex.par2psi(par.mle2,flag.mle2,mdl.mle2)

# run fitting
fit.mle2 <- sigex.mlefit(data,par.mle2,flag.mle2,mdl.mle2,"bfgs")

psi.mle2[flag.mle2==1] <- fit.mle2[[1]]$par
psi.mle2 <- psi.mle2 + 1i*flag.mle2
hess2 <- fit.mle2[[1]]$hessian
par.mle2 <- fit.mle2[[2]]
resid.mle2 <- sigex.resid(psi.mle2,mdl.mle2,data)[[1]]
acf(t(resid.mle2),lag.max=40)

print(eigen(hess)$values)
test.glr <- sigex.glr(data,psi.mle2,psi.mle,mdl.mle2,mdl.mle)
pchisq(test.glr[1],df=test.glr[2])

tstats <- sigex.tstats(mdl.mle2,psi.mle2,hess2)
stderrs <- sigex.psi2par(tstats,mdl.mle2,data)
print(tstats)

# bundle for default span
analysis.mle <- sigex.bundle(data,transform,mdl.mle2,psi.mle2)

}


####################################################
######################## Begin Signal Extraction

### MLE (skip)
#data <- analysis.mle[[1]]
#mdl <- analysis.mle[[3]]
#psi <- analysis.mle[[4]]
#param <- sigex.psi2par(psi,mdl,data)

## Rough: reduced MOM model
data <- analysis.mom[[1]]
mdl <- analysis.mom[[3]]
psi <- analysis.mom[[4]]
param <- sigex.psi2par(psi,mdl,data)

#######################################################################
########################## METHOD 1: DIRECT MATRIX APPROACH ############
############################### SKIP ###############################

signal.trendann <- sigex.signal(data,param,mdl,1)
signal.seas.week1 <- sigex.signal(data,param,mdl,2)
signal.seas.week2 <- sigex.signal(data,param,mdl,3)
signal.seas.week3 <- sigex.signal(data,param,mdl,4)
signal.seas.week <- sigex.signal(data,param,mdl,c(2,3,4))
signal.sa <- sigex.signal(data,param,mdl,c(1,5))

extract.trendann <- sigex.extract(data,signal.trendann,mdl,param)
extract.seas.week1 <- sigex.extract(data,signal.seas.week1,mdl,param)
extract.seas.week2 <- sigex.extract(data,signal.seas.week2,mdl,param)
extract.seas.week3 <- sigex.extract(data,signal.seas.week3,mdl,param)
extract.seas.week <- sigex.extract(data,signal.seas.week,mdl,param)
extract.sa <- sigex.extract(data,signal.sa,mdl,param)

####################################################################
############################# METHOD 2: FORECASTING and WK SIGEX
############## RECOMMENDED METHOD #####################

#grid <- 70000	# high accuracy, close to method 1
#grid <- 700	# low accuracy, but pretty fast
grid <- 7000	# need grid > filter length
window <- 200
horizon <- 2000
#leads <- c(-rev(seq(0,window-1)),seq(1,T),seq(T+1,T+window))
#data.ext <- t(sigex.cast(psi,mdl,data,leads,TRUE))
target <- array(diag(N),c(N,N,1))
 
extract.trendann <- sigex.wkextract2(psi,mdl,data,1,target,grid,window,horizon,FALSE)
extract.seas.week <- sigex.wkextract2(psi,mdl,data,c(2,3,4),target,grid,window,horizon,FALSE)
extract.seas.week1 <- sigex.wkextract2(psi,mdl,data,2,target,grid,window,horizon,NULL,FALSE)
extract.seas.week2 <- sigex.wkextract2(psi,mdl,data,3,target,grid,window,horizon,NULL,FALSE)
extract.seas.week3 <- sigex.wkextract2(psi,mdl,data,4,target,grid,window,horizon,NULL,FALSE)
extract.sa <- sigex.wkextract2(psi,mdl,data,c(1,5),target,grid,window,horizon,NULL,FALSE)
extract.irr <- sigex.wkextract2(psi,mdl,data,5,target,grid,window,horizon,NULL,FALSE)


##################################################################
#################### LP splitting of trend and cycle #################

cutoff <- pi/365
trunc <- 50000	# appropriate for mu = pi/(365)

extract.trend <- sigex.lpfiltering(mdl,data,1,NULL,psi,cutoff,grid,window,trunc,TRUE)
extract.seas.ann <- sigex.lpfiltering(mdl,data,1,NULL,psi,cutoff,grid,window,trunc,FALSE)
extract.trendirreg <- sigex.lpfiltering(mdl,data,1,5,psi,cutoff,grid,window,trunc,TRUE)
 


#########################################
### get fixed effects

reg.trend <- NULL


################################# PLOTS

## time series Plots

trendcol <- "tomato"
cyccol <- "orchid"
seascol <- "seagreen"
sacol <- "navyblue"
fade <- 60

subseries <- 1

plot(data[,subseries],xlab="Year",ylab="",ylim=c(2,9),lwd=1)
#sigex.graph(extract.sa,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trendirreg,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trend,reg.trend,begin,period,subseries,0,trendcol,10)
sigex.graph(extract.seas.ann,NULL,begin,period,subseries,5,seascol,10)
sigex.graph(extract.seas.week,NULL,begin,period,subseries,3,cyccol,fade)

plot(data[,subseries],xlab="Year",ylab="",ylim=c(0,9),lwd=1)
sigex.graph(extract.seas.week1,NULL,begin,period,subseries,5,cyccol,fade)
sigex.graph(extract.seas.week2,NULL,begin,period,subseries,3,cyccol,fade)
sigex.graph(extract.seas.week3,NULL,begin,period,subseries,1,cyccol,fade)


subseries <- 2

plot(data[,subseries],xlab="Year",ylab="",ylim=c(2,9),lwd=1)
#sigex.graph(extract.sa,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trendirreg,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trend,reg.trend,begin,period,subseries,0,trendcol,10)
sigex.graph(extract.seas.ann,NULL,begin,period,subseries,5,seascol,10)
sigex.graph(extract.seas.week,NULL,begin,period,subseries,3,cyccol,fade)

plot(data[,subseries],xlab="Year",ylab="",ylim=c(0,9),lwd=1)
sigex.graph(extract.seas.week1,NULL,begin,period,subseries,5,cyccol,fade)
sigex.graph(extract.seas.week2,NULL,begin,period,subseries,3,cyccol,fade)
sigex.graph(extract.seas.week3,NULL,begin,period,subseries,1,cyccol,fade)


subseries <- 3

plot(data[,subseries],xlab="Year",ylab="",ylim=c(2,9),lwd=1)
#sigex.graph(extract.sa,reg.trend,begin,period,subseries,0,sacol)
sigex.graph(extract.trendirreg,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trend,reg.trend,begin,period,subseries,0,trendcol,10)
sigex.graph(extract.seas.ann,NULL,begin,period,subseries,5,seascol,10)
sigex.graph(extract.seas.week,NULL,begin,period,subseries,3,cyccol,fade)

plot(data[,subseries],xlab="Year",ylab="",ylim=c(0,9),lwd=1)
sigex.graph(extract.seas.week1,NULL,begin,period,subseries,5,cyccol,fade)
sigex.graph(extract.seas.week2,NULL,begin,period,subseries,3,cyccol,fade)
sigex.graph(extract.seas.week3,NULL,begin,period,subseries,1,cyccol,fade)


subseries <- 4

plot(data[,subseries],xlab="Year",ylab="",ylim=c(2,9),lwd=1)
#sigex.graph(extract.sa,reg.trend,begin,period,subseries,0,sacol)
sigex.graph(extract.trendirreg,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trend,reg.trend,begin,period,subseries,0,trendcol,10)
sigex.graph(extract.seas.ann,NULL,begin,period,subseries,5,seascol,10)
sigex.graph(extract.seas.week,NULL,begin,period,subseries,3,cyccol,fade)

plot(data[,subseries],xlab="Year",ylab="",ylim=c(0,9),lwd=1)
sigex.graph(extract.seas.week1,NULL,begin,period,subseries,5,cyccol,fade)
sigex.graph(extract.seas.week2,NULL,begin,period,subseries,3,cyccol,fade)
sigex.graph(extract.seas.week3,NULL,begin,period,subseries,1,cyccol,fade)


subseries <- 5

plot(data[,subseries],xlab="Year",ylab="",ylim=c(2,7),lwd=1)
#sigex.graph(extract.sa,reg.trend,begin,period,subseries,0,sacol)
sigex.graph(extract.trendirreg,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trend,reg.trend,begin,period,subseries,0,trendcol,10)
sigex.graph(extract.seas.ann,NULL,begin,period,subseries,4,seascol,10)
sigex.graph(extract.seas.week,NULL,begin,period,subseries,3,cyccol,fade)

plot(data[,subseries],xlab="Year",ylab="",ylim=c(0,7),lwd=1)
sigex.graph(extract.seas.week1,NULL,begin,period,subseries,3,cyccol,fade)
sigex.graph(extract.seas.week2,NULL,begin,period,subseries,2,cyccol,fade)
sigex.graph(extract.seas.week3,NULL,begin,period,subseries,1,cyccol,fade)


subseries <- 6

plot(data[,subseries],xlab="Year",ylab="",ylim=c(2,7),lwd=1)
#sigex.graph(extract.sa,reg.trend,begin,period,subseries,0,sacol)
sigex.graph(extract.trendirreg,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trend,reg.trend,begin,period,subseries,0,trendcol,10)
sigex.graph(extract.seas.ann,NULL,begin,period,subseries,4,seascol,10)
sigex.graph(extract.seas.week,NULL,begin,period,subseries,3,cyccol,fade)

plot(data[,subseries],xlab="Year",ylab="",ylim=c(0,7),lwd=1)
sigex.graph(extract.seas.week1,NULL,begin,period,subseries,3,cyccol,fade)
sigex.graph(extract.seas.week2,NULL,begin,period,subseries,2,cyccol,fade)
sigex.graph(extract.seas.week3,NULL,begin,period,subseries,1,cyccol,fade)


#################################
### Seasonality Diagnostics

sigex.seascheck(extract.trend[[1]],7,.04,1)
sigex.seascheck(extract.seas.ann[[1]],7,.04,1)
sigex.seascheck(extract.seas.week[[1]],7,.04,0)
sigex.seascheck(extract.trendirreg[[1]],7,.04,1)

sigex.seascheck(extract.trend[[1]],365,.04,1)
sigex.seascheck(extract.seas.ann[[1]],365,.04,1)
sigex.seascheck(extract.seas.week[[1]],365,.04,0)
sigex.seascheck(extract.trendirreg[[1]],365,.04,1)

###### Signal extraction diagnostics

round(sigex.signalcheck(extract.trendann[[1]][(horizon+1):(horizon+T),],param,mdl,1,100),digits=3)
round(sigex.signalcheck(extract.seas.week[[1]],param,mdl,c(2,3,4),100),digits=3)

### Spectral density plots

subseries <- 1

## spectral plots
week.freq <- 365/7
month.freq <- 365/(365/12)
ann.freq <- 1

spec.ar(ts(extract.trend[[1]][,subseries],frequency=period),main="")
abline(v=ann.freq,col=4)
abline(v=week.freq,col=2)
abline(v=2*week.freq,col=2)
abline(v=3*week.freq,col=2)
abline(v=month.freq,col=3)

spec.ar(ts(extract.sa[[1]][,subseries],frequency=period),main="")
abline(v=ann.freq,col=4)
abline(v=week.freq,col=2)
abline(v=2*week.freq,col=2)
abline(v=3*week.freq,col=2)
abline(v=month.freq,col=3)

spec.ar(ts(extract.seas.week1[[1]][,subseries],frequency=period))
abline(v=week.freq,col=2)

spec.ar(ts(extract.seas.week2[[1]][,subseries],frequency=period))
abline(v=2*week.freq,col=2)

spec.ar(ts(extract.seas.week3[[1]][,subseries],frequency=period))
abline(v=3*week.freq,col=2)

##################################
#	filter analysis

grid <- 700
frf.trendann <- sigex.getfrf(data,param,mdl,1,TRUE,grid)
frf.week1 <- sigex.getfrf(data,param,mdl,2,TRUE,grid)
frf.week2 <- sigex.getfrf(data,param,mdl,3,TRUE,grid)
frf.week3 <- sigex.getfrf(data,param,mdl,4,TRUE,grid)
frf.weeks <- sigex.getfrf(data,param,mdl,c(2,3,4),TRUE,grid)
frf.irr <- sigex.getfrf(data,param,mdl,5,TRUE,grid)
frf.sa <- sigex.getfrf(data,param,mdl,c(1,5),TRUE,grid)

len <- 500
wk.trendann <- sigex.getwk(data,param,mdl,1,TRUE,grid,len)
wk.week1 <- sigex.getwk(data,param,mdl,2,TRUE,grid,len)
wk.week2 <- sigex.getwk(data,param,mdl,3,TRUE,grid,len)
wk.week3 <- sigex.getwk(data,param,mdl,4,TRUE,grid,len)
wk.weeks <- sigex.getwk(data,param,mdl,c(2,3,4),TRUE,grid,len)
wk.irr <- sigex.getwk(data,param,mdl,6,TRUE,grid,len)







##### SCRAP

