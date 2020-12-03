#######################################
###  Script for Housing Starts Data
#######################################

## wipe
rm(list=ls())

library(devtools)
library(Rcpp)

# suppose directory is set to where sigex is located, e.g.
#setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex")
load_all(".")
root.dir <- getwd()
setwd(paste(root.dir,"/tests/Starts",sep=""))

######################
### Part I: load data

# automatic

#############################################################
### Part II: Metadata Specifications and Exploratory Analysis

start.date = c(1964,1)
period <- 12

## create ts object and plot
dataALL.ts <- sigex.load(starts,start.date,period,
                         c("South","West","NE","MW"),TRUE)

#############################
## select span and transforms

## recent span with no transform
transform <- "none"
aggregate <- FALSE
subseries <- c(1,2,3,4)
begin.date <- start.date
end.date <- end(dataALL.ts)
range <- list(begin.date,end.date)
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)

#######################
## spectral exploratory

## levels
par(mfrow=c(2,2))
for(i in subseries) {	sigex.specar(data.ts,FALSE,i,period) }
dev.off()

## growth rates
par(mfrow=c(2,2))
for(i in subseries) {	sigex.specar(data.ts,TRUE,i,period) }
dev.off()


###############################
### Part III: Model Declaration

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]

################
## Default Model

## model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"trend",c(1,-2,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"first seasonal",c(1,-sqrt(3),1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"second seasonal",c(1,-1,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"third seasonal",c(1,0,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"fourth seasonal",c(1,1,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"fifth seasonal",c(1,sqrt(3),1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"sixth seasonal",c(1,1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),NULL,"irregular",1)
# regressors:
mdl <- sigex.meaninit(mdl,data.ts,0)


#############################
### Part IV: Model Estimation

##############
## MOM fitting

mdl.mom <- mdl
constraint <- NULL
par.default <- sigex.default(mdl.mom,data.ts,constraint)
par.mom <- sigex.momfit(data.ts,par.default,mdl.mom)
psi.mom <- sigex.par2psi(par.mom,mdl.mom)
sigex.lik(psi.mom,mdl.mom,data.ts,FALSE)
# divergence  6329.107

## residual analysis
resid.mom <- sigex.resid(psi.mom,mdl.mom,data.ts)[[1]]
resid.mom <- sigex.load(t(resid.mom),start(data.ts),
                        frequency(data.ts),colnames(data.ts),TRUE)
resid.acf <- acf(resid.mom,lag.max=4*period,plot=TRUE)$acf

## examine condition numbers
log(sigex.conditions(data.ts,psi.mom,mdl.mom))

## compute a reduced rank model
thresh <- -6.22
reduced.mom <- sigex.reduce(data.ts,par.mom,mdl.mom,thresh,TRUE)
mdl.mom <- reduced.mom[[1]]
par.mom <- reduced.mom[[2]]
psi.mom <- sigex.par2psi(par.mom,mdl.mom)
sigex.lik(psi.mom,mdl.mom,data.ts,FALSE)
# div  6292.841

## residual analysis
resid.mom <- sigex.resid(psi.mom,mdl.mom,data.ts)[[1]]
resid.mom <- sigex.load(t(resid.mom),start(data.ts),
                        frequency(data.ts),colnames(data.ts),TRUE)
resid.acf <- acf(resid.mom,lag.max=4*period,plot=TRUE)$acf

## examine condition numbers
log(sigex.conditions(data.ts,psi.mom,mdl.mom))

## model checking
sigex.portmanteau(resid.mom,4*period,length(psi.mom))
sigex.gausscheck(resid.mom)

## bundle
analysis.mom <- sigex.bundle(data.ts,transform,mdl.mom,psi.mom)


#############
# MLE Fitting

## load up the MOM model
data.ts <- analysis.mom[[1]]
mdl <- analysis.mom[[3]]
psi.mom <- analysis.mom[[4]]
par.mom <- sigex.psi2par(psi.mom,mdl,data.ts)

#  Initialize with MOM estimates
constraint <- NULL
psi.mle <- sigex.par2psi(par.mom,mdl)

## run fitting: can be commented out, this takes a while
#fit.mle <- sigex.mlefit(data.ts,par.mom,constraint,mdl,"bfgs",debug=TRUE)

## input parameter from previous fit (MLE on entire span)
#  divergence:  6152.607
psi.mle <- c(0.430604758988307, 0.138843971938441, 0.211034753655249, 0.122497565015475,
             0.273961955347927, 0.83829841362335, -2.19718581249183, -4.85906104866881,
             -6.69696771334971, -6.13282806492413, 1.12175825784647, 0.289882305425559,
             0.873677616319397, -0.745293003402034, -0.754732471766407, -0.498376140705462,
             -2.42443529609759, -4.72941611158525, -6.52941280216317, 0.459400674188665,
             0.208306224143506, 1.24487637498583, -1.388883213613, -0.526433078980047,
             -3.90916616951099, -5.82829964996699, -0.169794036259361, 0.22214398081366,
             0.152889428361731, -0.00908575947899202, -0.00318543112970812,
             0.0884748269673841, -2.73579660752933, -3.66216778737483, -4.49261642825817,
             -0.112193904428172, -0.203557855469594, -0.164528907998198, -0.680285482020193,
             0.244662042202419, -3.28673706196714, -6.93954034182091, 1.05554783833862,
             -0.16328976555996, 0.311473085983565, -4.85381927837337, -0.244998596862816,
             0.0508029104468634, 0.194372016161389, 1.14475061810928, -0.253386348374559,
             -3.61493340437869, -6.83641878016894, 0.0921157226442518, 0.0283709627836627,
             -0.0389611692363272, 0.0109437749519264, 0.0358792010736619,
             0.245040148693726, 2.69734505562902, 1.07979552299194, -0.0671443528547337,
             0.302809384452077, 0.000724801053757517, 0.000779888484070904,
             -2.68663079304438e-05, 0.000373610017630947)
par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

## residual analysis
resid.mle <- sigex.resid(psi.mle,mdl,data.ts)[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),
                        frequency(data.ts),colnames(data.ts),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*period,plot=TRUE)$acf

## examine condition numbers
log(sigex.conditions(data.ts,psi.mle,mdl))

## model checking
sigex.portmanteau(resid.mle,4*period,length(psi.mle))
sigex.gausscheck(resid.mle)

# bundle for default span
analysis.mle <- sigex.bundle(data.ts,transform,mdl,psi.mle)



##########################################
### Part V: Signal Extraction based on fit

## load up the MLE fit for signal extraction
data.ts <- analysis.mle[[1]]
mdl <- analysis.mle[[3]]
psi <- analysis.mle[[4]]
param <- sigex.psi2par(psi,mdl,data.ts)


#######################################################################
######################### METHOD 1: DIRECT MATRIX APPROACH ############

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
  reg.trend <- cbind(reg.trend,sigex.fixed(data.ts,mdl,i,param,"Trend")) }

## plotting
trendcol <- "tomato"
seascol <- "seagreen"
sacol <- "navyblue"
fade <- 60
#pdf(file="StartsSignals.pdf")
par(mfrow=c(2,2))
for(i in 1:N)
{
  plot(data.ts[,i],xlab="Year",ylab="",ylim=c(min(data.ts[,i])-25,max(data.ts[,i])),lwd=1)
  sigex.graph(extract.sa,reg.trend,begin.date,period,i,0,sacol,fade)
  sigex.graph(extract.trend,reg.trend,begin.date,period,i,0,trendcol,fade)
  sigex.graph(extract.seas,NULL,begin.date,period,i,min(data.ts[,i])-10,seascol,fade)
}
dev.off()

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
dev.off()
#pdf(file="StartsSeasfrf.pdf")
frf.seas <- sigex.getfrf(data.ts,param,mdl,seq(2,7),TRUE,grid)
dev.off()
#pdf(file="StartsSAfrf.pdf")
frf.sa <- sigex.getfrf(data.ts,param,mdl,c(1,8),TRUE,grid)
dev.off()

## filter analysis
len <- 50
target <- array(diag(N),c(N,N,1))
#pdf(file="StartsTrendwk.pdf")
wk.trend <- sigex.wk(data.ts,param,mdl,1,target,TRUE,grid,len)
dev.off()
#pdf(file="StartsSeaswk.pdf")
wk.seas <- sigex.wk(data.ts,param,mdl,seq(2,7),target,TRUE,grid,len)
dev.off()
#pdf(file="StartsSAwk.pdf")
wk.sa <- sigex.wk(data.ts,param,mdl,c(1,8),target,TRUE,grid,len)
dev.off()


####################################################################
############################# METHOD 2: FORECASTING and WK SIGEX

grid <- 7000
window <- 50
horizon <- 0
target <- array(diag(N),c(N,N,1))

extract.trend2 <- sigex.wkextract(psi,mdl,data.ts,1,target,grid,window,horizon,TRUE)
extract.seas2 <- sigex.wkextract(psi,mdl,data.ts,seq(2,7),target,grid,window,horizon,TRUE)
extract.sa2 <- sigex.wkextract(psi,mdl,data.ts,c(1,8),target,grid,window,horizon,TRUE)

## root mse plots: trend
#pdf(file="StartsTrendMSE50.pdf")
par(mfrow=c(2,2))
for(subseries in 1:N)
{
  ex.true <- (extract.trend[[1]][,subseries]-extract.trend[[3]][,subseries])/2
  ex.approx <- (extract.trend2[[1]][,subseries]-extract.trend2[[3]][,subseries])/2
  plot(ts(ex.true,start=begin.date,frequency=period),ylab="Root MSE",xlab="Year",ylim=c(0,max(ex.true)))
  lines(ts(ex.approx,start=begin.date,frequency=period),col=trendcol)
}
dev.off()

## root mse plots: sa
#pdf(file="StartsSAMSE50.pdf")
par(mfrow=c(2,2))
for(subseries in 1:N)
{
  ex.true <- (extract.sa[[1]][,subseries]-extract.sa[[3]][,subseries])/2
  ex.approx <- (extract.sa2[[1]][,subseries]-extract.sa2[[3]][,subseries])/2
  plot(ts(ex.true,start=begin.date,frequency=period),ylab="Root MSE",xlab="Year",ylim=c(0,max(ex.true)))
  lines(ts(ex.approx,start=begin.date,frequency=period),col=sacol)
}
dev.off()


##########################################
### Part VI: Trend Growth Rate

#######################################################################
########################## METHOD 1: DIRECT MATRIX APPROACH ############

gr.poly <- c(1,-1)
p <- length(gr.poly)-1
gr.mat <- t(rev(gr.poly)) %x% diag(N)
mse.trend <- array(signal.trend[[2]],c(N,T,N,T))
mse.trend.gr <- array(0,c(N,N,T-p))
trend.gr <- array(0,c(T-p,N))
trend.est <- extract.trend[[1]] + reg.trend
for(t in (p+1):T)
{
  trend.gr[t-p,] <- gr.mat %*% matrix(t(trend.est[(t-p):t,]),ncol=1)
  mse.trend.gr[,,t-p] <- gr.mat %*% matrix(mse.trend[,(t-p):t,,(t-p):t],c(N*(p+1),N*(p+1))) %*% t(gr.mat)
}
trend.gr <- ts(rbind(array(NA,c(p,N)),trend.gr),
               start=begin.date,frequency=period,names=colnames(data.ts))
extract.trendgr <- list(trend.gr,trend.gr,trend.gr)
for(k in 1:N)
{
  extract.trendgr[[2]][,k] <- extract.trendgr[[2]][,k] + 2*sqrt(c(rep(NA,p),mse.trend.gr[k,k,]))
  extract.trendgr[[3]][,k] <- extract.trendgr[[3]][,k] - 2*sqrt(c(rep(NA,p),mse.trend.gr[k,k,]))
}

## plot
par(mfrow=c(2,2))
for(i in 1:N)
{
  plot(data.ts[,i],xlab="Year",ylab="",ylim=c(min(extract.trendgr[[3]][,i],na.rm=TRUE),
                                              max(extract.trendgr[[2]][,i],na.rm=TRUE)),lwd=1,col=0)
  sigex.graph(extract.trendgr,NULL,begin.date,period,i,0,trendcol,fade)
}

####################################################################
############################# METHOD 2: FORECASTING and WK SIGEX

grid <- 70000	  # need large grid value to get accuracy
window <- 100
horizon <- 0
gr.array <- array(t(gr.poly) %x% diag(N),c(N,N,p+1))
reg.gr <- array(0,c(T,N))
for(k in 1:N) { reg.gr[,k] <- filter(reg.trend[,k],gr.poly,method="convolution",sides=1) }

## Trend case
extract.trendgr2 <- sigex.wkextract(psi,mdl,data.ts,1,gr.array,grid,window,horizon,TRUE)
## plot
#pdf(file="StartsTrendGR.pdf")
par(mfrow=c(2,2))
for(i in 1:N)
{
  plot(data.ts[,i],xlab="Year",ylab="",ylim=c(min(extract.trendgr2[[3]][,i]+reg.gr[,i],na.rm=TRUE),
                                              max(extract.trendgr2[[2]][,i]+reg.gr[,i],na.rm=TRUE)),lwd=1,col=0)
  sigex.graph(extract.trendgr2,reg.gr,begin.date,period,i,0,trendcol,fade)
}
dev.off()

# SA case
extract.sagr2 <- sigex.wkextract(psi,mdl,data.ts,c(1,8),gr.array,grid,window,horizon,TRUE)
## plot
#pdf(file="StartsSAGR.pdf")
par(mfrow=c(2,2))
for(i in 1:N)
{
  plot(data.ts[,i],xlab="Year",ylab="",ylim=c(min(extract.sagr2[[3]][,i]+reg.gr[,i],na.rm=TRUE),
                                              max(extract.sagr2[[2]][,i]+reg.gr[,i],na.rm=TRUE)),lwd=1,col=0)
  sigex.graph(extract.sagr2,reg.gr,begin.date,period,i,0,sacol,fade)
}
dev.off()




