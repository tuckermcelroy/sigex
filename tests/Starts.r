#######################################
###  Script for Housing Starts Data
#########################################

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

start.date = c(1964,1)
period <- 12

## create ts object and plot
dataALL.ts <- sigex.load(starts,start.date,period,c("South","West","NE","MW"),TRUE)

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


#############################
### Part IV: Model Estimation

##############
## MOM fitting
mdl.mom <- mdl
constraint <- NULL
par.default <- sigex.default(mdl.mom,data.ts,constraint)
par.mom <- sigex.momfit(data.ts,par.default,mdl.mom)
psi.mom <- sigex.par2psi(par.mom,mdl.mom)
resid.mom <- sigex.resid(psi.mom,mdl.mom,data.ts)[[1]]
resid.mom <- sigex.load(t(resid.mom),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)

## compute a reduced rank model
thresh <- -6.22
reduced.mom <- sigex.reduce(data.ts,par.mom,mdl.mom,thresh,FALSE)
mdl.mom <- reduced.mom[[1]]
par.mom <- reduced.mom[[2]]
psi.mom <- sigex.par2psi(par.mom,mdl.mom)
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


#############
# MLE Fitting

#  Initialize with MOM estimates
constraint <- NULL
psi.mle <- sigex.par2psi(par.mom,mdl)

## run fitting: can be commented out, this takes a while
fit.mle <- sigex.mlefit(data.ts,par.mom,constraint,mdl,"bfgs",debug=TRUE)

## input parameter from previous fit (MLE on entire span)
#  divergence:  6789.661
psi.mle <- c(0.49358609305695, 0.17848725859254, 0.34121739912571, 0.39917727415425,
             0.84832530464264, 0.68306879252262, -2.3494687111314, -5.47534663726587,
             -6.69385117951384, -6.08364145983965, 0.87510015081027, 0.22197127114861,
             0.50086675920103, 0.3406250169841, 0.7910378054958, 0.98544026276858,
             -2.52890913740106, -4.29524634814519, -5.98519527750281, -4.88659954275053,
             0.09574663273149, 0.20131335062649, 0.8493518091576, 0.48420520104336,
             0.62643997675928, 1.13945063379914, -4.04217214895869, -4.68919816059416,
             -4.73313805629826, -4.0627015759002, 0.9234957516084, -0.39606729445073,
             0.24466504619404, -0.36570474542918, 0.36399571873663, 0.75871517273776,
             -3.05567431351817, -4.74337970092605, -4.96364133429136, -5.06144086942249,
             0.26296368360579, -0.18159940066192, 0.14979583325899, -0.10599164910036,
             0.21503766242974, -0.14164986104397, -2.07489346121933, -3.64302004053168,
             -5.69277788172285, -5.3689470753418, 1.40718934367933, -0.00854528787477,
             -0.21988633727394, 0.0283662345071, 1.23786259577472, 0.19983413521575,
             -4.53336362894347, -4.70016052568401, -7.07530853221777, -6.03054443735399,
             -0.09955060405249, 0.11660784869795, 0.15789980223364, -0.03631849815476,
             0.18385749297074, 0.32935147758533, -2.1377604820296, -3.62882764786239,
             -5.11279846492415, -3.62475631527416, 0.12430528614515, 0.02925079204219,
             -0.08733491948454, 0.17897776431614, 0.48438912873225, 0.26583597642199,
             1.87566939226944, 0.1445002084775, -1.34264222816582, -0.30536763401493,
             -0.00488431480035, -0.00094565956468, -0.00106126820173, -0.00041365883889
)

# bundle for default span
analysis.mle <- sigex.bundle(data.ts,transform,mdl,psi.mle)




##########################################
### Part V: Signal Extraction based on fit

## load up the MOM fit for signal extraction
data.ts <- analysis.mom[[1]]
mdl <- analysis.mom[[3]]
psi <- analysis.mom[[4]]
param <- sigex.psi2par(psi,mdl,data.ts)


#######################################################################
########################## METHOD 1: DIRECT MATRIX APPROACH ############

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
#pdf(file="StartsSignals.pdf")
par(mfrow=c(2,2))
for(i in 1:N)
{
  plot(data.ts[,i],xlab="Year",ylab="",ylim=c(min(data.ts[,i])-25,max(data.ts[,i])),lwd=1)
  sigex.graph(extract.sa,reg.trend,begin.date,period,i,0,sacol,fade)
  sigex.graph(extract.trend,reg.trend,begin.date,period,i,0,trendcol,fade)
  sigex.graph(extract.seas,NULL,begin.date,period,i,min(data.ts[,i])-10,seascol,fade)
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


####################################################################
############################# METHOD 2: FORECASTING and WK SIGEX

#grid <- 70000	# high accuracy, close to method 1
#grid <- 700	# low accuracy, but pretty fast
grid <- 7000	# need grid > filter length
window <- 10
#window <- 50
horizon <- 0
target <- array(diag(N),c(N,N,1))

extract.trend2 <- sigex.wkextract(psi,mdl,data.ts,1,target,grid,window,horizon,TRUE)
extract.seas2 <- sigex.wkextract(psi,mdl,data.ts,seq(2,7),target,grid,window,horizon,TRUE)
extract.sa2 <- sigex.wkextract(psi,mdl,data.ts,c(1,8),target,grid,window,horizon,TRUE)


## root mse plots: trend
#pdf(file="StartsTrendMSE50.pdf")
#pdf(file="StartsTrendMSE10.pdf")
par(mfrow=c(2,2))
for(subseries in 1:N)
{
  ex.true <- (extract.trend[[1]][,subseries]-extract.trend[[3]][,subseries])/2
  ex.approx <- (extract.trend2[[1]][,subseries]-extract.trend2[[3]][,subseries])/2
  plot(ts(ex.true,start=begin.date,frequency=period),ylab="Root MSE",xlab="Year",ylim=c(0,max(ex.true)))
  lines(ts(ex.approx,start=begin.date,frequency=period),col=trendcol)
}
#dev.off()

## root mse plots: sa
#pdf(file="StartsSAMSE50.pdf")
#pdf(file="StartsSAMSE10.pdf")
par(mfrow=c(2,2))
for(subseries in 1:N)
{
  ex.true <- (extract.sa[[1]][,subseries]-extract.sa[[3]][,subseries])/2
  ex.approx <- (extract.sa2[[1]][,subseries]-extract.sa2[[3]][,subseries])/2
  plot(ts(ex.true,start=begin.date,frequency=period),ylab="Root MSE",xlab="Year",ylim=c(0,max(ex.true)))
  lines(ts(ex.approx,start=begin.date,frequency=period),col=sacol)
}
#dev.off()

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
for(k in 1:N)
{
  reg.gr[,k] <- filter(reg.trend[,k],gr.poly,method="convolution",sides=1)
}

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
#dev.off()

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
#dev.off()




