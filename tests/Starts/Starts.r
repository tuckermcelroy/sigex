#######################################
###  Script for Housing Starts Data
#########################################

## wipe
rm(list=ls())

library(devtools)

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
#fit.mle <- sigex.mlefit(data.ts,par.mom,constraint,mdl,"bfgs",debug=TRUE)

## input parameter from previous fit (MLE on entire span)
#  divergence:  6185.784
psi.mle <- c(0.423542257365999, 0.134501458946877, 0.225735894562736, 0.125500806127999,
             0.277059966369363, 0.840944149182132, -2.20367055605147, -4.85935947712897,
             -6.69729008511858, -6.11816927571922, 1.04269048770371, 0.250317339091603,
             0.809753277323495, -0.720313773007307, -0.616754630646598, -0.386827886912457,
             -2.40614000203074, -4.72919770852222, -6.53509436667913, -9.2861412575583,
             0.48380855905291, 0.237291064979205, 1.21270650674356, -1.3752518573378,
             -0.53671889759918, -0.0210044890364326, -3.89972527892544, -5.83039542759658,
             -11.0486991042118, -9.47623410286397, -0.172067074460752, 0.212808975179007,
             0.165690017103288, -0.00637244141186747, -0.00214477902773064,
             0.0933228648769178, -2.73819472549497, -3.6650836810548, -4.49590923175739,
             -10.0074671892842, -0.114866000887492, -0.210563103438697, -0.182426980405228,
             -0.678921798038978, 0.244592711560399, 0.000162440817641524,
             -3.29062077278078, -6.94067047408281, -13.8631958930152, -12.7740076647554,
             1.01976693974461, -0.176728058418555, 0.381748041673584, 7.76293530769596,
             -42.5224142483839, 0.000307656113372132, -4.86346606363989, -11.0656550501414,
             -12.3938515712309, -9.26711017917879, -0.245307834874516, 0.0814078360607122,
             0.156989749491884, 1.13373834073925, -0.245523810285003, -0.00181549624522899,
             -3.63022671080751, -6.84081933795189, -12.4239069581581, -11.3221468947418,
             0.0266707764314609, 0.0281723045059473, 0.000994959872918981,
             0.0659040290710181, 0.021912146700469, 0.280000953632169, 2.6066131275945,
             1.06766682710889, -0.0913293751220125, 0.0660748780411686, 0.00387853583107558,
             0.00175337864601441, 0.000505786944522091, 0.00104654719895988)
par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

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

# bundle for default span
analysis.mle <- sigex.bundle(data.ts,transform,mdl,psi.mle)



##########################################
### Part V: Signal Extraction based on fit

## load up the MOM fit for signal extraction
data.ts <- analysis.mle[[1]]
mdl <- analysis.mle[[3]]
psi <- analysis.mle[[4]]
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
for(i in 1:N) { reg.trend <- cbind(reg.trend,param[[4]][i]*mdl[[4]][[i]]) }

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




