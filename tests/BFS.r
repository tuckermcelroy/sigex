#########################################
###  Script for Weekly BFS Data
#########################################

## wipe
rm(list=ls())

library(devtools)

setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex")
load_all(".")

###########################################
### Part I: load data and special functions

load("C:\\Users\\neide\\OneDrive\\Documents\\Research\\bfs.RData")

ubgenerator <- function(period,trunc.len,m)
{
  
  ceps2wold <- function(ceps,q)
  {
    m <- length(ceps)
    if(q > m) {	ceps <- c(ceps,rep(0,q-m)) }
    wold <- 1
    wolds <- wold
    for(j in 1:q)
    {
      wold <- sum(seq(1,j)*ceps[1:j]*wolds[j:1])/j
      wolds <- c(wolds,wold)
    }
    return(wolds)
  }
  
  half.len <- floor(period/2)
  if(length(trunc.len)==0) { trunc.len <- half.len }
  ceps <- rep(0,m)
  
  for(ell in 1:m)
  {
    ceps[ell] <- -2*sum(cos(2*pi*ell*seq(1,trunc.len)/period))/ell
  }
  wolds <- ceps2wold(ceps,2*trunc.len)
  
  return(wolds)
}



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

## define trend and SA weekly filters
week.period <- 365.25/7
half.len <- floor(week.period/2)
p.seas <- 3
trend.filter <- ubgenerator(week.period,NULL,1000)
trend.filter <- trend.filter/sum(trend.filter)
#plot.ts(trend.filter)
detrend.filter <- c(rep(0,half.len),1,rep(0,half.len)) - trend.filter
seas.filter <- 0
for(j in 1:p.seas)
{
  week.periodj <- j*week.period
  half.lenj <- floor(week.periodj/2)
  seas.filterj <- ubgenerator(week.periodj,half.lenj-1,1000)
  seas.filterj <- polymult(seas.filterj,c(1,0,-1))
  seas.filter <- c(seas.filter,rep(0,length(seas.filterj)-length(seas.filter)))
  seas.filter <- seas.filter + seas.filterj
}
seas.filter <- c(rep(0,length(seas.filter)-1),seas.filter)
seas.filter <- seas.filter + rev(seas.filter)
seas.filter <- c(rep(0,(length(seas.filter)-1)/2),1,rep(0,(length(seas.filter)-1)/2)) - seas.filter/(2*p.seas+1)
#plot.ts(seas.filter)
sa.filter <- polymult(detrend.filter,seas.filter)
shift <- (length(sa.filter)-1)/2
sa.filter <- c(1,rep(0,shift)) - rev(sa.filter[1:(shift+1)])
sa.filter <- c(rev(sa.filter),sa.filter[-1])
#plot.ts(sa.filter)

trend.filter <- array(trend.filter,c(1,1,length(trend.filter)))
sa.filter <- array(sa.filter,c(1,1,length(sa.filter)))
trend.comp <- sigex.adhocextract(psi,mdl,data.ts,trend.filter,half.len,0,TRUE)
sa.comp <- sigex.adhocextract(psi,mdl,data.ts,sa.filter,shift,0,TRUE)

 
## get fixed effects
reg.trend <- NULL
reg.trend <- cbind(reg.trend,param[[4]]*rep(1,T))

## plotting
trendcol <- "tomato"
cyccol <- "orchid"
seascol <- "seagreen"
sacol <- "navyblue"
fade <- 60

#pdf(file="  .pdf",height=8,width=10)
plot(data.ts)
sigex.graph(trend.comp,reg.trend,start(data.ts),
            period,1,0,trendcol,fade)
sigex.graph(sa.comp,reg.trend,start(data.ts),
            period,1,0,sacol,fade)
#dev.off()
 

## spectral diagnostics: seasonal adjustment
sigex.specar(sa.comp[[1]],FALSE,1,period)
#dev.off()
