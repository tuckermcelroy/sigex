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

load("C:\\Users\\neide\\OneDrive\\Documents\\Research\\Casting\\bfs.RData")

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
end <- c(2020,27)
period <- 52

## create ts object and plot
dataALL.ts <- sigex.load(bfs,begin,period,"bfs",FALSE)

#############################
## select span and transforms

transform <- "log"
aggregate <- FALSE
subseries <- 1
range <- NULL
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,FALSE)


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

# code to get calendar date for Sunday of the first week
first.day <- 1
start.year <- begin[1]
start.week <- begin[2]
day.lead <- day2week(c(1,1,start.year)) - first.day
if(day.lead < 0) { day.lead <- day.lead + 7 }
day.index <- 7*(start.week-1) - day.lead + 1
year.index <- start.year
if(day.index <= 0)
{
  day.index <- date2day(12,31,start.year-1) + day.index
  year.index <- year.index-1
}
start.date <- day2date(day.index-1,c(1,1,year.index))
end.date <- day2date(day.index-2 + 7*T,c(1,1,year.index))

##############################
## Generate holiday regressors

easter.dates <- read.table("data\\easter500.txt")
easter.reg <- gethol(easter.dates,7,0,start.date,end.date)

nyd.dates <- read.table("data\\newyear500.txt")
nyd.reg <- gethol(nyd.dates,7,0,start.date,end.date)

mlk.dates <- read.table("data\\mlk500.txt")
mlk.reg <- gethol(mlk.dates,7,0,start.date,end.date)

gw.dates <- read.table("data\\gw500.txt")
gw.reg <- gethol(gw.dates,7,0,start.date,end.date)

mem.dates <- read.table("data\\mem500.txt")
mem.reg <- gethol(mem.dates,7,0,start.date,end.date)

ind.dates <- read.table("data\\ind500.txt")
ind.reg <- gethol(ind.dates,7,0,start.date,end.date)

labor.dates <- read.table("data\\labor500.txt")
labor.reg <- gethol(labor.dates,7,0,start.date,end.date)

col.dates <- read.table("data\\columbus500.txt")
col.reg <- gethol(col.dates,7,0,start.date,end.date)

vet.dates <- read.table("data\\vet500.txt")
vet.reg <- gethol(vet.dates,7,0,start.date,end.date)

tg.dates <- read.table("data\\thanksgiving500.txt")
tg.reg <- gethol(tg.dates,7,0,start.date,end.date)

xmas.dates <- read.table("data\\xmas500.txt")
xmas.reg <- gethol(xmas.dates,7,0,start.date,end.date)

black.dates <- read.table("data\\black400.txt")
black.reg <- gethol(black.dates,7,0,start.date,end.date)

####################################
## Convert to weekly flow regressors

easter.reg <- sigex.daily2weekly(easter.reg,first.day,start.date)
easter.reg <- rowSums(easter.reg)/7

nyd.reg <- sigex.daily2weekly(nyd.reg,first.day,start.date)
nyd.reg <- rowSums(nyd.reg)/7

mlk.reg <- sigex.daily2weekly(mlk.reg,first.day,start.date)
mlk.reg <- rowSums(mlk.reg)/7

gw.reg <- sigex.daily2weekly(gw.reg,first.day,start.date)
gw.reg <- rowSums(gw.reg)/7

mem.reg <- sigex.daily2weekly(mem.reg,first.day,start.date)
mem.reg <- rowSums(mem.reg)/7

ind.reg <- sigex.daily2weekly(ind.reg,first.day,start.date)
ind.reg <- rowSums(ind.reg)/7

labor.reg <- sigex.daily2weekly(labor.reg,first.day,start.date)
labor.reg <- rowSums(labor.reg)/7

col.reg <- sigex.daily2weekly(col.reg,first.day,start.date)
col.reg <- rowSums(col.reg)/7

vet.reg <- sigex.daily2weekly(vet.reg,first.day,start.date)
vet.reg <- rowSums(vet.reg)/7

tg.reg <- sigex.daily2weekly(tg.reg,first.day,start.date)
tg.reg <- rowSums(tg.reg)/7

xmas.reg <- sigex.daily2weekly(xmas.reg,first.day,start.date)
xmas.reg <- rowSums(xmas.reg)/7

black.reg <- sigex.daily2weekly(black.reg,first.day,start.date)
black.reg <- rowSums(black.reg)/7


##############
## Basic Model

# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"sarma",c(1,1,1,1,52),list(1,1,1,1),"process",1)
mdl <- sigex.meaninit(mdl,data.ts,0)

reg <- ts(as.matrix(easter.reg),start=start(easter.reg),
          frequency=period,
          names="Easter")

mdl <- sigex.reg(mdl,1,reg)


mdl <- sigex.reg(mdl,1,ts(as.matrix(easter.reg),
                            start=start(easter.reg),
                            frequency=period,
                            names="Easter"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(nyd.reg),
                          start=start(nyd.reg),
                          frequency=period,
                          names="NewYearDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(mlk.reg),
                          start=start(mlk.reg),
                          frequency=period,
                          names="MLK"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(gw.reg),
                          start=start(gw.reg),
                          frequency=period,
                          names="GeorgeWashington"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(mem.reg),
                          start=start(mem.reg),
                          frequency=period,
                          names="MemorialDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(ind.reg),
                          start=start(ind.reg),
                          frequency=period,
                          names="IndependenceDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(labor.reg),
                          start=start(labor.reg),
                          frequency=period,
                          names="LaborDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(col.reg),
                          start=start(col.reg),
                          frequency=period,
                          names="ColumbusDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(vet.reg),
                          start=start(vet.reg),
                          frequency=period,
                          names="VeteransDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(tg.reg),
                          start=start(tg.reg),
                          frequency=period,
                          names="Thanksgiving"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(xmas.reg),
                          start=start(xmas.reg),
                          frequency=period,
                          names="Xmas"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(black.reg),
                          start=start(black.reg),
                          frequency=period,
                          names="BlackFriday"))


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
#  divergence:    -2076.881 lik
#psi.mle <- c(-3.82130660051201, 6.55294873414615, 3.80965936769506, 3.88542473367833,
#     1.62833607331469, 11.0528714439272)
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
