#########################################
###  Script for Weekly BFS Data
#########################################

## wipe
rm(list=ls())

library(devtools)

setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex")
load_all(".")

#####################
### Part I: load data

#load("C:\\Users\\neide\\OneDrive\\Documents\\Research\\Casting\\bfs.RData")




#############################################################
### Part II: Metadata Specifications and Exploratory Analysis

begin <- c(2006,1)
end <- c(2020,27)
period <- 52

## create ts object and plot
dataALL.ts <- sigex.load(bfs,begin,period,c("bfs-ba","bfs-hba","bfs-wba","bfs-cba"),FALSE)

#############################
## select span and transforms

# we focus on "bfs-ba", business applications of bfs
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




###########################################
### PART IV: Model Construction and Fitting

## (a) Initial Model: no holidays

# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"sarma",c(1,1,1,1,52),list(1,1,1,1),"process",1)
mdl <- sigex.meaninit(mdl,data.ts,0)

constraint <- NULL
par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## manage output
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]

## MLE fitting results, no holidays
#  divergence:    -2084.366 lik
#psi.mle <- c(-3.86830013490459, 6.08187233213583, 3.6806305674002, 3.98416618713543,
#             1.65616209571557, 10.9517270937398)
#par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)


## (b) Improved Model: add holidays

# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"sarma",c(1,1,1,1,52),list(1,1,1,1),"process",1)
mdl <- sigex.meaninit(mdl,data.ts,0)

# add regressors
mdl <- sigex.reg(mdl,1,ts(as.matrix(easter.reg),start=start(easter.reg),frequency=period,names="Easter"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(nyd.reg),start=start(nyd.reg),frequency=period,names="NewYearDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(mlk.reg),start=start(mlk.reg),frequency=period,names="MLK"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(gw.reg),start=start(gw.reg),frequency=period,names="GeorgeWashington"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(mem.reg),start=start(mem.reg),frequency=period,names="MemorialDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(ind.reg),start=start(ind.reg),frequency=period,names="IndependenceDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(labor.reg),start=start(labor.reg),frequency=period,names="LaborDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(col.reg),start=start(col.reg),frequency=period,names="ColumbusDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(vet.reg),start=start(vet.reg),frequency=period,names="VeteransDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(tg.reg),start=start(tg.reg),frequency=period,names="Thanksgiving"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(xmas.reg),start=start(xmas.reg),frequency=period,names="Xmas"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(black.reg),start=start(black.reg),frequency=period,names="BlackFriday"))
# Note:  IndependenceDay, VeteransDay, and Christmas are automatically removed

constraint <- NULL
par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## manage output
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]

## MLE fitting results, all holidays
#  divergence:    -2103.4 lik
#psi.mle <- c(-3.895857364407, 6.26245279980034, 3.67958193253905, 4.10627094233825,
#1.78483739260471, 10.9684488747207, -0.00360963369959194, -0.214866996734432,
#-0.250155217787393, 0.0850364492677804, 0.105441772823659, 0.13384379803163,
#-0.0815683136196976, -0.185303694782982, 0.0963338582928288)
#par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)


## (c) Final Model: retain holidays NewYears, MLK, and Labor Day,
#       and AO at time 314.

AO.times <- 314
dataNA.ts <- data.ts
dataNA.ts[AO.times] <- NA

# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"sarma",c(1,1,1,1,52),list(1,1,1,1),"process",1)
mdl <- sigex.meaninit(mdl,dataNA.ts,0)

# add regressors
mdl <- sigex.reg(mdl,1,ts(as.matrix(nyd.reg),start=start(nyd.reg),frequency=period,names="NewYearDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(mlk.reg),start=start(mlk.reg),frequency=period,names="MLK"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(labor.reg),start=start(labor.reg),frequency=period,names="LaborDay"))

constraint <- NULL
par.mle <- sigex.default(mdl,dataNA.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle <- sigex.mlefit(dataNA.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## manage output
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]

## MLE fitting results, three holidays
#  divergence:     -22352.565 lik
#psi.mle <- c(-4.2172882182669, 5.57674429648121, 3.40898364772876, 3.15499052371806,
#             1.0722723838464, 10.9310705588039, -0.405574584987, -0.230573115487113,
#             0.127008620952714)
#par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)



## (d) Model checking

# t statistics for parameters
sigex.tstats(mdl,psi.mle,hess,constraint)
#c(-80.1982272544983, 8.77130698395869, 11.8384768942071, 10.3717137333452,
#5.72147253489335, 76.6054586574405, -4.93403242525173, -2.86138888156494,
#1.42997308493168)

# outlier calculations
data.casts <- sigex.midcast(psi.mle,mdl,dataNA.ts,0)

# residual analysis
resid.mle <- sigex.resid(psi.mle,mdl,dataNA.ts)[[1]]
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

## get fixed effects
reg.trend <- as.matrix(param[[4]][1]*mdl[[4]][[1]][,1])
reg.nyd <- as.matrix(param[[4]][2]*mdl[[4]][[1]][,2])
reg.mlk <- as.matrix(param[[4]][3]*mdl[[4]][[1]][,3])
reg.labor <- as.matrix(param[[4]][4]*mdl[[4]][[1]][,4])
dataLIN.ts <- data.ts - (reg.trend + reg.nyd + reg.mlk + reg.labor)

## define trend and SA weekly filters
week.period <- 365.25/7
half.len <- floor(week.period/2)
p.seas <- 1
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
factor <- 2*p.seas + 1
seas.filter <- seas.filter/factor
seas.filter <- c(rep(0,(length(seas.filter)-1)/2),1,rep(0,(length(seas.filter)-1)/2)) - seas.filter
#plot.ts(seas.filter)
sa.filter <- polymult(detrend.filter,seas.filter)
shift <- (length(sa.filter)-1)/2
sa.filter <- c(1,rep(0,shift)) - rev(sa.filter[1:(shift+1)])
sa.filter <- c(rev(sa.filter),sa.filter[-1])
#plot.ts(sa.filter)
trend.filter <- array(trend.filter,c(1,1,length(trend.filter)))
sa.filter <- array(sa.filter,c(1,1,length(sa.filter)))

## compute extractions
trend.comp <- sigex.adhocextract(psi,mdl,dataNA.ts,trend.filter,half.len,0,TRUE)
sa.comp <- sigex.adhocextract(psi,mdl,dataNA.ts,sa.filter,shift,0,TRUE)
AO.errs <- dataLIN.ts[AO.times] - data.casts[[1]]
sa.comp[[1]][AO.times] <- sa.comp[[1]][AO.times] + AO.errs

## plotting
trendcol <- "tomato"
cyccol <- "orchid"
seascol <- "seagreen"
sacol <- "navyblue"
fade <- 60

#pdf(file="bfs-signal.pdf",height=8,width=10)
plot(data.ts)
#lines(data.ts-as.matrix(reg.nyd+reg.mlk+reg.labor),col=cyccol)
sigex.graph(trend.comp,reg.trend,start(data.ts),
            period,1,0,trendcol,fade)
sigex.graph(sa.comp,reg.trend,start(data.ts),
            period,1,0,sacol,fade)
#dev.off()

write(t(cbind(trend.comp[[1]][,1]+reg.trend,
              sa.comp[[1]][,1]+reg.trend)),
              file="signals.dat",ncol=2)

## spectral diagnostics: seasonal adjustment
sigex.specar(sa.comp[[1]],FALSE,1,period)
#dev.off()

