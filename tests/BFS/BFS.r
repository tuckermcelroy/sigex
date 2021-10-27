#########################################
###  Script for Weekly BFS Data
#########################################

## wipe
rm(list=ls())

library(devtools)
library(Rcpp)

# suppose directory is set to where sigex is located, e.g.
#setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex")
load_all(".")
root.dir <- getwd()
setwd(paste(root.dir,"/tests/BFS",sep=""))


#####################
### Part I: load data

# automatic

#############################################################
### Part II: Metadata Specifications and Exploratory Analysis

begin <- c(2006,1)
end <- c(2020,27)
period <- 52

## create ts object and plot
dataALL.ts <- sigex.load(bfs[,3:6],begin,period,c("bfs-ba","bfs-hba","bfs-wba","bfs-cba"),FALSE)

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
all.date <- weekly2date(first.day,begin,T)
start.date <- all.date[[1]]
end.date <- all.date[[2]]


##############################
## Generate holiday regressors

easter.dates <- read.table(paste(root.dir,"/data/easter500.txt",sep=""))
easter.reg <- gethol(easter.dates,7,0,start.date,end.date)

nyd.dates <- read.table(paste(root.dir,"/data/newyear500.txt",sep=""))
nyd.reg <- gethol(nyd.dates,7,0,start.date,end.date)

mlk.dates <- read.table(paste(root.dir,"/data/mlk500.txt",sep=""))
mlk.reg <- gethol(mlk.dates,7,0,start.date,end.date)

gw.dates <- read.table(paste(root.dir,"/data/gw500.txt",sep=""))
gw.reg <- gethol(gw.dates,7,0,start.date,end.date)

mem.dates <- read.table(paste(root.dir,"/data/mem500.txt",sep=""))
mem.reg <- gethol(mem.dates,7,0,start.date,end.date)

ind.dates <- read.table(paste(root.dir,"/data/ind500.txt",sep=""))
ind.reg <- gethol(ind.dates,7,0,start.date,end.date)

labor.dates <- read.table(paste(root.dir,"/data/labor500.txt",sep=""))
labor.reg <- gethol(labor.dates,7,0,start.date,end.date)

col.dates <- read.table(paste(root.dir,"/data/columbus500.txt",sep=""))
col.reg <- gethol(col.dates,7,0,start.date,end.date)

vet.dates <- read.table(paste(root.dir,"/data/vet500.txt",sep=""))
vet.reg <- gethol(vet.dates,7,0,start.date,end.date)

tg.dates <- read.table(paste(root.dir,"/data/thanksgiving500.txt",sep=""))
tg.reg <- gethol(tg.dates,7,0,start.date,end.date)

xmas.dates <- read.table(paste(root.dir,"/data/xmas500.txt",sep=""))
xmas.reg <- gethol(xmas.dates,7,0,start.date,end.date)

black.dates <- read.table(paste(root.dir,"/data/black400.txt",sep=""))
black.reg <- gethol(black.dates,7,0,start.date,end.date)

## Independence Day, Veteran's Day, and Christmas are purely seasonal
sum(ind.reg^2)
sum(vet.reg^2)
sum(xmas.reg^2)

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

delta.s <- ubgenerator(365.25/7,NULL,1000,1)
delta.full <- polymult(c(1,-1),delta.s)

## (a) Initial Model: no holidays

# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"sarma",c(2,1,0,1,365.25/7),NULL,"process",delta.full)
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
#  divergence:    -2076.881
#psi.mle <- c(-3.82130660051201, 6.55294873414615, 3.80965936769506, 3.88542473367833,
#   1.62833607331469, 11.0528714439272)
#par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)


## (b) Improved Model: add holidays

# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"sarma",c(2,1,0,1,52),NULL,"process",delta.full)
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
#  divergence:    -2089.925
#psi.mle <- c(-3.85211236205511, 10.5732757025848, 3.89562694233651, 4.08558530244727,
#1.80869692181588, 13.1014026340267, -0.000293736512980122, -0.222659717911844,
#-0.250177133568453, 0.0988001854155316, 0.0859847200965163, 0.139098491514792,
#-0.0819199138033471, -0.166531078048871, 0.0762222534500927)
#par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

## check on standard errors and get t statistics
print(eigen(hess)$values)

# residual analysis
resid.mle <- sigex.resid(psi.mle,mdl,data.ts)[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*53,plot=TRUE)$acf

## (c) Final Model: retain holidays NewYears, MLK, and Labor Day,
#       and AO at time 314.

AO.times <- 314
dataNA.ts <- data.ts
dataNA.ts[AO.times] <- NA

# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"sarma",c(3,1,0,1,52),NULL,"process",delta.full)
mdl <- sigex.meaninit(mdl,dataNA.ts,0)

# add regressors
#mdl <- sigex.reg(mdl,1,ts(as.matrix(nyd.reg),start=start(nyd.reg),frequency=period,names="NewYearDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(mlk.reg),start=start(mlk.reg),frequency=period,names="MLK"))

constraint <- NULL
par.mle <- sigex.default(mdl,dataNA.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle <- sigex.mlefit(dataNA.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## manage output
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]

## MLE fitting results, two holidays
#  divergence:     -2329.284
#psi.mle <- c(-4.1429305511094, 5.66956218425899, 3.31919328004112, 3.01604771505878,
#  1.03883169244118, 10.9941148479326, -0.42082779363966, -0.231945823550937)
#par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)



## (d) Model checking

# t statistics for parameters
sigex.tstats(mdl,psi.mle,hess,constraint)

# residual analysis
resid.mle <- sigex.resid(psi.mle,mdl,dataNA.ts)[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*53,plot=FALSE)$acf

## model checking
sigex.portmanteau(resid.mle,4*period,length(psi.mle))
sigex.gausscheck(resid.mle)

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
dev.off()

# outlier calculations
data.casts <- sigex.midcast(psi.mle,mdl,dataNA.ts,0)

# bundle
analysis.mle <- sigex.bundle(dataNA.ts,transform,mdl,psi.mle)


##########################################
### Part V: Signal Extraction

## load up the fitted model for signal extraction
dataNA.ts <- analysis.mle[[1]]
mdl <- analysis.mle[[3]]
psi <- analysis.mle[[4]]
param <- sigex.psi2par(psi,mdl,dataNA.ts)

## get fixed effects
reg.trend <- sigex.fixed(data.ts,mdl,1,param,"Trend")
reg.nyd <- sigex.fixed(data.ts,mdl,1,param,"NewYearDay")
reg.mlk <- sigex.fixed(data.ts,mdl,1,param,"MLK")
dataLIN.ts <- data.ts - ts(reg.trend + reg.nyd + reg.mlk,
                           start=start(data.ts),frequency=period)

## define trend and SA weekly filters
week.period <- 365.25/7
half.len <- floor(week.period/2)
x11.filters <- x11filters(week.period,1)
trend.filter <- x11.filters[[1]]
seas.filter <- x11.filters[[2]]
sa.filter <- x11.filters[[3]]
shift <- (dim(sa.filter)[3]-1)/2

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

#pdf(file="bfs-signal-ba.pdf",height=8,width=10)
plot(data.ts)
sigex.graph(trend.comp,reg.trend,start(data.ts),
            period,1,0,trendcol,fade)
sigex.graph(sa.comp,reg.trend,start(data.ts),
            period,1,0,sacol,fade)
dev.off()

## spectral diagnostics: seasonal adjustment
sigex.specar(sa.comp[[1]],FALSE,1,period)
dev.off()


