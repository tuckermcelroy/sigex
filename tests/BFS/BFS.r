#########################################
###  Script for Weekly BFS Data
#########################################

## wipe
rm(list=ls())

library(devtools)

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
bfs.dates <- bfs[,1:2]

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
#  divergence:    -2076.881
#psi.mle <- c(-3.82130660051201, 6.55294873414615, 3.80965936769506, 3.88542473367833,
#1.62833607331469, 11.0528714439272)
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
#  divergence:    -2089.925
#psi.mle <- c(-3.85211236205511, 10.5732757025848, 3.89562694233651, 4.08558530244727,
#1.80869692181588, 13.1014026340267, -0.000293736512980122, -0.222659717911844,
#-0.250177133568453, 0.0988001854155316, 0.0859847200965163, 0.139098491514792,
#-0.0819199138033471, -0.166531078048871, 0.0762222534500927)
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
#  divergence:     -2331.178
#psi.mle <- c(-4.14582235347973, 5.68312818296934, 3.32001033852795, 3.03830040237583,
#1.05421882755007, 10.9953892306492, -0.419201207247771, -0.232180383605472,
#0.1267872310395)
#par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)



## (d) Model checking

# t statistics for parameters
sigex.tstats(mdl,psi.mle,hess,constraint)
#c(-79.2742683431225, 6.2483970476832, 7.11364294700982, 9.72902490405569,
#5.47281717698996, 63.7788853761768, -4.91944353936497, -2.78142108144931,
#1.37750462582758)

# outlier calculations
data.casts <- sigex.midcast(psi.mle,mdl,dataNA.ts,0)

# residual analysis
resid.mle <- sigex.resid(psi.mle,mdl,dataNA.ts)[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
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

setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\Casting\\BFS")

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

#pdf(file="bfs-signal-ba.pdf",height=8,width=10)
plot(data.ts)
#lines(data.ts-as.matrix(reg.nyd+reg.mlk+reg.labor),col=cyccol)
sigex.graph(trend.comp,reg.trend,start(data.ts),
            period,1,0,trendcol,fade)
sigex.graph(sa.comp,reg.trend,start(data.ts),
            period,1,0,sacol,fade)
#dev.off()

#write(t(cbind(trend.comp[[1]][,1]+reg.trend,
#              sa.comp[[1]][,1]+reg.trend)),
#              file="signals-ba.dat",ncol=2)

## spectral diagnostics: seasonal adjustment
sigex.specar(sa.comp[[1]],FALSE,1,period)
#dev.off()



## Use this???
#leads <- c(-rev(seq(0,window-1)),seq(1,T),seq(T+1,T+window))
#data.casts <- sigex.cast(psi,mdl,data.ts,leads)


##  SCRAP

##  structural modeling

week.period <- 365.25/7
half.len <- floor(week.period/2)
delta.seas <- ubgenerator(week.period,half.len-1,1000)

mdl2 <- NULL
mdl2 <- sigex.add(mdl2,seq(1,N),"arma",c(2,0),list(c(1,1),NULL),"trend",c(1,-1))
mdl2 <- sigex.add(mdl2,seq(1,N),"sarma",c(0,2,0,1,52),list(NULL,c(1,1),NULL,1),
                  "seasonal",delta.seas)
mdl2 <- sigex.add(mdl2,seq(1,N),"arma",c(0,0),0,"irregular",1)
mdl2 <- sigex.meaninit(mdl2,data.ts,0)

constraint <- NULL
par.mle <- sigex.default(mdl2,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl2)

## run fitting:
fit.mle2 <- sigex.mlefit(data.ts,par.mle,constraint,mdl2,"bfgs",debug=TRUE)

