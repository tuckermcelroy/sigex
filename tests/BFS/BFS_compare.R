#########################################
###  Script for Weekly BFS Data
#########################################

## wipe
rm(list=ls())

library(devtools)
library(Rcpp)

# suppose directory is set to where sigex is located, e.g.
# setwd("C:\\Users\\neide\\OneDrive\\Documents\\GitHub\\sigex")
load_all(".")
root.dir <- "~\\GitHub\\sigex"
ex.dir <- file.path(root.dir, "tests/BFS")

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

# Define our differencing operator
delta.s <- ubgenerator(365.25/7,NULL,1000,1)
delta.full <- polymult(c(1,-1),delta.s)

# Define JD+ differencing operator
s.period <- 365.25/7
rho.s <- 1
s.div <- floor(s.period)
s.frac <- s.period - s.div
sar.op <- c(1, rep(0, s.div - 1), (s.frac - 1) * rho.s, -1 * s.frac * rho.s)

# ---- (a) Initial Model: no holidays ----

# ---- (a, JD+) ----

# model construction
mdl <- NULL
mdl <- sigex.add(mdl    = mdl,
                 vrank  = seq(1,N),
                 class  = "sarmaf",
                 order  = c(2,1,0,1,365.25/7),
                 bounds = NULL,
                 name   = "process",
                 delta  = sar.op )
mdl <- sigex.meaninit(mdl,data.ts,0)

constraint <- NULL
par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle_a_jd <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## manage output
# psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
# hess <- fit.mle[[1]]$hessian
# par.mle <- fit.mle[[2]]

# ---- (a, US) ----

# model construction
mdl <- NULL
mdl <- sigex.add(mdl    = mdl,
                 vrank  = seq(1,N),
                 class  = "sarma",
                 order  = c(2,1,0,1,365.25/7),
                 bounds = NULL,
                 name   = "process",
                 delta  = delta.full )
mdl <- sigex.meaninit(mdl,data.ts,0)

constraint <- NULL
par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle_a_us <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## manage output
# psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
# hess <- fit.mle[[1]]$hessian
# par.mle <- fit.mle[[2]]

# ---- (b) Improved Model: add holidays ----

# ---- (b, JD+) ----
# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"sarmaf",c(2,1,0,1,365.25/7),NULL,"process",sar.op)
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
fit.mle_b_jd <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

# ---- (b, US) ----
mdl$type[[1]][[1]] <- "sarma"
mdl$diffop[[1]] <- delta.full

constraint <- NULL
par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle_b_us <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

# ---- (c) Final Model ----
#   retain holidays NewYears, MLK, and Labor Day, and AO at time 314.

# ---- (c, JD+) ----

AO.times <- 314
dataNA.ts <- data.ts
dataNA.ts[AO.times] <- NA

# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"sarmaf",c(2,0,0,1,365.25/7),NULL,"process",sar.op)
mdl <- sigex.meaninit(mdl,dataNA.ts,0)

# add regressors
#mdl <- sigex.reg(mdl,1,ts(as.matrix(nyd.reg),start=start(nyd.reg),frequency=period,names="NewYearDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(mlk.reg),start=start(mlk.reg),frequency=period,names="MLK"))

constraint <- NULL
par.mle <- sigex.default(mdl,dataNA.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle_c_jd <- sigex.mlefit(dataNA.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE,thresh=20)

## manage output
# psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
# hess <- fit.mle[[1]]$hessian
# par.mle <- fit.mle[[2]]

## MLE fitting results, two holidays
#  divergence:     -2329.284
#psi.mle <- c(-4.1429305511094, 5.66956218425899, 3.31919328004112, 3.01604771505878,
#  1.03883169244118, 10.9941148479326, -0.42082779363966, -0.231945823550937)
#par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

# ---- (c, US) ----
mdl$type[[1]][[1]] <- "sarma"
mdl$diffop[[1]] <- delta.full

par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle_c_us <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

# ---- Save Image ----
save.image()



fit.mle_c_jd


