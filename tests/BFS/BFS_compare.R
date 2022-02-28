###  Script for Weekly BFS Data

# Preamble ----


library(devtools)
library(Rcpp)

## * set directories ----
# suppose directory is set to where sigex is located, e.g.
# setwd("C:\\Users\\neide\\OneDrive\\Documents\\GitHub\\sigex")
load_all(".")
root.dir <- "~\\GitHub\\sigex"
ex.dir <- file.path(root.dir, "tests/BFS")

# #### Part I: load data ######################################################

# automatic

# #### Part II: Metadata Specifications and Model Selection ##############

begin <- c(2006,1)
end <- c(2020,27)
period <- 52

# * ts object and plot ----
dataALL.ts <- sigex.load(bfs[,3:6],
                         begin,
                         period,
                         c("bfs-ba","bfs-hba","bfs-wba","bfs-cba"),
                         FALSE)

# * select span and transforms ----
# we focus on "bfs-ba", business applications of bfs
transform <- "log"
aggregate <- FALSE
subseries <- 1
range <- NULL
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)




# #### Part III: Model Declaration ############################################

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]

# code to get calendar date for Sunday of the first week
first.day <- 1
all.date <- weekly2date(first.day,begin,T)
start.date <- all.date[[1]]
end.date <- all.date[[2]]


# ---- * holiday regressors ----

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

# ---- * weekly flow regressors ----

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


# #### PART IV: Model Construction ############################################

# ---- * differencing operator ----
delta.s <- ubgenerator(365.25/7,NULL,1000,1)
delta.full <- polymult(c(1,-1), delta.s)

# ---- * Define JD+ differencing operator ----
s.period <- 365.25/7
rho.s <- 1
s.div <- floor(s.period)
s.frac <- s.period - s.div
sar.op <- c(1, rep(0, s.div - 1), (s.frac - 1) * rho.s, -1 * s.frac * rho.s)


# ---- * Model Selection ----
dataJD <- filter(data.ts, filter = sar.op, sides = 1)
dataUS <- filter(data.ts, filter = delta.full, sides = 1)

par(mfrow = c(1, 2))
acf(dataJD, na.action = na.pass, lag.max = 4 * 53)
acf(dataUS, na.action = na.pass, lag.max= 4 * 53)

pacf(dataJD, na.action = na.pass, lag.max = 4 * 53)
pacf(dataUS, na.action = na.pass, lag.max= 4 * 53)




# #### PART V: Model Fitting ##################################################

# ---- * (a) Initial Model: no holidays ----

# ---- * (a, JD+) ----

# model construction
mdl <- NULL
mdl <- sigex.add(mdl    = mdl,
                 vrank  = seq(1,N),
                 class  = "sarmaf",
                 order  = c(2, 1, 0, 1, 365.25/7),
                 bounds = NULL,
                 name   = "process",
                 delta  = sar.op )
mdl <- sigex.meaninit(mdl,data.ts,0)

constraint <- NULL
par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle_a_jd <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

#params
psi_a_jd <-           fit.mle_a_jd[[1]]$par
par_a_jd <- sigex.psi2par(psi_a_jd, mdl, data.ts)


# ---- * (a, US) ----

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


# params
psi_a_us <- fit.mle_a_us[[1]]$par
par_a_us <- sigex.psi2par(psi_a_us, mdl, data.ts)















# ---- * (b) Add holidays model ----

# ---- * (b, JD+) ----
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

#params
psi_b_jd <-           fit.mle_b_jd[[1]]$par
par_b_jd <- sigex.psi2par(psi_b_jd, mdl, data.ts)




# ---- * (b, US) ----
mdl$type[[1]][[1]] <- "sarma"
mdl$diffop[[1]] <- delta.full

constraint <- NULL
par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle_b_us <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

#params
psi_b_us <-           fit.mle_b_us[[1]]$par
par_b_us <- sigex.psi2par(psi_b_us, mdl, data.ts)








# ---- * (c) Refined Model ----
#   retain holidays NewYears, MLK, and Labor Day, and AO at time 314.

# ---- * (c, JD+) ----

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

#params
psi_c_jd <-           fit.mle_c_jd[[1]]$par
par_c_jd <- sigex.psi2par(psi_c_jd, mdl, data.ts)



# ---- * (c, US) ----
mdl$type[[1]][[1]] <- "sarma"
mdl$diffop[[1]] <- delta.full

par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle_c_us <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

#params
psi_c_us <-           fit.mle_c_us[[1]]$par
par_c_us <- sigex.psi2par(psi_c_us, mdl, data.ts)
























# ---- * (d) Visual Model ----
#   look at ACF and PACF of differenced series and pick SMA(1, 1) model
#   Regressors: NewYears, MLK, and Labor Day, and AO at time 314.

# ---- * (d, JD+) ----

# Setup AO as missing value
AO.times <- 314
dataNA.ts <- data.ts
dataNA.ts[AO.times] <- NA

# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"sarmaf",c(0,1,0,1,365.25/7),NULL,"process",sar.op)
mdl <- sigex.meaninit(mdl,dataNA.ts,0)

# add regressors
mdl <- sigex.reg(mdl,1,ts(as.matrix(nyd.reg),start=start(nyd.reg),frequency=period,names="NewYearDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(mlk.reg),start=start(mlk.reg),frequency=period,names="MLK"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(labor.reg),start=start(labor.reg),frequency=period,names="LaborDay"))

constraint <- NULL
par.mle <- sigex.default(mdl,dataNA.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle_d_jd <- sigex.mlefit(dataNA.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE,thresh=20)

#params
psi_d_jd <-           fit.mle_d_jd[[1]]$par
par_d_jd <- sigex.psi2par(psi_d_jd, mdl, data.ts)



# ---- * (d, US) ----
mdl$type[[1]][[1]] <- "sarma"
mdl$diffop[[1]] <- delta.full

par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle_d_us <- sigex.mlefit(dataNA.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

#params
psi_d_us <-           fit.mle_d_us[[1]]$par
par_d_us <- sigex.psi2par(psi_d_us, mdl, data.ts)











# ---- * (e) auto.arima Model ----
#   choose model based on the results of auto.arima()
#   this requires the data to be twice first differenced
#   Regressors: NewYears, MLK, and Labor Day, and AO at time 314.

library(forecast)
auto.arima(dataJD)
auto.arima(dataUS)
# > auto.arima(dataJD)
#
# Series: dataJD
# ARIMA(1,1,2)(0,0,1)[52]
#
# Coefficients:
#   ar1      ma1      ma2     sma1
# -0.0353  -0.9571  -0.0083  -0.4123
# s.e.   0.5686   0.5679   0.5498   0.0327
#
# sigma^2 estimated as 0.02403:  log likelihood=308.79
# AIC=-607.58   AICc=-607.5   BIC=-584.8
#
# > auto.arima(dataUS)
#
# Series: dataUS
# ARIMA(1,1,2)(0,0,1)[52]
#
# Coefficients:
#   ar1     ma1      ma2     sma1
# -0.3039  0.0266  -0.9504  -0.4080
# s.e.   0.0406  0.0146   0.0142   0.0335
#
# sigma^2 estimated as 0.02829:  log likelihood=250.24
# AIC=-490.47   AICc=-490.39   BIC=-467.69


# ---- * (e, JD+) ----

# Setup AO as missing value
AO.times <- 314
dataNA.ts <- data.ts
dataNA.ts[AO.times] <- NA

# New JD+ differencing operator after additional 1st diff
sar.op2 <- polymult(c(1,-1), sar.op)

# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"sarmaf",c(1,2,0,1,365.25/7),NULL,"process",sar.op2)
mdl <- sigex.meaninit(mdl,dataNA.ts,0)

# add regressors
mdl <- sigex.reg(mdl,1,ts(as.matrix(nyd.reg),start=start(nyd.reg),frequency=period,names="NewYearDay"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(mlk.reg),start=start(mlk.reg),frequency=period,names="MLK"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(labor.reg),start=start(labor.reg),frequency=period,names="LaborDay"))

constraint <- NULL
par.mle <- sigex.default(mdl,dataNA.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle_e_jd <- sigex.mlefit(dataNA.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE,thresh=20)


#params
psi_e_jd <-           fit.mle_e_jd[[1]]$par
par_e_jd <- sigex.psi2par(psi_e_jd, mdl, data.ts)


# ---- * (e, US) ----

# New differencing operator after additonal 1st diff
delta.full2 <- polymult(c(1, -1), delta.full)

mdl$type[[1]][[1]] <- "sarma"
mdl$diffop[[1]] <- delta.full2

par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle_e_us <- sigex.mlefit(dataNA.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

#params
psi_e_us <-           fit.mle_e_us[[1]]$par
par_e_us <- sigex.psi2par(psi_e_us, mdl, data.ts)








# ---- Fix Param ----

# here we force SMA param to be -.8
# We want to understand the negative seasonal lag in model residuals

sigex.lik(psi = psi_d_us, mdl = mdl, data.ts = dataNA.ts)
psi_d_us_zeroReg <- psi_d_us
psi_d_us_zeroReg[5:7] <- 0
sigex.lik(psi = psi_d_us_zeroReg, mdl = mdl, data.ts = dataNA.ts)

psi_d_us_zeroReg_fixedTheta <- psi_d_us_zeroReg
psi_d_us_zeroReg_fixedTheta[3] <- -2.1972245773
sigex.lik(psi = psi_d_us_zeroReg_fixedTheta, mdl = mdl, data.ts = dataNA.ts)


# Zero out regression effects
par_temp <- par_d_us
par_temp[[3]][[1]][2] <- -.8
par_temp[[3]][[1]][1] <- -.5
# par_temp[[4]][2:4] <- 0
par_temp
psi_temp <- sigex.par2psi(par_temp, mdl)

sigex.lik(psi = psi_temp, mdl = mdl, data.ts = dataNA.ts)
e <- sigex.resid(psi = psi_temp, mdl = mdl, data.ts = dataNA.ts)
acf(c(e[[1]]), lag.max = 52*4)

# model construction
mdl <- NULL
mdl <- sigex.add(mdl    = mdl,
                 vrank  = seq(1,N),
                 class  = "sarma",
                 order  = c(0,1,0,1,365.25/7),
                 bounds = NULL,
                 name   = "process",
                 delta  = delta.full )
mdl <- sigex.meaninit(mdl,data.ts,0)


constraint <- matrix(c(-2.197, 0, 0, 1, 0), nrow = 1) # Constraint making SMA = -.8
# where did -2.197 come from?


par.mle <- sigex.default(mdl, data.ts, constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting:
fit.mle_f_us <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

psi <- sigex.eta2psi(eta = fit.mle_f_us[[1]]$par,
                     constraint = constraint)
par <- sigex.psi2par(psi, mdl, data.ts)

sigex.lik(psi = psi, mdl = mdl, data.ts = dataNA.ts)









# ---- * Save Image ----
# save.image(file = "model_d_e.RData")



# #### Analysis of Results ####################################################

# ---- * Pick model to look at ----
fit.mle <- fit.mle_f_us

# ---- * manage output ----
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par, constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]

## check on standard errors and get t statistics
eigen(hess)$values

# residual analysis
resid.mle <- sigex.resid(psi.mle, mdl, data.ts)[[1]]
resid.mle <- sigex.load(t(resid.mle), start(data.ts), frequency(data.ts), colnames(data.ts), TRUE)
resid.acf <- acf(resid.mle, lag.max = 4 * 53, plot = TRUE)$acf

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
# dev.off()

# outlier calculations
data.casts <- sigex.midcast(psi.mle,mdl,dataNA.ts,0)


# #### Part VI: Signal Extraction #############################################

# bundle
analysis.mle <- sigex.bundle(dataNA.ts,transform,mdl,psi.mle)

# ---- * load fitted model ----
dataNA.ts <- analysis.mle[[1]]
mdl <- analysis.mle[[3]]
psi <- analysis.mle[[4]]
param <- sigex.psi2par(psi,mdl,dataNA.ts)

# ---- * get fixed effects ----
reg.trend <- sigex.fixed(data.ts,mdl,1,param,"Trend")
reg.nyd <- sigex.fixed(data.ts,mdl,1,param,"NewYearDay")
reg.mlk <- sigex.fixed(data.ts,mdl,1,param,"MLK")
dataLIN.ts <- data.ts - ts(reg.trend + reg.nyd + reg.mlk,
                           start=start(data.ts),frequency=period)

# ---- * trend and SA weekly filters ----
week.period <- 365.25/7
half.len <- floor(week.period/2)
x11.filters <- x11filters(week.period,1)
trend.filter <- x11.filters[[1]]
seas.filter <- x11.filters[[2]]
sa.filter <- x11.filters[[3]]
shift <- (dim(sa.filter)[3]-1)/2

# ---- * extractions -----
trend.comp <- sigex.adhocextract(psi,mdl,dataNA.ts,trend.filter,half.len,0,TRUE)
sa.comp <- sigex.adhocextract(psi,mdl,dataNA.ts,sa.filter,shift,0,TRUE)
AO.errs <- dataLIN.ts[AO.times] - data.casts[[1]]
sa.comp[[1]][AO.times] <- sa.comp[[1]][AO.times] + AO.errs

# ---- * plotting ----
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

# ---- * spectral diagnostics: seasonal adjustment ----
sigex.specar(sa.comp[[1]],FALSE,1,period)
dev.off()



