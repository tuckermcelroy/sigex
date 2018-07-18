library(mvtnorm)
library(normalQuadProb)

# ---- Simulate Data ----------------------------------------------------------
N = 3
T = 500
t = 1:T
Phi = diag(N)
sig1 = toeplitz(.9^(0:(N-1)))
s1 = gen_trendComp(n = T, Ndim = N, Sig = sig1)
s2 = gen_monthComp(n = T, Ndim = N, Sig = diag(N), Phi = sig1)
s0 = rmvnorm(n = T, mean = rep(0,N), sigma = diag(N))

data = s1+s2+s0
plot(ts(s1))

spec.ar(s1[,3])
abline(v=1/24)

# ---- Modeling ---------------------------------------------------------------
transform = "none"
agg <- FALSE	# set TRUE to aggregate
series <- 1:3
range <- seq(1,dim(data)[1])
x <- t(data)
N <- dim(x)[1]
T <- dim(x)[2]

# ---- Default Model ----------------------------------------------------------

# stochastic effects
delta.trend <- c(1,-1)
def <- c(0,1,0,1)

mdl <- NULL
# trend-cycle (annual cycle)
mdl <- sigex.add(mdl,seq(1,N),"wn",delta.trend,def)
# first atomic monthly seas
mdl <- sigex.add(mdl,seq(1,N),"wn",c(1,-2*cos(2*pi/12),1),def)
# second atomic monthly seas
#mdl <- sigex.add(mdl,seq(1,N),"wn",c(1,-2*cos(4*pi/12),1),def)
# third atomic monthly seas
#mdl <- sigex.add(mdl,seq(1,N),"wn",c(1,-2*cos(6*pi/12),1),def)
# irregular
mdl <- sigex.add(mdl,seq(1,N),"wn",1,def)
# fixed effects [mean effect]
mdl <- sigex.meaninit(mdl,data,0)

par.default <- sigex.default(mdl,data)[[1]]
flag.default <- sigex.default(mdl,data)[[2]]
psi.default <- sigex.par2psi(par.default,flag.default,mdl)

# ---- MOM estimation and reduced specification -------------------------------
mdl.mom <- mdl
par.mom <- sigex.momfit(data,par.default,mdl.mom)
psi.mom <- sigex.par2psi(par.mom,flag.default,mdl.mom)
resid.mom <- sigex.resid(psi.mom,mdl.mom,data)

#thresh <- -6.22
#thresh <- -3.92
thresh <- -1.66

if(N > 1) {
  reduced.mom <- sigex.reduce(data,par.mom,flag.default,mdl.mom,thresh,FALSE)
  mdl.mom <- reduced.mom[[1]]
  par.mom <- reduced.mom[[2]]
  flag.mom <- sigex.default(mdl.mom,data)[[2]]
  psi.mom <- sigex.par2psi(par.mom,flag.mom,mdl.mom)
  #resid.mom <- sigex.resid(psi.mom,mdl.mom,data)
}

log(sigex.conditions(data,psi.mom,mdl.mom))


# model checking
sigex.portmanteau(t(resid.mom),48,length(psi.mom))
sigex.gausscheck(t(resid.mom))

# bundle for default span
analysis.mom <- sigex.bundle(data,transform,mdl.mom,psi.mom)

## Rough: reduced MOM model
data <- analysis.mom[[1]]
mdl <- analysis.mom[[3]]
psi <- analysis.mom[[4]]
param <- sigex.psi2par(psi,mdl,data)

# ---- METHOD 1: DIRECT MATRIX APPROACH ---------------------------------------
signal.trendann <- sigex.signal(data,param,mdl,1)
signal.seas.week1 <- sigex.signal(data,param,mdl,2)
signal.seas.week2 <- sigex.signal(data,param,mdl,3)
signal.seas.week3 <- sigex.signal(data,param,mdl,4)
signal.seas.week <- sigex.signal(data,param,mdl,c(2,3,4))
signal.sa <- sigex.signal(data,param,mdl,c(1,5))

extract.trendann <- sigex.extract(data,signal.trendann,mdl,param)
extract.seas.week1 <- sigex.extract(data,signal.seas.week1,mdl,param)
extract.seas.week2 <- sigex.extract(data,signal.seas.week2,mdl,param)
extract.seas.week3 <- sigex.extract(data,signal.seas.week3,mdl,param)
extract.seas.week <- sigex.extract(data,signal.seas.week,mdl,param)
extract.sa <- sigex.extract(data,signal.sa,mdl,param)

# ---- PLOTS ------------------------------------------------------------------
trendcol <- "tomato"
cyccol <- "orchid"
seascol <- "seagreen"
sacol <- "navyblue"
fade <- 60

subseries <- 1

begin = 1
period = 12

plot(data[,subseries],xlab="Year",ylab="",lwd=1, type="l")
sigex.graph(extract.sa, NULL, begin, period, subseries, 0, sacol, fade)
#sigex.graph(extract.trendirreg, NULL,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trendann, reg.trend,begin,period,subseries,0,trendcol,10)
sigex.graph(extract.seas.ann,NULL,begin,period,subseries,5,seascol,10)
sigex.graph(extract.seas.week,NULL,begin,period,subseries,3,cyccol,fade)

plot(data[,subseries],xlab="Year",ylab="", type="l",lwd=1)
sigex.graph(extract.seas.week1,NULL,begin,period,subseries,5,cyccol,fade)
sigex.graph(extract.seas.week2,NULL,begin,period,subseries,3,cyccol,fade)
sigex.graph(extract.seas.week3,NULL,begin,period,subseries,1,cyccol,fade)
