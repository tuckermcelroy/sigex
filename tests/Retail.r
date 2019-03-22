##########################################
###  Script for Daily Retail Palantir Data
###########################################

## wipe
rm(list=ls())

library(devtools)

setwd("C:\\Users\\Tucker\\Documents\\GitHub\\sigex")
load_all(".")
 

######################
### Part I: load data
 
# automatic


#############################################################
### Part II: Metadata Specifications and Exploratory Analysis

start.date <- c(10,1,2012)
end.date <- day2date(dim(retail)[1]-1,start.date)
period <- 365

# calendar calculations
start.day <- date2day(start.date[1],start.date[2],start.date[3])
end.day <- date2day(end.date[1],end.date[2],end.date[3])
begin <- c(start.date[3],start.day) 
end <- c(end.date[3],end.day)

## create ts object and plot
dataALL.ts <- sigex.load(retail,begin,period,colnames(retail),TRUE)
 
#############################
## select span and transforms 

## all data with no transform
transform <- "none"
aggregate <- FALSE
subseries <- 5
range <- NULL
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)

#######################
## spectral exploratory

## levels
par(mfrow=c(1,1))
for(i in 1:length(subseries))
{
	sigex.specar(data.ts,FALSE,i,7)
}
dev.off()

## growth rates
par(mfrow=c(1,1))
for(i in 1:length(subseries))
{
	sigex.specar(data.ts,TRUE,i,7)
}
dev.off()


###############################
### Part III: Model Declaration

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]


##############
## Basic Model
 
## setup holiday regressors

# Easter Day
easter.reg1 <- gethol(easter.dates,0,0,start.date,end.date)
# pre-Easter
easter.reg2 <- gethol(easter.dates,8,-1,start.date,end.date)
# Cyber Monday
cyber.reg <- gethol(cyber.dates,0,0,start.date,end.date)
# Black Friday
black.reg <- gethol(black.dates,0,0,start.date,end.date)
# Superbowl Sunday
super.reg <- gethol(super.dates,0,0,start.date,end.date)
# Labor Day
labor.reg <- gethol(labor.dates,0,0,start.date,end.date)
# Chinese New Year
cny.reg <- gethol(cny.dates,0,0,start.date,end.date)

## setup regressor for data weakness
weak.reg <- rep(0,T)
weak.reg[187:199] <- 1

## week effect regressors
 
temp <- (day2week(start.date)-1 + seq(1,T)) %% 7
# Sunday effect
tdsun.reg <- rep(0,T)
tdsun.reg[temp==1] <- 1
# Monday effect
tdmon.reg <- rep(0,T)
tdmon.reg[temp==2] <- 1
# Tuesday effect
tdtue.reg <- rep(0,T)
tdtue.reg[temp==3] <- 1
# Wednesday effect
tdwed.reg <- rep(0,T)
tdwed.reg[temp==4] <- 1
# Thursday effect
tdthu.reg <- rep(0,T)
tdthu.reg[temp==5] <- 1
# Friday effect
tdfri.reg <- rep(0,T)
tdfri.reg[temp==6] <- 1
# Saturday effect
tdsat.reg <- rep(0,T)
tdsat.reg[temp==0] <- 1

# normalize
tdsun.reg <- tdsun.reg - 1/7
tdmon.reg <- tdmon.reg - 1/7
tdtue.reg <- tdtue.reg - 1/7
tdwed.reg <- tdwed.reg - 1/7
tdthu.reg <- tdthu.reg - 1/7
tdfri.reg <- tdfri.reg - 1/7
tdsat.reg <- tdsat.reg - 1/7


## direct model construction
#mdl <- NULL
#mdl <- sigex.add(mdl,seq(1,N),"sarma",c(8,0,1,0,period),0,"all",1)		
#mdl <- sigex.add(mdl,seq(1,N),"sarma",c(0,7,1,0,period),0,"all",1)		



## model construction
mdl <- NULL
#mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"trend",c(1,-1))		
mdl <- sigex.add(mdl,seq(1,N),"bw",1,c(.8,1,1/7,3/7),"first weekly seasonal",1)      
mdl <- sigex.add(mdl,seq(1,N),"bw",1,c(.8,1,3/7,5/7),"second weekly seasonal",1)      
mdl <- sigex.add(mdl,seq(1,N),"bw",1,c(.8,1,5/7,1),"third weekly seasonal",1)      
#mdl <- sigex.add(mdl,seq(1,N),"bw",1,c(.5,1,1/period,3/period),"annual cycle",1)
mdl <- sigex.add(mdl,seq(1,N),"sarma",c(0,0,1,0,364),0,"annual cycle",1)
mdl <- sigex.add(mdl,seq(1,N),"arma",c(1,0),0,"transient",1)			
#mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"irregular",1)	


# regressors:
mdl <- sigex.meaninit(mdl,data.ts,0)		
mdl <- sigex.reg(mdl,1,ts(as.matrix(easter.reg1),
	start=begin,frequency=period,names="Easter Day"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(easter.reg2),
	start=begin,frequency=period,names="pre-Easter"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(black.reg),
	start=begin,frequency=period,names="Black Friday"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(labor.reg),
	start=begin,frequency=period,names="Labor Day"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(weak.reg),
	start=begin,frequency=period,names="Weak Span"))

## parameter initialization and checks
par.init <- sigex.default(mdl,data.ts)[[1]]
par.init[[2]][[1]] <- -9
par.init[[2]][[2]] <- -9.3
par.init[[2]][[3]] <- -10.7
par.init[[2]][[4]] <- 0
par.init[[2]][[5]] <- -1
par.init[[3]][[1]] <- c(.999,2/7)
par.init[[3]][[2]] <- c(.998,4/7)
par.init[[3]][[3]] <- c(.998,6/7)
par.init[[3]][[4]] <- c(.5)
par.init[[3]][[5]] <- c(.5)

flag.init <- sigex.default(mdl,data.ts)[[2]]
psi.init <- sigex.par2psi(par.init,flag.init,mdl)

resid.init <- sigex.resid(psi.init,mdl,data.ts)
resid.init <- sigex.load(t(resid.init),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
acf(resid.init,lag.max=2*period)


###########################
### Part IV: MLE Estimation

## setup, and fix parameters as desired
mdl.mle <- mdl
psi.mle <- psi.init
# psi.mle <- rnorm(length(psi.init))	# random initialization
flag.mle <- Im(psi.mle)
par.mle <- sigex.psi2par(psi.mle,mdl.mle,data.ts)

## run fitting
fit.mle <- sigex.mlefit(data.ts,par.mle,flag.mle,mdl.mle,"bfgs")

## manage output
psi.mle[flag.mle==1] <- fit.mle[[1]]$par 
psi.mle <- psi.mle + 1i*flag.mle
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]

##  model checking 
resid.mle <- sigex.resid(psi.mle,mdl.mle,data.ts)
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
sigex.portmanteau(resid.mle,period,length(psi.mle))
sigex.gausscheck(resid.mle)
acf(resid.mle,lag.max=2*period)
      
## check on standard errors and condition numbers
print(eigen(hess)$values)
taus <- log(sigex.conditions(data.ts,psi.mle,mdl.mle))
print(taus)
tstats <- sigex.tstats(mdl.mle,psi.mle,hess)
stderrs <- sigex.psi2par(tstats,mdl,data.ts)
print(tstats)

## bundle  
analysis.mle <- sigex.bundle(data.ts,transform,mdl.mle,psi.mle)

# HERE : notes model needs an annual component due to residual correlation at lag 365

##########################################
### Part V: Signal Extraction based on MLE

## load up the MLE fit for signal extraction
data.ts <- analysis.mle[[1]]
mdl <- analysis.mle[[3]]
psi <- analysis.mle[[4]]
param <- sigex.psi2par(psi,mdl,data.ts)

## get signal filters
signal.trend <- sigex.signal(data.ts,param,mdl,1)
signal.cycle <- sigex.signal(data.ts,param,mdl,2)
signal.irr <- sigex.signal(data.ts,param,mdl,3)
signal.noncyc <- sigex.signal(data.ts,param,mdl,c(1,3))
signal.hp <- sigex.signal(data.ts,param,mdl,c(2,3))
signal.lp <- sigex.signal(data.ts,param,mdl,c(1,2))

## get extractions 
extract.trend <- sigex.extract(data.ts,signal.trend,mdl,param)
extract.cycle <- sigex.extract(data.ts,signal.cycle,mdl,param)
extract.irr <- sigex.extract(data.ts,signal.irr,mdl,param)
extract.noncyc <- sigex.extract(data.ts,signal.noncyc,mdl,param)
extract.hp <- sigex.extract(data.ts,signal.hp,mdl,param)
extract.lp <- sigex.extract(data.ts,signal.lp,mdl,param)
  
## get fixed effects
reg.trend <- NULL
for(i in 1:N) {
reg.trend <- cbind(reg.trend,sigex.fixed(data.ts,mdl,i,param,"Trend")) }

## plotting
trendcol <- "tomato"
cyccol <- "navyblue"
fade <- 40
par(mfrow=c(2,1),mar=c(5,4,4,5)+.1)
for(i in 1:N)
{
	plot(data.ts[,i],xlab="Year",ylab="Trend",ylim=c(min(data.ts[,i])-2,max(data.ts[,i])),
		lwd=2,yaxt="n",xaxt="n")
	sigex.graph(extract.trend,reg.trend,begin.date,period,i,0,trendcol,fade)
	sigex.graph(extract.cycle,NULL,begin.date,period,i,10,cyccol,fade)
	axis(1,cex.axis=1)
	axis(2,cex.axis=1)
	axis(4,cex.axis=1,at=c(9.5,10,10.5,11,11.5),labels=c(9.5,10,10.5,11,11.5)-10)
	mtext("Cycle", side=4, line=3)
}
dev.off()

## transfer function analysis
grid <- 200
frf.trend <- sigex.getfrf(data.ts,param,mdl,1,TRUE,grid)
frf.cycle <- sigex.getfrf(data.ts,param,mdl,2,TRUE,grid)

## filter analysis
len <- 50
wk.trend <- sigex.wk(data.ts,param,mdl,1,TRUE,grid,len)
wk.cycle <- sigex.wk(data.ts,param,mdl,2,TRUE,grid,len)
















##################################################
################################################
################ SCRAP


phi2psi <- function(phi)
{
	p <- length(phi)
	pacfs <- phi[p]
	if(p > 1)
	{
		phi <- as.vector(phi[-p])
		for(j in p:2)
		{
			A.mat <- diag(j-1) - pacfs[1]*diag(j-1)[,(j-1):1,drop=FALSE]
			phi <- solve(A.mat,phi)
			pacfs <- c(phi[j-1],pacfs)
			phi <- phi[-(j-1)]
		}
	}
	psi <- log(1+pacfs) - log(1-pacfs)
	return(psi)
}


## all data with no transform
transform <- "none"
aggregate <- FALSE
subseries <- 5
range <- NULL
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)
trend.reg <- rep(1,T)

X.mat <- cbind(trend.reg,easter.reg1,easter.reg2,black.reg,labor.reg,weak.reg,
	tdmon.reg,tdtue.reg,tdwed.reg,tdthu.reg,tdfri.reg,tdsat.reg)
beta <- solve( t(X.mat) %*% X.mat ) %*% t(X.mat) %*% data.ts
data.ts <- data.ts - X.mat %*% beta 
 
 
plot(data.ts)
acf(data.ts,lag.max=400)
pacf(data.ts,lag.max=400)
x.acf <- acf(data.ts,lag.max=800,type="covariance")$acf


#plot(diff(data.ts,lag=365))
#acf(diff(data.ts,lag=365),lag.max=400)
#x.acf <- acf(diff(data.ts,lag=365),lag.max=800,type="covariance")$acf

#data.diff <- ts(filter(data.ts,c(1,rep(0,363),-3/4,-1/4),method="convolution",sides=1)[366:length(data.ts)])
#plot(data.diff)
#acf(data.diff,lag.max=400)
#x.acf <- acf(data.diff,lag.max=800,type="covariance")$acf

p.order <- 366
phi.ar <- solve(toeplitz(x.acf[1:p.order]),x.acf[2:(p.order+1)])
kappa <- phi2psi(phi.ar)
my.inds <- which(abs(kappa)>.2)
#phi.ar[365] <- .45
#my.inds <- which(abs(phi.ar) > .1)
#phi.ar <- solve(toeplitz(x.acf[my.inds]),x.acf[my.inds+1])
#my.filter <- rep(0,p.order)
#my.filter[my.inds] <- phi.ar
#my.filter <- c(1,-1*my.filter)

#my.filter <- polymult(my.filter,c(1,rep(0,364),-1))

#my.filter <- polymult(c(1,-1*phi.ar),c(1,0,rep(0,362),-.6,-.3))
#my.filter <- polymult(c(1,-1*phi.ar),c(1,rep(0,363),-.3))
rho <- 1
omega <- 365
phi.ar <- c(2*rho*cos(2*pi/omega),-rho^2)
my.filter <- c(1,-1*phi.ar)

ent.ts <- filter(data.ts,my.filter,method="convolution",sides=1)[length(my.filter):length(data.ts)]
plot(ts(ent.ts))
acf(ent.ts,lag.max=p.order)
pacf(ent.ts,lag.max=p.order) 
 
phi.mat <- rbind(phi.ar,diag(p.order)[1:(p.order-1),])
sort(1/Mod(eigen(phi.mat)$values))
spec.ar(data.ts,order=p.order)

ent.acf <- ARMAauto(phi.ar,NULL,800)
plot(ts(ent.acf))


## AR gap model
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"arma",c(p.order,0),0,"process",1)


## alternate model construction
mdl <- NULL
#mdl <- sigex.add(mdl,seq(1,N),"arma",c(2,0),0,"trend",c(1,-1))		
mdl <- sigex.add(mdl,seq(1,N),"sarma",c(0,0,1,0,7),0,"weekly cycle",1)
#mdl <- sigex.add(mdl,seq(1,N),"sarma",c(0,2,1,0,364),0,"annual cycle",1)
mdl <- sigex.add(mdl,seq(1,N),"sarma",c(3,3,1,0,365),0,"annual cycle",1)
mdl <- sigex.add(mdl,seq(1,N),"arma",c(3,0),0,"transient",1)			
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"irregular",1)	


# regressors:
mdl <- sigex.meaninit(mdl,data.ts,0)		
mdl <- sigex.reg(mdl,1,ts(as.matrix(easter.reg1),
	start=begin,frequency=period,names="Easter Day"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(easter.reg2),
	start=begin,frequency=period,names="pre-Easter"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(black.reg),
	start=begin,frequency=period,names="Black Friday"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(labor.reg),
	start=begin,frequency=period,names="Labor Day"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(weak.reg),
	start=begin,frequency=period,names="Weak Span"))
#mdl <- sigex.reg(mdl,1,ts(as.matrix(tdsun.reg),
#	start=begin,frequency=period,names="Sunday Weekly"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(tdmon.reg),
	start=begin,frequency=period,names="Monday Weekly"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(tdtue.reg),
	start=begin,frequency=period,names="Tuesday Weekly"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(tdwed.reg),
	start=begin,frequency=period,names="Wednesday Weekly"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(tdthu.reg),
	start=begin,frequency=period,names="Thursday Weekly"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(tdfri.reg),
	start=begin,frequency=period,names="Friday Weekly"))
mdl <- sigex.reg(mdl,1,ts(as.matrix(tdsat.reg),
	start=begin,frequency=period,names="Saturday Weekly"))



## parameter initialization and checks
par.init <- sigex.default(mdl,data.ts)[[1]]
par.init[[4]] <- beta

#data.mat <- matrix(data.ts[1:1288],nrow=7)
#temp <- rowMeans(data.mat)
#temp <- temp - mean(temp)
#par.init[[2]][[1]] <- -3.25
#par.init[[3]][[1]] <- c(-.5,.6)
#par.init[[4]]  <- c(0.9800019653, -0.9501984918,  0.1729073488,  3.4112734771,
#  0.1906021777, -0.7336086723, -0.0085941787, -0.0388358646, -0.0005719893,
#  0.0537246269,  0.3449043921,  0.7922742353)


flag.init <- sigex.default(mdl,data.ts)[[2]]
psi.init <- sigex.par2psi(par.init,flag.init,mdl)

resid.init <- sigex.resid(psi.init,mdl,data.ts)
resid.init <- sigex.load(t(resid.init),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
acf(resid.init,lag.max=2*period)


## setup, and fix parameters as desired
mdl.mle <- mdl
psi.mle <- psi.init
# psi.mle <- rnorm(length(psi.init))	# random initialization
flag.mle <- Im(psi.mle)
flag.mle[1+setdiff(seq(1,p.order),my.inds)] <- 0
par.mle <- sigex.psi2par(psi.mle,mdl.mle,data.ts)

## run fitting
fit.mle <- sigex.mlefit(data.ts,par.mle,flag.mle,mdl.mle,"bfgs",whittle=TRUE)



## manage output
psi.mle[flag.mle==1] <- fit.mle[[1]]$par 
psi.mle <- psi.mle + 1i*flag.mle
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]

##  model checking 
resid.mle <- sigex.resid(psi.mle,mdl.mle,data.ts)
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
sigex.portmanteau(resid.mle,period,length(psi.mle))
sigex.gausscheck(resid.mle)
acf(resid.mle,lag.max=2*period)


 
####

par.last <- sigex.psi2par(psi.last,mdl.mle,data.ts)
resid.last <- sigex.resid(psi.last,mdl.mle,data.ts)
resid.last <- sigex.load(t(resid.last),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
acf(resid.last,lag.max=2*period)



#####


## wipe
rm(list=ls())

library(devtools)

setwd("C:\\Users\\Tucker\\Documents\\GitHub\\sigex")
load_all(".")
 


start.date <- c(10,1,2012)
end.date <- day2date(dim(retail)[1]-1,start.date)
period <- 365

# calendar calculations
start.day <- date2day(start.date[1],start.date[2],start.date[3])
end.day <- date2day(end.date[1],end.date[2],end.date[3])
begin <- c(start.date[3],start.day) 
end <- c(end.date[3],end.day)

## create ts object and plot
dataALL.ts <- sigex.load(retail,begin,period,colnames(retail),TRUE)
 

## all data with no transform
transform <- "none"
aggregate <- FALSE
subseries <- 5
range <- NULL
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)
trend.reg <- rep(1,T)

X.mat <- cbind(trend.reg,easter.reg1,easter.reg2,black.reg,labor.reg,weak.reg,
	tdmon.reg,tdtue.reg,tdwed.reg,tdthu.reg,tdfri.reg,tdsat.reg)
beta <- solve( t(X.mat) %*% X.mat ) %*% t(X.mat) %*% data.ts
data.ts <- data.ts - X.mat %*% beta 
T <- dim(data.ts)[1]
 

#### seasonal vector form


data.mat <- t(matrix(data.ts[1:1288],nrow=7))
#acf(data.mat,lag.max=60)

## VARMA model


N <- dim(data.mat)[2]
T <- dim(data.mat)[1]

mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"varma",c(1,0),0,"process",1)
mdl <- sigex.meaninit(mdl,data.mat,0)	

par.init <- sigex.default(mdl,data.mat)[[1]]
flag.init <- sigex.default(mdl,data.mat)[[2]]
psi.init <- sigex.par2psi(par.init,flag.init,mdl)

resid.init <- sigex.resid(psi.init,mdl,data.ts)
resid.init <- sigex.load(t(resid.init),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
acf(resid.init,lag.max=2*period)





x.acf <- acf(data.mat[,1],lag.max=100,type="covariance")$acf
p.order <- 53
phi.ar <- solve(toeplitz(x.acf[1:p.order]),x.acf[2:(p.order+1)])
plot(ts(phi.ar))
kappa <- phi2psi(phi.ar)
my.inds <- which(abs(kappa)>.2)
my.inds







ent.ts <- filter(data.ts,my.filter,method="convolution",sides=1)[length(my.filter):length(data.ts)]
plot(ts(ent.ts))
acf(ent.ts,lag.max=p.order)
pacf(ent.ts,lag.max=p.order) 
 



 
 
