########################
#### Script for NDC Data
########################

## wipe
rm(list=ls())

library(devtools)

# suppose directory is set to where sigex is located, e.g.
#setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex")
load_all(".")
root.dir <- getwd()
setwd(paste(root.dir,"/tests/NDC",sep=""))

######################
### Part I: load data

# automatic

#############################################################
### Part II: Metadata Specifications and Exploratory Analysis

start.date <- c(1992,1)
end.date <- c(2020,5)
period <- 12

## create ts object and plot
dataALL.ts <- sigex.load(ndc,start.date,period,
                         c("Shipments","NewOrders"),TRUE)


#############################
## select span and transforms

## all data with log transform
transform <- "none"
aggregate <- FALSE
subseries <- c(1,2)
begin.date <- start(dataALL.ts)
end.date <- end(dataALL.ts)
range <- NULL
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)


###############################
### Part III: Model Declaration

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]

###################
## Basic Model: VAR

## preliminary analysis
ar.fit <- ar.yw(diff(ts(ndc[2:T,])))
p.order <- ar.fit$order
par.yw <- aperm(ar.fit$ar,c(2,3,1))
covmat.yw <- getGCD(ar.fit$var.pred,2)
var.out <- var.par2pre(par.yw)
psi.init <- as.vector(c(covmat.yw[[1]][2,1],log(covmat.yw[[2]]),
                        var.out,colMeans(diff(ts(ndc[2:T,])))))

## model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"varma",c(p.order,0),
                 list(delta.vec,NULL),"process",c(1,-1))
# regressors:
mdl <- sigex.meaninit(mdl,data.ts,0)


#############################
### Part IV: Model Estimation

## parameter initialization
constraint <- NULL
psi.mle <- psi.init
par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

## run fitting:
fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## MLE fitting results
#  divergence:    11290.23
#psi.mle <- c(1.53454631317219, 14.2420529755789, 17.0326943325739, -1.01484531295215,
#-1.49659230756268, -1.65715595164053, -3.69255666212682, -0.148667723328847,
#-1.68122085524361, -2.4142027430775, 6.85456315483752, 0.818866719387369,
#-2.74235824783242, -2.33303451873513, -1.51075101680441, 3.55513179255684,
#-3.99871784332864, -4.47259883776234, -1.46113313214773, -2.56118230968804,
#-4.8212531368411, -1.50035665350046, -1.66518796626803, -0.849141980158845,
#-4.58898275038447, -4.64767408889132, -23.6742445003945, -0.432191733306375,
#-5.91174154691496, -5.77568558118038, 1.43216029852601, 76.0344509428473,
#77.7796646303267)

## manage output
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]

## residual analysis
resid.mle <- sigex.resid(psi.mle,mdl,data.ts)[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),
                        colnames(data.ts),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*period,plot=TRUE)$acf

## examine condition numbers
log(sigex.conditions(data.ts,psi.mle,mdl))

## model checking
sigex.portmanteau(resid.mle,4*period,length(psi.mle))
sigex.gausscheck(resid.mle)

## check on standard errors and get t statistics
print(eigen(hess)$values)
hess.plus <- hess + diag(length(psi.mle))*max(0,-(1+10e-8)*
                                  min(eigen(hess)$values))
tstats <- sigex.tstats(mdl,psi.mle,hess.plus,constraint)
print(tstats)

## bundle
analysis.mle <- sigex.bundle(data.ts,transform,mdl,psi.mle)


##########################################
### Part V: Casting

## load up the fitted model for casting
data.ts <- analysis.mle[[1]]
mdl <- analysis.mle[[3]]
psi <- analysis.mle[[4]]
param <- sigex.psi2par(psi,mdl,data.ts)

## Generate aftcasts and forecasts with uncertainty
window <- 50
data.casts <- sigex.midcast(psi,mdl,data.ts,window)
extract.casts <- sigex.castextract(data.ts,data.casts,mdl,window,param)

## display
castcol <- "black"
fade <- 60
dataPad.ts <- rbind(matrix(NA,nrow=window,ncol=N),data.ts,matrix(NA,nrow=window,ncol=N))
#pdf(file="NdcCasts.pdf",height=8,width=10)
par(mfrow=c(2,1))
for(i in 1:N)
{
  plot(ts(dataPad.ts[,i],start=start.date,frequency=period),
       xlab="Year",ylab="",lwd=1,col=1)
  sigex.graph(extract.casts,NULL,start.date,period,i,0,castcol,fade)
}
dev.off()


