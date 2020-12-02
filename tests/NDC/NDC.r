########################
#### Script for NDC Data
########################

## wipe
rm(list=ls())

library(devtools)
library(Rcpp)

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
var.out <- var2.par2pre(par.yw)
psi.init <- as.vector(c(covmat.yw[[1]][2,1],log(covmat.yw[[2]]),
                        var.out,colMeans(diff(ts(ndc[2:T,])))))

## model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"varma",c(p.order,0),NULL,"process",c(1,-1))
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
#  divergence:    11290.17
#psi.mle <- c(1.53390034116153, 14.2420315924048, 17.0324361912539, -0.471755451431251,
#0.489727692311396, -0.017677234835485, -0.392499692057969, -0.39909419297684,
#-0.0522002118651599, 0.161306450376942, -0.300580812314135, 0.173188276541793,
#0.366507618868933, 0.181796113107285, -0.0644480542794434, 0.136113149623089,
#0.496532863295004, 0.013062477268366, -0.0818707058555948, 0.0102111694634749,
#0.445693008411153, 0.0898405276107301, -0.268826276410836, 0.0903359864828166,
#0.000320560369556658, 0.0484856870030809, -0.144308229871243,
#-0.0089375652660538, -0.0907709494584039, 0.0516025826426896,
#-0.0326238274507487, 75.9102025269533, 77.8748319106466)

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
tstats <- sigex.tstats(mdl,psi.mle,hess,constraint)
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


