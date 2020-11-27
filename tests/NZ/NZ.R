#########################################
###  Script for Daily Immigration Data
#########################################

## wipe
rm(list=ls())

library(devtools)
library(Rcpp)

# suppose directory is set to where sigex is located, e.g.
setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex")
#setwd("/home/tucker/Documents/GitHub/sigex")
root.dir <- getwd()
#setwd(paste(root.dir,"/tests",sep=""))
setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\SigExNew")
sourceCpp('autoVARMA.cpp')
setwd(root.dir)
load_all(".")
setwd(paste(root.dir,"/tests/NZ",sep=""))


######################
### Part I: load data

# automatic: raw data

# processing

n.months <- dim(imm)[1]/32
imm <- imm[-seq(1,n.months)*32,]	# strip out every 32nd row (totals)
imm <- matrix(imm[imm != 0],ncol=6) # strip out 31st False days

# enter regressors
NZregs <- read.table(paste(root.dir,"/data/NZregressors.txt",sep=""))


#############################################################
### Part II: Metadata Specifications and Exploratory Analysis

start.date <- c(9,1,1997)
end.date <- day2date(dim(imm)[1]-1,start.date)
#end.date <- c(7,31,2012)
period <- 365

# calendar calculations
start.day <- date2day(start.date[1],start.date[2],start.date[3])
end.day <- date2day(end.date[1],end.date[2],end.date[3])
begin <- c(start.date[3],start.day)
end <- c(end.date[3],end.day)

## create ts object and plot
dataALL.ts <- sigex.load(imm,begin,period,c("NZArr","NZDep","VisArr","VisDep","PLTArr","PLTDep"),TRUE)


#############################
## select span and transforms

## first series with log transform
transform <- "log"
aggregate <- FALSE
subseries <- 1
range <- NULL
dataONE.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)


#######################
## spectral exploratory

## levels
par(mfrow=c(1,1))
for(i in 1:length(subseries))
{
  sigex.specar(dataONE.ts,FALSE,i,7)
}
dev.off()

## growth rates
par(mfrow=c(1,1))
for(i in 1:length(subseries))
{
  sigex.specar(dataONE.ts,TRUE,i,7)
}
dev.off()


###########################
## embed as a weekly series

first.day <- 1
data.ts <- sigex.daily2weekly(dataONE.ts,first.day,start.date)
plot(data.ts)



###############################
### Part III: Model Declaration

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]

##########################
## Load holiday regressors
##
## NOTES: easter is based on Easter day and day before Easter
##        school1 is beginning of first school holiday,
##          with window for day of and day after.
##        School2 and school3 are analogous for 2nd and 3rd holidays
##        school1e is end of first school holiday,
##          with window for day of and day before.
##        School2e and school3e are analogous for 2nd and 3rd holidays

easter.reg <- NZregs[,1]
school1.reg <- NZregs[,2]
school1e.reg <- NZregs[,3]
school2.reg <- NZregs[,4]
school2e.reg <- NZregs[,5]
school3.reg <- NZregs[,6]
school3e.reg <- NZregs[,7]

###########################
## Embed holiday regressors

easter.reg <- sigex.daily2weekly(easter.reg,first.day,start.date)
school1.reg <- sigex.daily2weekly(school1.reg,first.day,start.date)
school1e.reg <- sigex.daily2weekly(school1e.reg,first.day,start.date)
school2.reg <- sigex.daily2weekly(school2.reg,first.day,start.date)
school2e.reg <- sigex.daily2weekly(school2e.reg,first.day,start.date)
school3.reg <- sigex.daily2weekly(school3.reg,first.day,start.date)
school3e.reg <- sigex.daily2weekly(school3e.reg,first.day,start.date)

# replace ragged NA with zero
easter.reg[is.na(easter.reg)] <- 0
school1.reg[is.na(school1.reg)] <- 0
school1e.reg[is.na(school1e.reg)] <- 0
school2.reg[is.na(school2.reg)] <- 0
school2e.reg[is.na(school2e.reg)] <- 0
school3.reg[is.na(school3.reg)] <- 0
school3e.reg[is.na(school3e.reg)] <- 0


##############
## Basic Model

# SVARMA
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"svarma",c(1,3,1,0,52),NULL,"process",c(1,-1))
mdl <- sigex.meaninit(mdl,data.ts,0)


# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"svarma",c(1,0,1,0,52),NULL,"process",c(1,-1))
mdl <- sigex.meaninit(mdl,data.ts,0)

for(i in 1:N) {
  mdl <- sigex.reg(mdl,i,ts(as.matrix(easter.reg[,i]),
                            start=start(easter.reg),
                            frequency=frequency(easter.reg),
                            names="Easter-day"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school1.reg[,i]),
                            start=start(school1.reg),
                            frequency=frequency(school1.reg),
                            names="School1-Start"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school1e.reg[,i]),
                            start=start(school1e.reg),
                            frequency=frequency(school1e.reg),
                            names="School1-End"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school2.reg[,i]),
                            start=start(school2.reg),
                            frequency=frequency(school2.reg),
                            names="School2-Start"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school2e.reg[,i]),
                            start=start(school2e.reg),
                            frequency=frequency(school2e.reg),
                            names="School2-End"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school3.reg[,i]),
                            start=start(school3.reg),
                            frequency=frequency(school3.reg),
                            names="School3-Start"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school3e.reg[,i]),
                            start=start(school3e.reg),
                            frequency=frequency(school3e.reg),
                            names="School3-End"))
}


##################################
### PART IV: Model Fitting

constraint <- NULL

# model constraints
par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)
psi.len <- length(psi.mle)
constraint <- diag(psi.len)
index.mat <- matrix(seq(1,N^2),nrow=N,ncol=N)
index.constraint <- setdiff(matrix(index.mat,ncol=1),diag(index.mat))
offset <- N*(N+1)/2
index.constraints <- offset + index.constraint
offset <- offset + N^2
index.constraints <- c(index.constraints,offset + index.constraint)
offset <- offset + N^2
index.constraints <- c(index.constraints,offset + index.constraint)
offset <- offset + N^2
index.constraints <- c(index.constraints,offset + index.constraint)
offset <- offset + N^2
index.constraints <- c(index.constraints,offset + index.constraint)
constraint <- constraint[index.constraints,,drop=FALSE]
constraint <- cbind(rep(0,length(index.constraints)),constraint)

# regression constraints
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(2,2,2,2,2,2,2),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(3,3,3,3,3,3,3),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(4,4,4,4,4,4,4),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(5,5,5,5,5,5,5),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(6,6,6,6,6,6,6),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(7,7,7,7,7,7,7),NULL))

par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

#psi.mle <- psi.init
#par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

## run fitting: commented out, this took 2 weeks!
fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## manage output
#psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
#hess <- fit.mle[[1]]$hessian
#par.mle <- fit.mle[[2]]

## MLE fitting results
# divergence -16986.08
psi.mle <- c(0.47517469995357, 0.36556846772533, 0.33254203944139, 0.3082308605032,
             0.373055224939302, 0.417588352519522, 0.713598330367237, 0.521464162565376,
             0.349334254457817, 0.264266755532234, 0.238286390817609, 0.509671491935296,
             0.45890702860529, 0.409346821825133, 0.400588635461024, 0.75659352001764,
             0.675179639292339, 0.555761562607551, 0.566544368439753, 0.47305509420043,
             0.658380006327245, -4.15761912020721, -3.89862500467477, -4.16038016620892,
             -4.27712507803133, -4.15898092216139, -4.25444979859734, -4.02872571339124,
             -0.350782673917829, 0.00986758523639334, -0.175329578308252,
             -0.19183263485956, -0.263055965952826, -0.318524581972687, -0.337054236431742,
             0.0523827460265609, -0.564552820022245, -0.134757517222772, -0.130205089798105,
             -0.0664722359231971, -0.0994146125779702, -0.10345126913527,
             0.181827588149573, 0.12745586929405, -0.416423536825781, -0.0134510372029961,
             -0.103134191756003, -0.119517796908138, -0.0383620173387325,
             0.0862878359007305, 0.0787774257545998, 0.0879625683893564, -0.513728816607305,
             -0.104293875416977, -0.166596991460797, -0.0887883606300169,
             0.0205321116689622, 0.125030906196229, 0.191953122890016, 0.190161546936267,
             -0.422574457218557, -0.041974354617824, -0.0211591620781137,
             0.254348963234842, 0.0448957722531037, -0.0677036441234219, -0.0730853798256886,
             -0.0484342221725277, -0.693919455253131, -0.234235719372947,
             0.45000830677086, 0.391073260345561, 0.339575953559084, 0.320501513429975,
             0.312519726025935, 0.399172504364575, -0.15487468456843, 0.809332733670788,
             0.577418609186509, 0.297597104400988, 0.221699973198504, 0.135266134155234,
             0.0617460848429357, 0.0904708295504799, 0.109699163396578, 0.391892409463993,
             0.32711972062965, 0.0841213649460174, 0.0689840450445127, -0.0171325599345873,
             -0.0388640919586531, -0.0663904938247799, 0.0411362274498638,
             0.216622632247402, 0.243358318737632, 0.116153922772715, 0.100850380759862,
             0.137624417482269, 0.0220820456667854, -0.0203464331056514, 0.0684624214281305,
             0.152186268755705, 0.277999721019927, 0.281410952696035, 0.0956062132761577,
             -0.0630103328965702, -0.115472368150293, -0.121869738055998,
             0.0143810471530451, 0.222345622359049, 0.255152472454988, 0.179935938537309,
             0.0331107633936892, -0.150142586075399, -0.116184455538622, -0.0518205826338443,
             -0.0889433534224612, -0.0162212323304854, 0.26109675735318, -0.0838255641679354,
             0.110394779393117, 0.151972218541225, 0.134578527090108, 0.129584622312092,
             0.217545007487204, 0.209105087870079, 0.00108629064422719, -0.254999295619355,
             -0.127256204283817, 0.0831758105426799, -0.00358150947718402,
             0.034686365553574, 0.0055439782280033, 0.0470067193051992, 0.00189076791099345,
             -0.254999295619355, -0.127256204283817, 0.0831758105426799, -0.00358150947718402,
             0.034686365553574, 0.00554397822800332, 0.166613285358371, 0.00130962676899953,
             -0.254999295619355, -0.127256204283817, 0.08317581054268, -0.00358150947718403,
             0.034686365553574, 0.00554397822800326, -0.181664825297958, 0.000692005459341965,
             -0.254999295619355, -0.127256204283817, 0.0831758105426799, -0.00358150947718403,
             0.034686365553574, 0.0055439782280033, -0.0935055531659118, -9.09882234482756e-05,
             -0.254999295619355, -0.127256204283817, 0.0831758105426799, -0.003581509477184,
             0.034686365553574, 0.00554397822800331, 0.0652111039775332, -0.000470436916659368,
             -0.254999295619355, -0.127256204283817, 0.0831758105426798, -0.00358150947718402,
             0.034686365553574, 0.0055439782280033, -0.0514675037891409, -0.000260167371337442,
             -0.254999295619355, -0.127256204283817, 0.0831758105426798, -0.00358150947718402,
             0.034686365553574, 0.0055439782280033, 0.0414210579401776)
par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

##  model checking
resid.mle <- sigex.resid(psi.mle,mdl,data.ts)[[1]]
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
dev.off()


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



# HERE









#################################
########################## METHOD 1: DIRECT MATRIX APPROACH ############
############################### SKIP ###############################

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

####################################################################
############################# METHOD 2: FORECASTING and WK SIGEX
############## RECOMMENDED METHOD #####################

#grid <- 70000	# high accuracy, close to method 1
#grid <- 700	# low accuracy, but pretty fast
grid <- 7000	# need grid > filter length
window <- 200
horizon <- 2000
#leads <- c(-rev(seq(0,window-1)),seq(1,T),seq(T+1,T+window))
#data.ext <- t(sigex.cast(psi,mdl,data,leads,TRUE))
target <- array(diag(N),c(N,N,1))

extract.trendann <- sigex.wkextract2(psi,mdl,data,1,target,grid,window,horizon,FALSE)
extract.seas.week <- sigex.wkextract2(psi,mdl,data,c(2,3,4),target,grid,window,horizon,FALSE)
extract.seas.week1 <- sigex.wkextract2(psi,mdl,data,2,target,grid,window,horizon,NULL,FALSE)
extract.seas.week2 <- sigex.wkextract2(psi,mdl,data,3,target,grid,window,horizon,NULL,FALSE)
extract.seas.week3 <- sigex.wkextract2(psi,mdl,data,4,target,grid,window,horizon,NULL,FALSE)
extract.sa <- sigex.wkextract2(psi,mdl,data,c(1,5),target,grid,window,horizon,NULL,FALSE)
extract.irr <- sigex.wkextract2(psi,mdl,data,5,target,grid,window,horizon,NULL,FALSE)


##################################################################
#################### LP splitting of trend and cycle #################

cutoff <- pi/365
trunc <- 50000	# appropriate for mu = pi/(365)

extract.trend <- sigex.lpfiltering(mdl,data,1,NULL,psi,cutoff,grid,window,trunc,TRUE)
extract.seas.ann <- sigex.lpfiltering(mdl,data,1,NULL,psi,cutoff,grid,window,trunc,FALSE)
extract.trendirreg <- sigex.lpfiltering(mdl,data,1,5,psi,cutoff,grid,window,trunc,TRUE)



#########################################
### get fixed effects

reg.trend <- NULL


################################# PLOTS

## time series Plots

trendcol <- "tomato"
cyccol <- "orchid"
seascol <- "seagreen"
sacol <- "navyblue"
fade <- 60

subseries <- 1

plot(data[,subseries],xlab="Year",ylab="",ylim=c(2,9),lwd=1)
#sigex.graph(extract.sa,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trendirreg,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trend,reg.trend,begin,period,subseries,0,trendcol,10)
sigex.graph(extract.seas.ann,NULL,begin,period,subseries,5,seascol,10)
sigex.graph(extract.seas.week,NULL,begin,period,subseries,3,cyccol,fade)

plot(data[,subseries],xlab="Year",ylab="",ylim=c(0,9),lwd=1)
sigex.graph(extract.seas.week1,NULL,begin,period,subseries,5,cyccol,fade)
sigex.graph(extract.seas.week2,NULL,begin,period,subseries,3,cyccol,fade)
sigex.graph(extract.seas.week3,NULL,begin,period,subseries,1,cyccol,fade)


subseries <- 2

plot(data[,subseries],xlab="Year",ylab="",ylim=c(2,9),lwd=1)
#sigex.graph(extract.sa,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trendirreg,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trend,reg.trend,begin,period,subseries,0,trendcol,10)
sigex.graph(extract.seas.ann,NULL,begin,period,subseries,5,seascol,10)
sigex.graph(extract.seas.week,NULL,begin,period,subseries,3,cyccol,fade)

plot(data[,subseries],xlab="Year",ylab="",ylim=c(0,9),lwd=1)
sigex.graph(extract.seas.week1,NULL,begin,period,subseries,5,cyccol,fade)
sigex.graph(extract.seas.week2,NULL,begin,period,subseries,3,cyccol,fade)
sigex.graph(extract.seas.week3,NULL,begin,period,subseries,1,cyccol,fade)


subseries <- 3

plot(data[,subseries],xlab="Year",ylab="",ylim=c(2,9),lwd=1)
#sigex.graph(extract.sa,reg.trend,begin,period,subseries,0,sacol)
sigex.graph(extract.trendirreg,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trend,reg.trend,begin,period,subseries,0,trendcol,10)
sigex.graph(extract.seas.ann,NULL,begin,period,subseries,5,seascol,10)
sigex.graph(extract.seas.week,NULL,begin,period,subseries,3,cyccol,fade)

plot(data[,subseries],xlab="Year",ylab="",ylim=c(0,9),lwd=1)
sigex.graph(extract.seas.week1,NULL,begin,period,subseries,5,cyccol,fade)
sigex.graph(extract.seas.week2,NULL,begin,period,subseries,3,cyccol,fade)
sigex.graph(extract.seas.week3,NULL,begin,period,subseries,1,cyccol,fade)


subseries <- 4

plot(data[,subseries],xlab="Year",ylab="",ylim=c(2,9),lwd=1)
#sigex.graph(extract.sa,reg.trend,begin,period,subseries,0,sacol)
sigex.graph(extract.trendirreg,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trend,reg.trend,begin,period,subseries,0,trendcol,10)
sigex.graph(extract.seas.ann,NULL,begin,period,subseries,5,seascol,10)
sigex.graph(extract.seas.week,NULL,begin,period,subseries,3,cyccol,fade)

plot(data[,subseries],xlab="Year",ylab="",ylim=c(0,9),lwd=1)
sigex.graph(extract.seas.week1,NULL,begin,period,subseries,5,cyccol,fade)
sigex.graph(extract.seas.week2,NULL,begin,period,subseries,3,cyccol,fade)
sigex.graph(extract.seas.week3,NULL,begin,period,subseries,1,cyccol,fade)


subseries <- 5

plot(data[,subseries],xlab="Year",ylab="",ylim=c(2,7),lwd=1)
#sigex.graph(extract.sa,reg.trend,begin,period,subseries,0,sacol)
sigex.graph(extract.trendirreg,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trend,reg.trend,begin,period,subseries,0,trendcol,10)
sigex.graph(extract.seas.ann,NULL,begin,period,subseries,4,seascol,10)
sigex.graph(extract.seas.week,NULL,begin,period,subseries,3,cyccol,fade)

plot(data[,subseries],xlab="Year",ylab="",ylim=c(0,7),lwd=1)
sigex.graph(extract.seas.week1,NULL,begin,period,subseries,3,cyccol,fade)
sigex.graph(extract.seas.week2,NULL,begin,period,subseries,2,cyccol,fade)
sigex.graph(extract.seas.week3,NULL,begin,period,subseries,1,cyccol,fade)


subseries <- 6

plot(data[,subseries],xlab="Year",ylab="",ylim=c(2,7),lwd=1)
#sigex.graph(extract.sa,reg.trend,begin,period,subseries,0,sacol)
sigex.graph(extract.trendirreg,reg.trend,begin,period,subseries,0,sacol,fade)
sigex.graph(extract.trend,reg.trend,begin,period,subseries,0,trendcol,10)
sigex.graph(extract.seas.ann,NULL,begin,period,subseries,4,seascol,10)
sigex.graph(extract.seas.week,NULL,begin,period,subseries,3,cyccol,fade)

plot(data[,subseries],xlab="Year",ylab="",ylim=c(0,7),lwd=1)
sigex.graph(extract.seas.week1,NULL,begin,period,subseries,3,cyccol,fade)
sigex.graph(extract.seas.week2,NULL,begin,period,subseries,2,cyccol,fade)
sigex.graph(extract.seas.week3,NULL,begin,period,subseries,1,cyccol,fade)


#################################
### Seasonality Diagnostics

sigex.seascheck(extract.trend[[1]],7,.04,1)
sigex.seascheck(extract.seas.ann[[1]],7,.04,1)
sigex.seascheck(extract.seas.week[[1]],7,.04,0)
sigex.seascheck(extract.trendirreg[[1]],7,.04,1)

sigex.seascheck(extract.trend[[1]],365,.04,1)
sigex.seascheck(extract.seas.ann[[1]],365,.04,1)
sigex.seascheck(extract.seas.week[[1]],365,.04,0)
sigex.seascheck(extract.trendirreg[[1]],365,.04,1)

###### Signal extraction diagnostics

round(sigex.signalcheck(extract.trendann[[1]][(horizon+1):(horizon+T),],param,mdl,1,100),digits=3)
round(sigex.signalcheck(extract.seas.week[[1]],param,mdl,c(2,3,4),100),digits=3)

### Spectral density plots

subseries <- 1

## spectral plots
week.freq <- 365/7
month.freq <- 365/(365/12)
ann.freq <- 1

spec.ar(ts(extract.trend[[1]][,subseries],frequency=period),main="")
abline(v=ann.freq,col=4)
abline(v=week.freq,col=2)
abline(v=2*week.freq,col=2)
abline(v=3*week.freq,col=2)
abline(v=month.freq,col=3)

spec.ar(ts(extract.sa[[1]][,subseries],frequency=period),main="")
abline(v=ann.freq,col=4)
abline(v=week.freq,col=2)
abline(v=2*week.freq,col=2)
abline(v=3*week.freq,col=2)
abline(v=month.freq,col=3)

spec.ar(ts(extract.seas.week1[[1]][,subseries],frequency=period))
abline(v=week.freq,col=2)

spec.ar(ts(extract.seas.week2[[1]][,subseries],frequency=period))
abline(v=2*week.freq,col=2)

spec.ar(ts(extract.seas.week3[[1]][,subseries],frequency=period))
abline(v=3*week.freq,col=2)

##################################
#	filter analysis

grid <- 700
frf.trendann <- sigex.getfrf(data,param,mdl,1,TRUE,grid)
frf.week1 <- sigex.getfrf(data,param,mdl,2,TRUE,grid)
frf.week2 <- sigex.getfrf(data,param,mdl,3,TRUE,grid)
frf.week3 <- sigex.getfrf(data,param,mdl,4,TRUE,grid)
frf.weeks <- sigex.getfrf(data,param,mdl,c(2,3,4),TRUE,grid)
frf.irr <- sigex.getfrf(data,param,mdl,5,TRUE,grid)
frf.sa <- sigex.getfrf(data,param,mdl,c(1,5),TRUE,grid)

len <- 500
wk.trendann <- sigex.getwk(data,param,mdl,1,TRUE,grid,len)
wk.week1 <- sigex.getwk(data,param,mdl,2,TRUE,grid,len)
wk.week2 <- sigex.getwk(data,param,mdl,3,TRUE,grid,len)
wk.week3 <- sigex.getwk(data,param,mdl,4,TRUE,grid,len)
wk.weeks <- sigex.getwk(data,param,mdl,c(2,3,4),TRUE,grid,len)
wk.irr <- sigex.getwk(data,param,mdl,6,TRUE,grid,len)










#####  SCRAP

## fit 7 univariate models

# Sunday  SARMA

i <- 1
mdl.sun <- NULL
mdl.sun <- sigex.add(mdl.sun,1,"sarma",c(0,2,4,0,52),NULL,"process",c(1,-1))
mdl.sun <- sigex.meaninit(mdl.sun,data.ts[,i,drop=FALSE],0)
mdl.sun <- sigex.reg(mdl.sun,1,ts(as.matrix(easter.reg[,i]),
                            start=start(easter.reg),
                            frequency=frequency(easter.reg),
                            names="Easter-day"))
mdl.sun <- sigex.reg(mdl.sun,1,ts(as.matrix(school1.reg[,i]),
                            start=start(school1.reg),
                            frequency=frequency(school1.reg),
                            names="School1-Start"))
mdl.sun <- sigex.reg(mdl.sun,1,ts(as.matrix(school1e.reg[,i]),
                            start=start(school1e.reg),
                            frequency=frequency(school1e.reg),
                            names="School1-End"))
mdl.sun <- sigex.reg(mdl.sun,1,ts(as.matrix(school2.reg[,i]),
                            start=start(school2.reg),
                            frequency=frequency(school2.reg),
                            names="School2-Start"))
mdl.sun <- sigex.reg(mdl.sun,1,ts(as.matrix(school2e.reg[,i]),
                            start=start(school2e.reg),
                            frequency=frequency(school2e.reg),
                            names="School2-End"))
mdl.sun <- sigex.reg(mdl.sun,1,ts(as.matrix(school3.reg[,i]),
                            start=start(school3.reg),
                            frequency=frequency(school3.reg),
                            names="School3-Start"))
mdl.sun <- sigex.reg(mdl.sun,1,ts(as.matrix(school3e.reg[,i]),
                            start=start(school3e.reg),
                            frequency=frequency(school3e.reg),
                            names="School3-End"))

constraint <- NULL
par.mle <- sigex.default(mdl.sun,data.ts[,i,drop=FALSE],constraint)
psi.mle <- sigex.par2psi(par.mle,mdl.sun)

fit.mle <- sigex.mlefit(data.ts[,i,drop=FALSE],par.mle,constraint,mdl.sun,"bfgs",debug=TRUE)

psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]
tstats <- sigex.tstats(mdl.sun,psi.mle,hess,constraint)
print(tstats)

resid.mle <- sigex.resid(psi.mle,mdl.sun,data.ts[,i,drop=FALSE])[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),
                        colnames(data.ts[,i,drop=FALSE]),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*frequency(data.ts),plot=TRUE)$acf

# -2313.506
psi.sun <- psi.mle
psi.sun <- c(-4.07102751905699, 2.44060964365952, 0.504283072867637, 2.40960579766798,
0.284471757019792, 0.376982011579344, -0.170354318249828, 0.000821039182642259,
-0.356693054572869, -0.210338500206163, 0.178053589641864, -0.182968103386536,
0.22141697063521, -0.130994499688667, 0.263268549888182)

# Monday  SARMA

i <- 2
mdl.mon <- NULL
mdl.mon <- sigex.add(mdl.mon,1,"sarma",c(0,2,4,0,52),NULL,"process",c(1,-1))
mdl.mon <- sigex.meaninit(mdl.mon,data.ts[,i,drop=FALSE],0)
mdl.mon <- sigex.reg(mdl.mon,1,ts(as.matrix(easter.reg[,i]),
                                  start=start(easter.reg),
                                  frequency=frequency(easter.reg),
                                  names="Easter-day"))
mdl.mon <- sigex.reg(mdl.mon,1,ts(as.matrix(school1.reg[,i]),
                                  start=start(school1.reg),
                                  frequency=frequency(school1.reg),
                                  names="School1-Start"))
mdl.mon <- sigex.reg(mdl.mon,1,ts(as.matrix(school1e.reg[,i]),
                                  start=start(school1e.reg),
                                  frequency=frequency(school1e.reg),
                                  names="School1-End"))
mdl.mon <- sigex.reg(mdl.mon,1,ts(as.matrix(school2.reg[,i]),
                                  start=start(school2.reg),
                                  frequency=frequency(school2.reg),
                                  names="School2-Start"))
mdl.mon <- sigex.reg(mdl.mon,1,ts(as.matrix(school2e.reg[,i]),
                                  start=start(school2e.reg),
                                  frequency=frequency(school2e.reg),
                                  names="School2-End"))
mdl.mon <- sigex.reg(mdl.mon,1,ts(as.matrix(school3.reg[,i]),
                                  start=start(school3.reg),
                                  frequency=frequency(school3.reg),
                                  names="School3-Start"))
mdl.mon <- sigex.reg(mdl.mon,1,ts(as.matrix(school3e.reg[,i]),
                                  start=start(school3e.reg),
                                  frequency=frequency(school3e.reg),
                                  names="School3-End"))

constraint <- NULL
par.mle <- sigex.default(mdl.mon,data.ts[,i,drop=FALSE],constraint)
psi.mle <- sigex.par2psi(par.mle,mdl.mon)

fit.mle <- sigex.mlefit(data.ts[,i,drop=FALSE],par.mle,constraint,mdl.mon,"bfgs",debug=TRUE)

psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]
tstats <- sigex.tstats(mdl.mon,psi.mle,hess,constraint)
print(tstats)

resid.mle <- sigex.resid(psi.mle,mdl.mon,data.ts[,i,drop=FALSE])[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),
                        colnames(data.ts[,i,drop=FALSE]),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*frequency(data.ts),plot=TRUE)$acf

# -2058.759
psi.mon <- psi.mle
psi.mon <- c(-3.71673255276367, 2.81665989393415, 0.455821930160383, 1.95849931377876,
0.382618404217175, 0.347620963878928, -0.105333354153538, 0.00108330881916517,
0.1544093642847, 0.0571421511614232, 0.226684037278746, -0.337666113897167,
0.0955607194603728, 0.0871900840586311, 0.0355048167930293)

# Tuesday  SARMA

i <- 3
mdl.tue <- NULL
mdl.tue <- sigex.add(mdl.tue,1,"sarma",c(0,2,4,0,52),NULL,"process",c(1,-1))
mdl.tue <- sigex.meaninit(mdl.tue,data.ts[,i,drop=FALSE],0)
mdl.tue <- sigex.reg(mdl.tue,1,ts(as.matrix(easter.reg[,i]),
                                  start=start(easter.reg),
                                  frequency=frequency(easter.reg),
                                  names="Easter-day"))
mdl.tue <- sigex.reg(mdl.tue,1,ts(as.matrix(school1.reg[,i]),
                                  start=start(school1.reg),
                                  frequency=frequency(school1.reg),
                                  names="School1-Start"))
mdl.tue <- sigex.reg(mdl.tue,1,ts(as.matrix(school1e.reg[,i]),
                                  start=start(school1e.reg),
                                  frequency=frequency(school1e.reg),
                                  names="School1-End"))
mdl.tue <- sigex.reg(mdl.tue,1,ts(as.matrix(school2.reg[,i]),
                                  start=start(school2.reg),
                                  frequency=frequency(school2.reg),
                                  names="School2-Start"))
mdl.tue <- sigex.reg(mdl.tue,1,ts(as.matrix(school2e.reg[,i]),
                                  start=start(school2e.reg),
                                  frequency=frequency(school2e.reg),
                                  names="School2-End"))
mdl.tue <- sigex.reg(mdl.tue,1,ts(as.matrix(school3.reg[,i]),
                                  start=start(school3.reg),
                                  frequency=frequency(school3.reg),
                                  names="School3-Start"))
mdl.tue <- sigex.reg(mdl.tue,1,ts(as.matrix(school3e.reg[,i]),
                                  start=start(school3e.reg),
                                  frequency=frequency(school3e.reg),
                                  names="School3-End"))

constraint <- NULL
par.mle <- sigex.default(mdl.tue,data.ts[,i,drop=FALSE],constraint)
psi.mle <- sigex.par2psi(par.mle,mdl.tue)

fit.mle <- sigex.mlefit(data.ts[,i,drop=FALSE],par.mle,constraint,mdl.tue,"bfgs",debug=TRUE)

psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]
tstats <- sigex.tstats(mdl.tue,psi.mle,hess,constraint)
print(tstats)

resid.mle <- sigex.resid(psi.mle,mdl.tue,data.ts[,i,drop=FALSE])[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),
                        colnames(data.ts[,i,drop=FALSE]),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*frequency(data.ts),plot=TRUE)$acf

# -2072.44
psi.tue <- psi.mle
psi.tue <- c(-3.73650814903223, 3.38356662731394, 0.466866273929686, 1.78191727599679,
             0.55344678224801, 0.422070587987303, -0.123367710789483, 0.00105292031797421,
             0.328220766528598, -0.425224905861729, 0.461160570456221, -0.231132142135402,
             -0.204556232341297, 0.106365657597149, -1.08026453182959)


# Wednesday  SARMA

i <- 4
mdl.wed <- NULL
mdl.wed <- sigex.add(mdl.wed,1,"sarma",c(0,2,4,0,52),NULL,"process",c(1,-1))
mdl.wed <- sigex.meaninit(mdl.wed,data.ts[,i,drop=FALSE],0)
mdl.wed <- sigex.reg(mdl.wed,1,ts(as.matrix(easter.reg[,i]),
                                  start=start(easter.reg),
                                  frequency=frequency(easter.reg),
                                  names="Easter-day"))
mdl.wed <- sigex.reg(mdl.wed,1,ts(as.matrix(school1.reg[,i]),
                                  start=start(school1.reg),
                                  frequency=frequency(school1.reg),
                                  names="School1-Start"))
mdl.wed <- sigex.reg(mdl.wed,1,ts(as.matrix(school1e.reg[,i]),
                                  start=start(school1e.reg),
                                  frequency=frequency(school1e.reg),
                                  names="School1-End"))
mdl.wed <- sigex.reg(mdl.wed,1,ts(as.matrix(school2.reg[,i]),
                                  start=start(school2.reg),
                                  frequency=frequency(school2.reg),
                                  names="School2-Start"))
mdl.wed <- sigex.reg(mdl.wed,1,ts(as.matrix(school2e.reg[,i]),
                                  start=start(school2e.reg),
                                  frequency=frequency(school2e.reg),
                                  names="School2-End"))
mdl.wed <- sigex.reg(mdl.wed,1,ts(as.matrix(school3.reg[,i]),
                                  start=start(school3.reg),
                                  frequency=frequency(school3.reg),
                                  names="School3-Start"))
mdl.wed <- sigex.reg(mdl.wed,1,ts(as.matrix(school3e.reg[,i]),
                                  start=start(school3e.reg),
                                  frequency=frequency(school3e.reg),
                                  names="School3-End"))

constraint <- NULL
par.mle <- sigex.default(mdl.wed,data.ts[,i,drop=FALSE],constraint)
psi.mle <- sigex.par2psi(par.mle,mdl.wed)

fit.mle <- sigex.mlefit(data.ts[,i,drop=FALSE],par.mle,constraint,mdl.wed,"bfgs",debug=TRUE)

psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]
tstats <- sigex.tstats(mdl.wed,psi.mle,hess,constraint)
print(tstats)

resid.mle <- sigex.resid(psi.mle,mdl.wed,data.ts[,i,drop=FALSE])[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),
                        colnames(data.ts[,i,drop=FALSE]),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*frequency(data.ts),plot=TRUE)$acf

# -2187.547
psi.wed <- psi.mle
psi.wed <- c(-3.88752052496473, 2.57991874774606, 0.381616148947387, 1.7688876355792,
0.476209520549579, 0.566399196309179, -0.00854481257841716, 0.000893904427356007,
0.0492484190910288, 0.0166551032816581, -1.09965928066817, 0.507858990853311,
-0.645837739465141, 0.490057288297561, -0.642952148261199)


#  Thursday  SARMA

i <- 5
mdl.thu <- NULL
mdl.thu <- sigex.add(mdl.thu,1,"sarma",c(0,2,4,0,52),NULL,"process",c(1,-1))
mdl.thu <- sigex.meaninit(mdl.thu,data.ts[,i,drop=FALSE],0)
mdl.thu <- sigex.reg(mdl.thu,1,ts(as.matrix(easter.reg[,i]),
                                  start=start(easter.reg),
                                  frequency=frequency(easter.reg),
                                  names="Easter-day"))
mdl.thu <- sigex.reg(mdl.thu,1,ts(as.matrix(school1.reg[,i]),
                                  start=start(school1.reg),
                                  frequency=frequency(school1.reg),
                                  names="School1-Start"))
mdl.thu <- sigex.reg(mdl.thu,1,ts(as.matrix(school1e.reg[,i]),
                                  start=start(school1e.reg),
                                  frequency=frequency(school1e.reg),
                                  names="School1-End"))
mdl.thu <- sigex.reg(mdl.thu,1,ts(as.matrix(school2.reg[,i]),
                                  start=start(school2.reg),
                                  frequency=frequency(school2.reg),
                                  names="School2-Start"))
mdl.thu <- sigex.reg(mdl.thu,1,ts(as.matrix(school2e.reg[,i]),
                                  start=start(school2e.reg),
                                  frequency=frequency(school2e.reg),
                                  names="School2-End"))
mdl.thu <- sigex.reg(mdl.thu,1,ts(as.matrix(school3.reg[,i]),
                                  start=start(school3.reg),
                                  frequency=frequency(school3.reg),
                                  names="School3-Start"))
mdl.thu <- sigex.reg(mdl.thu,1,ts(as.matrix(school3e.reg[,i]),
                                  start=start(school3e.reg),
                                  frequency=frequency(school3e.reg),
                                  names="School3-End"))

constraint <- NULL
par.mle <- sigex.default(mdl.thu,data.ts[,i,drop=FALSE],constraint)
psi.mle <- sigex.par2psi(par.mle,mdl.thu)

fit.mle <- sigex.mlefit(data.ts[,i,drop=FALSE],par.mle,constraint,mdl.thu,"bfgs",debug=TRUE)

psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]
tstats <- sigex.tstats(mdl.thu,psi.mle,hess,constraint)
print(tstats)

resid.mle <- sigex.resid(psi.mle,mdl.thu,data.ts[,i,drop=FALSE])[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),
                        colnames(data.ts[,i,drop=FALSE]),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*frequency(data.ts),plot=TRUE)$acf

# -2057.218
psi.thu <- psi.mle
psi.thu <- c(-3.72998751169167, 2.50300536739333, 0.383291814679505, 1.99917886503278,
             0.518392117669723, 0.407488351143634, -0.103163118274746, 0.000640049582275312,
             -0.0573437351433279, 0.497956950408837, 1.50904671713997, 0.311572498937319,
             0.656559136345757, 0.055313659778188, 0.672494154045533)

#  Friday  SARMA

i <- 6
mdl.fri <- NULL
mdl.fri <- sigex.add(mdl.fri,1,"sarma",c(0,2,4,0,52),NULL,"process",c(1,-1))
mdl.fri <- sigex.meaninit(mdl.fri,data.ts[,i,drop=FALSE],0)
mdl.fri <- sigex.reg(mdl.fri,1,ts(as.matrix(easter.reg[,i]),
                                  start=start(easter.reg),
                                  frequency=frequency(easter.reg),
                                  names="Easter-day"))
mdl.fri <- sigex.reg(mdl.fri,1,ts(as.matrix(school1.reg[,i]),
                                  start=start(school1.reg),
                                  frequency=frequency(school1.reg),
                                  names="School1-Start"))
mdl.fri <- sigex.reg(mdl.fri,1,ts(as.matrix(school1e.reg[,i]),
                                  start=start(school1e.reg),
                                  frequency=frequency(school1e.reg),
                                  names="School1-End"))
mdl.fri <- sigex.reg(mdl.fri,1,ts(as.matrix(school2.reg[,i]),
                                  start=start(school2.reg),
                                  frequency=frequency(school2.reg),
                                  names="School2-Start"))
mdl.fri <- sigex.reg(mdl.fri,1,ts(as.matrix(school2e.reg[,i]),
                                  start=start(school2e.reg),
                                  frequency=frequency(school2e.reg),
                                  names="School2-End"))
mdl.fri <- sigex.reg(mdl.fri,1,ts(as.matrix(school3.reg[,i]),
                                  start=start(school3.reg),
                                  frequency=frequency(school3.reg),
                                  names="School3-Start"))
mdl.fri <- sigex.reg(mdl.fri,1,ts(as.matrix(school3e.reg[,i]),
                                  start=start(school3e.reg),
                                  frequency=frequency(school3e.reg),
                                  names="School3-End"))

constraint <- NULL
par.mle <- sigex.default(mdl.fri,data.ts[,i,drop=FALSE],constraint)
psi.mle <- sigex.par2psi(par.mle,mdl.fri)

fit.mle <- sigex.mlefit(data.ts[,i,drop=FALSE],par.mle,constraint,mdl.fri,"bfgs",debug=TRUE)

psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]
tstats <- sigex.tstats(mdl.fri,psi.mle,hess,constraint)
print(tstats)

resid.mle <- sigex.resid(psi.mle,mdl.fri,data.ts[,i,drop=FALSE])[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),
                        colnames(data.ts[,i,drop=FALSE]),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*frequency(data.ts),plot=TRUE)$acf

# -2074.888
psi.fri <- psi.mle
psi.fri <- c(-3.75835853837978, 2.92661574783167, 0.152738271758969, 2.13559818245616,
             0.463643958101823, 0.528739568613193, -0.0597676412636782, 0.000326919006301636,
             0.0750839533564237, -0.382807515259472, 1.14385613449601, -0.436345776605543,
             0.0393416191070808, 0.310748871253615, 0.153824128601056)

# Saturday  SARMA

i <- 7
mdl.sat <- NULL
mdl.sat <- sigex.add(mdl.sat,1,"sarma",c(0,2,4,0,52),NULL,"process",c(1,-1))
mdl.sat <- sigex.meaninit(mdl.sat,data.ts[,i,drop=FALSE],0)
mdl.sat <- sigex.reg(mdl.sat,1,ts(as.matrix(easter.reg[,i]),
                                  start=start(easter.reg),
                                  frequency=frequency(easter.reg),
                                  names="Easter-day"))
mdl.sat <- sigex.reg(mdl.sat,1,ts(as.matrix(school1.reg[,i]),
                                  start=start(school1.reg),
                                  frequency=frequency(school1.reg),
                                  names="School1-Start"))
mdl.sat <- sigex.reg(mdl.sat,1,ts(as.matrix(school1e.reg[,i]),
                                  start=start(school1e.reg),
                                  frequency=frequency(school1e.reg),
                                  names="School1-End"))
mdl.sat <- sigex.reg(mdl.sat,1,ts(as.matrix(school2.reg[,i]),
                                  start=start(school2.reg),
                                  frequency=frequency(school2.reg),
                                  names="School2-Start"))
mdl.sat <- sigex.reg(mdl.sat,1,ts(as.matrix(school2e.reg[,i]),
                                  start=start(school2e.reg),
                                  frequency=frequency(school2e.reg),
                                  names="School2-End"))
mdl.sat <- sigex.reg(mdl.sat,1,ts(as.matrix(school3.reg[,i]),
                                  start=start(school3.reg),
                                  frequency=frequency(school3.reg),
                                  names="School3-Start"))
mdl.sat <- sigex.reg(mdl.sat,1,ts(as.matrix(school3e.reg[,i]),
                                  start=start(school3e.reg),
                                  frequency=frequency(school3e.reg),
                                  names="School3-End"))

constraint <- NULL
par.mle <- sigex.default(mdl.sat,data.ts[,i,drop=FALSE],constraint)
psi.mle <- sigex.par2psi(par.mle,mdl.sat)

fit.mle <- sigex.mlefit(data.ts[,i,drop=FALSE],par.mle,constraint,mdl.sat,"bfgs",debug=TRUE)

psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]
tstats <- sigex.tstats(mdl.sat,psi.mle,hess,constraint)
print(tstats)

resid.mle <- sigex.resid(psi.mle,mdl.sat,data.ts[,i,drop=FALSE])[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),
                        colnames(data.ts[,i,drop=FALSE]),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*frequency(data.ts),plot=TRUE)$acf

#  -2050.405
psi.sat <- psi.mle
psi.sat <- c(-3.72122242545419, 2.84637162942386, 0.481698020325866, 1.95540213883444,
             0.54614410843796, 0.513339642154714, 0.0526384416356759, 0.000932076031292542,
             -0.245759680535988, -0.199815969680003, 0.435888749475343, -0.0960892559637839,
             0.399847261433419, -0.116661022378742, 0.434501022976733)


##
psi.days <- NULL
psi.days <- cbind(psi.days,psi.sun)
psi.days <- cbind(psi.days,psi.mon)
psi.days <- cbind(psi.days,psi.tue)
psi.days <- cbind(psi.days,psi.wed)
psi.days <- cbind(psi.days,psi.thu)
psi.days <- cbind(psi.days,psi.fri)
psi.days <- cbind(psi.days,psi.sat)
psi.days <- psi.days[1:7,]

psi.days <- structure(c(-4.07102751905699, 2.44060964365952, 0.504283072867637,
                        2.40960579766798, 0.284471757019792, 0.376982011579344, -0.170354318249828,
                        -3.71673255276367, 2.81665989393415, 0.455821930160383, 1.95849931377876,
                        0.382618404217175, 0.347620963878928, -0.105333354153538, -3.73650814903223,
                        3.38356662731394, 0.466866273929686, 1.78191727599679, 0.55344678224801,
                        0.422070587987303, -0.123367710789483, -3.88752052496473, 2.57991874774606,
                        0.381616148947387, 1.7688876355792, 0.476209520549579, 0.566399196309179,
                        -0.00854481257841716, -3.72998751169167, 2.50300536739333, 0.383291814679505,
                        1.99917886503278, 0.518392117669723, 0.407488351143634, -0.103163118274746,
                        -3.75835853837978, 2.92661574783167, 0.152738271758969, 2.13559818245616,
                        0.463643958101823, 0.528739568613193, -0.0597676412636782, -3.72122242545419,
                        2.84637162942386, 0.481698020325866, 1.95540213883444, 0.54614410843796,
                        0.513339642154714, 0.0526384416356759), .Dim = c(7L, 7L), .Dimnames = list(
                          NULL, c("psi.sun", "psi.mon", "psi.tue", "psi.wed", "psi.thu",
                                  "psi.fri", "psi.sat")))

# SARMA
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"sarma",c(0,2,4,0,52),NULL,"process",c(1,-1))
mdl <- sigex.meaninit(mdl,data.ts,0)
for(i in 1:N) {
  mdl <- sigex.reg(mdl,i,ts(as.matrix(easter.reg[,i]),
                            start=start(easter.reg),
                            frequency=frequency(easter.reg),
                            names="Easter-day"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school1.reg[,i]),
                            start=start(school1.reg),
                            frequency=frequency(school1.reg),
                            names="School1-Start"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school1e.reg[,i]),
                            start=start(school1e.reg),
                            frequency=frequency(school1e.reg),
                            names="School1-End"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school2.reg[,i]),
                            start=start(school2.reg),
                            frequency=frequency(school2.reg),
                            names="School2-Start"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school2e.reg[,i]),
                            start=start(school2e.reg),
                            frequency=frequency(school2e.reg),
                            names="School2-End"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school3.reg[,i]),
                            start=start(school3.reg),
                            frequency=frequency(school3.reg),
                            names="School3-Start"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school3e.reg[,i]),
                            start=start(school3e.reg),
                            frequency=frequency(school3e.reg),
                            names="School3-End"))
}

## parameter initialization
constraint <- NULL

# regression constraints
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(2,2,2,2,2,2,2),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(3,3,3,3,3,3,3),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(4,4,4,4,4,4,4),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(5,5,5,5,5,5,5),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(6,6,6,6,6,6,6),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(7,7,7,7,7,7,7),NULL))

par.init <- sigex.default(mdl,data.ts,constraint)
psi.init <- sigex.par2psi(par.init,mdl)

psi.init[22:28] <- psi.days[1,]
psi.init[29:42] <- matrix(psi.days[2:3,],ncol=1)
psi.init[43:70] <- matrix(psi.days[4:7,],ncol=1)
par.mle <- sigex.psi2par(psi.init,mdl,data.ts)


sigex.lik(psi.init,mdl,data.ts)

fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)




###  getting constrained SVARMA


# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"svarma",c(0,2,4,0,52),NULL,"process",c(1,-1))
mdl <- sigex.meaninit(mdl,data.ts,0)

for(i in 1:N) {
  mdl <- sigex.reg(mdl,i,ts(as.matrix(easter.reg[,i]),
                            start=start(easter.reg),
                            frequency=frequency(easter.reg),
                            names="Easter-day"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school1.reg[,i]),
                            start=start(school1.reg),
                            frequency=frequency(school1.reg),
                            names="School1-Start"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school1e.reg[,i]),
                            start=start(school1e.reg),
                            frequency=frequency(school1e.reg),
                            names="School1-End"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school2.reg[,i]),
                            start=start(school2.reg),
                            frequency=frequency(school2.reg),
                            names="School2-Start"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school2e.reg[,i]),
                            start=start(school2e.reg),
                            frequency=frequency(school2e.reg),
                            names="School2-End"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school3.reg[,i]),
                            start=start(school3.reg),
                            frequency=frequency(school3.reg),
                            names="School3-Start"))
  mdl <- sigex.reg(mdl,i,ts(as.matrix(school3e.reg[,i]),
                            start=start(school3e.reg),
                            frequency=frequency(school3e.reg),
                            names="School3-End"))
}


constraint <- NULL

# model constraints and initial values
par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)
psi.len <- length(psi.mle)
constraint <- diag(psi.len)
index.mat <- matrix(seq(1,N^2),nrow=N,ncol=N)
index.constraint <- setdiff(matrix(index.mat,ncol=1),diag(index.mat))
offset <- N*(N+1)/2
psi.mle[offset + diag(index.mat)] <- psi.days[2,]
index.constraints <- offset + index.constraint
offset <- offset + N^2
psi.mle[offset + diag(index.mat)] <- psi.days[3,]
index.constraints <- c(index.constraints,offset + index.constraint)
offset <- offset + N^2
psi.mle[offset + diag(index.mat)] <- psi.days[4,]
index.constraints <- c(index.constraints,offset + index.constraint)
offset <- offset + N^2
psi.mle[offset + diag(index.mat)] <- psi.days[5,]
index.constraints <- c(index.constraints,offset + index.constraint)
offset <- offset + N^2
psi.mle[offset + diag(index.mat)] <- psi.days[6,]
index.constraints <- c(index.constraints,offset + index.constraint)
offset <- offset + N^2
psi.mle[offset + diag(index.mat)] <- psi.days[7,]
index.constraints <- c(index.constraints,offset + index.constraint)
constraint <- constraint[index.constraints,,drop=FALSE]
constraint <- cbind(rep(0,length(index.constraints)),constraint)

# regression constraints
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(2,2,2,2,2,2,2),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(3,3,3,3,3,3,3),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(4,4,4,4,4,4,4),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(5,5,5,5,5,5,5),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(6,6,6,6,6,6,6),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(7,7,7,7,7,7,7),NULL))

fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)











