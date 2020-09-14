#########################################
###  Script for Daily Immigration Data
#########################################

## wipe
rm(list=ls())

library(devtools)

# suppose directory is set to where sigex is located, e.g.
#setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex")
load_all(".")
root.dir <- getwd()
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
offset <- N*(N+1)/2 + N^2
index.constraints <- offset + index.constraint
offset <- offset + N^2
index.constraints <- c(index.constraints,offset + index.constraint)
offset <- offset + N^2
index.constraints <- c(index.constraints,offset + index.constraint)
constraint <- constraint[index.constraints,,drop=FALSE]
constraint <- cbind(constraint,rep(0,length(index.constraints)))

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







##### SCRAP

series <- 7

# SARMA
mdl2 <- NULL
mdl2 <- sigex.add(mdl2,1,"sarma",c(0,3,3,0,52),NULL,"process",c(1,-1))
mdl2 <- sigex.meaninit(mdl2,data.ts[,series],0)


par.mle <- sigex.default(mdl2,data.ts[,series,drop=FALSE],constraint)
psi.mle <- sigex.par2psi(par.mle,mdl2)

fit.mle2 <- sigex.mlefit(data.ts[,series,drop=FALSE],par.mle,constraint,mdl2,"bfgs",debug=TRUE)

## manage output
psi.mle <- sigex.eta2psi(fit.mle2[[1]]$par,constraint)
hess <- fit.mle2[[1]]$hessian
par.mle <- fit.mle2[[2]]


resid.mle <- sigex.resid(psi.mle,mdl2,data.ts[,series,drop=FALSE])[[1]]
resid.acf <- acf(t(resid.mle),lag.max=4*53,plot=TRUE)$acf


print(eigen(hess)$values)
tstats <- sigex.tstats(mdl,psi.mle,hess,constraint)
print(tstats)





