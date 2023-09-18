#########################################
###  Script for Daily Immigration Data
#########################################

## wipe
rm(list=ls())

library(devtools)
library(Rcpp)

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
dataALL.ts <- sigex.load(imm,begin,period,
                         c("NZArr","NZDep","VisArr","VisDep","PLTArr","PLTDep"),TRUE)


#############################
## select span and transforms

## first series with log transform
transform <- "log"
aggregate <- FALSE
subseries <- 1
range <- list(c(2008,1),end)
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

times <- seq((range[[1]][1]-begin[1])*period+(range[[1]][2]-begin[2])+1,
             (range[[2]][1]-begin[1])*period+(range[[2]][2]-begin[2])+1,1)

easter.reg <- NZregs[times,1]
school1.reg <- NZregs[times,2]
school1e.reg <- NZregs[times,3]
school2.reg <- NZregs[times,4]
school2e.reg <- NZregs[times,5]
school3.reg <- NZregs[times,6]
school3e.reg <- NZregs[times,7]

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
mdl <- sigex.add(mdl,seq(1,N),"svarma",c(0,1,1,0,52),NULL,"process",c(1,rep(0,51),-1))
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

# regression constraints
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(2,2,2,2,2,2,2),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(3,3,3,3,3,3,3),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(4,4,4,4,4,4,4),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(5,5,5,5,5,5,5),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(6,6,6,6,6,6,6),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(7,7,7,7,7,7,7),NULL))

par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting: commented out, this took a long time
#fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## input parameter from previous fit (MLE on reduced span)
# divergence -4129.609
psi.mle <- c(0.195543136661102, 0.210896311604314, 0.224884839193806, 0.297710366252992,
             0.318076964322756, 0.288167112704856, 0.130848490409544, 0.17299456406149,
             0.213720675984967, 0.0780486278573135, 0.188624612762698, 0.467926733192453,
             0.434759495960342, 0.314822770386916, 0.280937951321511, 0.484722644882768,
             0.413444260303994, 0.193778981158949, 0.550345260075344, 0.621238984502903,
             0.563243057038342, -3.99208014979805, -4.19642305454813, -3.8415818035668,
             -4.32856255035426, -3.95499509133263, -3.96279211490641, -4.03911188215587,
             -0.369778097197869, -0.0318374973193093, -0.0459307709775452,
             -0.138903858423006, -0.10135335156786, -0.2142068950075, -0.112796685064521,
             -0.140796608637905, -0.237750837664742, -0.194243215055276, 0.0219933280900727,
             -0.0185736546749142, 0.139055815995136, 0.185830059881372, 0.0297317536468041,
             -0.0155812148508496, -0.0671396386183212, 0.0476284405203113,
             -0.0570660189385589, -0.155367001953893, -0.198330717820409,
             -0.0295063266207021, -0.073765478114995, -0.0829227752708303,
             -0.187804918878073, -0.0422869655675992, -0.0251387557650662,
             0.0656914410289482, -0.0840920843255883, -0.145935243196311,
             -0.0687523983229335, -0.117014892923535, -0.144653364980642,
             -0.0790943505687107, 0.0611829711400193, -0.104109662533708,
             0.112359996172669, -0.124553026055195, -0.163431409482439, -0.0283323592215797,
             -0.0312605011210402, -0.178946106450384, -0.266581746889425,
             -0.164720047360709, 0.00959770743708815, 0.159928359498417, -0.0132411806278071,
             -0.0642101491748187, -0.143487318464752, -0.086691994085389,
             0.0814467471639927, -0.0123455935133012, 0.00908782055333108,
             -0.0441428053848971, -0.0194126102991047, -0.00524075702955244,
             -0.0105279430631441, -0.00420628604690744, 0.0789699162821441,
             -0.015850906929479, -0.043499933859938, -0.0062803013365077,
             -0.0207414587537556, -0.0459309858718756, 0.0445914128582114,
             -0.0946876132299148, 0.00588335338136514, 0.00922656520387984,
             -0.0140866870525241, -0.0157842160720905, -0.0319375488373468,
             0.0514126088016708, 0.0141484170525309, -0.124101882618575, -0.00774709495827467,
             0.0255995670388491, -0.0364850207049338, 0.0230572817360244,
             0.0220446070005591, -0.000478476269822034, -0.0168616246820038,
             -0.102427194747844, 0.00112763406309242, -0.0259316390382219,
             -0.010758329862067, -0.00198961225224573, -0.0464823007621206,
             -0.0165414789720724, 0.0104091697972505, -0.14827440476309, 0.0505648610397286,
             -0.0183900737977979, 0.0311301005304512, 0.00518380184796121,
             0.0104592851480601, -0.0109208667318225, -0.0637859237889965,
             -0.0656126331052389, 0.00033010004554818, -0.191533893764455,
             -0.157111537080312, 0.141669136870683, -0.0225792281088362, 0.157100658532765,
             -0.08441684545665, 0.1847632577291, 0.000561604448005, -0.191533893764455,
             -0.157111537080311, 0.141669136870683, -0.0225792281088361, 0.157100658532765,
             -0.0844168454566499, 0.183411740706148, 0.000312708939566006,
             -0.191533893764455, -0.157111537080312, 0.141669136870683, -0.0225792281088362,
             0.157100658532765, -0.08441684545665, -0.00163601810705088, 0.00100594355939265,
             -0.191533893764455, -0.157111537080311, 0.141669136870683, -0.0225792281088362,
             0.157100658532765, -0.08441684545665, 0.0156601320968592, 9.68417444876468e-05,
             -0.191533893764455, -0.157111537080312, 0.141669136870683, -0.0225792281088361,
             0.157100658532765, -0.08441684545665, 0.00328023110252021, 0.000461006050610564,
             -0.191533893764455, -0.157111537080312, 0.141669136870683, -0.0225792281088361,
             0.157100658532765, -0.08441684545665, 0.00430421751741802, 0.000496876409071807,
             -0.191533893764455, -0.157111537080311, 0.141669136870683, -0.0225792281088362,
             0.157100658532765, -0.08441684545665, 0.0023532873975757)
par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

##  model checking
resid.mle <- sigex.resid(psi.mle,mdl,data.ts)[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*53,plot=FALSE)$acf

## examine condition numbers
log(sigex.conditions(data.ts,psi.mle,mdl))

## model checking
sigex.portmanteau(resid.mle,2*52,length(psi.mle))
sigex.gausscheck(resid.mle)

#pdf(file="nzResidAcf.pdf",height=10,width=10)
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

## bundle
analysis.mle <- sigex.bundle(data.ts,transform,mdl,psi.mle)


##########################################
### Part V: Signal Extraction

## load up the fitted model for signal extraction
data.ts <- analysis.mle[[1]]
mdl <- analysis.mle[[3]]
psi <- analysis.mle[[4]]
param <- sigex.psi2par(psi,mdl,data.ts)

## embed daily SA filter as a weekly filter
sa.hifilter <- c(1,rep(2,365),1)/(2*365)
len <- 183
hi.freq <- 7
low.freq <- 1
shift.hi <- len
out <- sigex.hi2low(sa.hifilter,hi.freq,low.freq,shift.hi)
sa.lowfilter <- out[[1]]
shift.low <- out[[2]]

sa.low <- sigex.adhocextract(psi,mdl,data.ts,sa.lowfilter,shift.low,0,TRUE)
sa.hi.daily <- list()
sa.hi.daily[[1]] <- sigex.weekly2daily(ts(sa.low[[1]],start=start(dataONE.ts),
                                           frequency=frequency(dataONE.ts)),first.day)
sa.hi.daily[[2]] <- sigex.weekly2daily(ts(sa.low[[2]],start=start(dataONE.ts),
                                           frequency=frequency(dataONE.ts)),first.day)
sa.hi.daily[[3]] <- sigex.weekly2daily(ts(sa.low[[3]],start=start(dataONE.ts),
                                           frequency=frequency(dataONE.ts)),first.day)

## embed daily TD filter as a weekly filter
td.hifilter <- rep(1,7)/7
len <- 3
hi.freq <- 7
low.freq <- 1
shift.hi <- len
out <- sigex.hi2low(td.hifilter,hi.freq,low.freq,shift.hi)
td.lowfilter <- out[[1]]
shift.low <- out[[2]]

td.low <- sigex.adhocextract(psi,mdl,data.ts,td.lowfilter,shift.low,0,TRUE)
td.hi.daily <- list()
td.hi.daily[[1]] <- sigex.weekly2daily(ts(td.low[[1]],start=start(dataONE.ts),
                                           frequency=frequency(dataONE.ts)),first.day)
td.hi.daily[[2]] <- sigex.weekly2daily(ts(td.low[[2]],start=start(dataONE.ts),
                                           frequency=frequency(dataONE.ts)),first.day)
td.hi.daily[[3]] <- sigex.weekly2daily(ts(td.low[[3]],start=start(dataONE.ts),
                                           frequency=frequency(dataONE.ts)),first.day)

## get fixed effects
reg.td <- NULL
for(i in 1:N)
{
  reg.td <- cbind(reg.td,sigex.fixed(data.ts,mdl,i,param,"Trend"))
}
reg.td <- rbind(rep(0,N),reg.td)
reg.td <- ts(sigex.weekly2daily(reg.td,first.day),
                start=start(dataONE.ts),frequency=period)
reg.trend <- stats::filter(reg.td,td.hifilter,method="convolution",sides=1)
reg.trend <- as.matrix(reg.trend[8:length(reg.td)])


## plotting
trendcol <- "tomato"
cyccol <- "orchid"
seascol <- "seagreen"
sacol <- "navyblue"
fade <- 60

#pdf(file="nz-signals.pdf",height=8,width=10)
plot(dataONE.ts,xlab="Year")
sigex.graph(sa.hi.daily,reg.trend,start(sa.hi.daily[[1]]),
            period,1,0,trendcol,fade)
sigex.graph(td.hi.daily,reg.trend,start(td.hi.daily[[1]]),
            period,1,0,sacol,fade)
dev.off()

## spectral diagnostics: seasonal adjustment
sigex.specar(sa.hi.daily[[1]],FALSE,1,7)
dev.off()

## spectral diagnostics: non-weekly effect
sigex.specar(td.hi.daily[[1]],FALSE,1,7)
dev.off()




