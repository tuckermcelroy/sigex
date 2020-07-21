#########################################
###  Script for Daily Immigration Data
#########################################

## wipe
rm(list=ls())

library(devtools)

#setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex")
setwd("/home/tucker/Documents/GitHub/sigex")
load_all(".")


######################
### Part I: load data

# automatic: raw data

# processing

n.months <- dim(imm)[1]/32
imm <- imm[-seq(1,n.months)*32,]	# strip out every 32nd row (totals)
imm <- matrix(imm[imm != 0],ncol=6) # strip out 31st False days

# enter regressors
NZregs <- read.table("C:\\Users\\neide\\Documents\\GitHub\\sigex\\data\\NZregressors.dat")



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

## all data with log transform
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

# model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"svarma",c(1,0,1,0,52),list(-1,1,1,1),"process",c(1,-1))
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
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(2,2,2,2,2,2,2),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(3,3,3,3,3,3,3),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(4,4,4,4,4,4,4),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(5,5,5,5,5,5,5),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(6,6,6,6,6,6,6),NULL))
constraint <- rbind(constraint,sigex.constrainreg(mdl,data.ts,list(7,7,7,7,7,7,7),NULL))

par.mle <- sigex.default(mdl,data.ts,constraint)
psi.mle <- sigex.par2psi(par.mle,mdl)

## run fitting: commented out, this took 2 weeks!
fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

## fit from first model
psi.mle <- c(0.51760596525827, 0.439939433340027, 0.341971516664744, 0.388139046943854,
             0.418853891103159, 0.464968608915972, 0.709312914960205, 0.455461460979668,
             0.339370325959773, 0.25187135366609, 0.230731337783255, 0.480592700293915,
             0.461178712098393, 0.420665953504732, 0.407603454396325, 0.654165574391287,
             0.616559331219179, 0.45821247032649, 0.559087776741675, 0.501500545560296,
             0.624756459479047, -4.11411535390305, -3.86649122997175, -4.12282006467438,
             -4.14352251305745, -4.09256331893617, -4.26412180312553, -3.99327709096135,
             0.45125721213338, 0.32353976331698, 0.422326213267438, 0.416978158837141,
             0.257484666829598, -0.0142433335792283, 0.306524358015199, 0.33034741287189,
             0.000434708380751397, 0.335362522632065, 0.313715126963322, 0.361026511518399,
             -0.00887460370005155, 0.41695106042505, -0.0555470870701567,
             1.58499541202741, -0.410063873377494, 0.612008573274832, 2.00048468352641,
             -0.538258408076739, 1.38642829853969, -0.7497661081984, -0.832273525977879,
             -1.21633288748606, -2.82554272049041, -2.86416045493807, -2.22595183707097,
             -3.03433319752487, -0.989660885819353, 2.14262898292818, -0.18441515172129,
             -1.82892437934759, 1.94909433314714, 0.128364625306125, -1.43869287117603,
             1.52886606616434, 1.87180695085432, -2.88271632918305, 2.03328110157552,
             0.429369463208369, 2.86434929128214, 1.83409415248455, -1.57900145879439,
             2.26586325285055, 0.891976969340423, 0.0337802195054125, 2.08630748199078,
             0.457086510632105, 2.1429409355955, 0.679651496181021, 0.329100762213215,
             0.239004399634381, 0.0969786787699651, 0.0242212597057822, 0.0640763528019264,
             0.921253414189022, 0.456510060844654, 0.307433277145931, 0.202710370503894,
             -0.226408189975761, 0.917582287859971, 0.101753218280991, 0.081507792065844,
             0.546418370349773, 1.15517525089724, 1.17422164159534, 0.726734392762665,
             0.948405703934301, 0.533468840769047, 0.654425583446069, -0.59163629141574,
             -1.8152232389266, -3.06728487190952, -2.67699952108507, -2.82598344743467,
             -3.35455145947087, -2.59503943987494, -0.0687073844054769, 0.119918448031456,
             -0.20694627785148, -0.00805952041715179, 0.315412704624931, -0.278332547979233,
             -0.430511041741036, 0.357189887132695, -0.0405823036457195, -0.546600023580311,
             0.387228814389113, 0.0227289586461886, -0.164772632575552, 0.628780413349036,
             -0.431156499527014, -0.185412934756364, -1.05055876957074, 0.722660829455804,
             -0.0394027990246995, -0.207072466047046, -0.163698549978204,
             0.00215074003032701, -0.266297522612305, -0.137093567965138,
             0.0851546535686134, -0.0149444051946342, 0.0358090120843427,
             -7.65592755685606e-05, 0.0364386379329101, 0.00169685678096456,
             -0.266297522612304, -0.137093567965138, 0.0851546535686135, -0.0149444051946342,
             0.0358090120843427, -7.6559275568536e-05, 0.15642984342176, 0.00070436907534173,
             -0.266297522612305, -0.137093567965139, 0.0851546535686135, -0.0149444051946342,
             0.0358090120843427, -7.65592755685972e-05, -0.172455161593865,
             -0.00233158509725544, -0.266297522612305, -0.137093567965138,
             0.0851546535686134, -0.0149444051946342, 0.0358090120843427,
             -7.65592755685624e-05, -0.0166107413493935, -0.00296386973878176,
             -0.266297522612305, -0.137093567965138, 0.0851546535686134, -0.0149444051946342,
             0.0358090120843427, -7.65592755685521e-05, 0.0976909319888783,
             -0.00559000446619553, -0.266297522612305, -0.137093567965138,
             0.0851546535686134, -0.0149444051946342, 0.0358090120843427,
             -7.6559275568555e-05, -0.0768525179813562, -0.00368827372534544,
             -0.266297522612305, -0.137093567965138, 0.0851546535686133, -0.0149444051946342,
             0.0358090120843427, -7.65592755685606e-05, 0.0489376995896597
)
par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

## fit from 2nd model so far...
psi.mle <- c(0.445988564485456, 0.315737077279683, 0.319065844461127, 0.265417327908842,
             0.293896378255691, 0.291346057216799, 0.70299159411062, 0.422290481160688,
             0.307524424107244, 0.212232156395402, 0.183559409008681, 0.476372259158625,
             0.472443239864719, 0.399499785255757, 0.358261834646677, 0.678581108779599,
             0.583392273244843, 0.399097122136101, 0.555071225260371, 0.427259073347545,
             0.633160903677843, -3.76755333451245, -3.84656330133457, -4.11690453089444,
             -4.15666647609728, -4.13286503651802, -4.29150753241301, -4.03089835155497,
             1.67092226019338, 0.679063091300303, 0.694118644440229, -0.596839619515064,
             0.0420262434430208, -1.38943494122707, 0.372613520795414, 0.84387011119455,
             -0.07133847754067, 0.450957680010647, 0.537170280127437, 0.436669426950197,
             -0.140872539774387, 0.314459395199935, 0.0224523585219865, 1.9327703903551,
             -0.483930132457992, -0.284457599293032, 2.07878838610121, -0.523581768536943,
             0.463638715839605, -3.0763164717996, -1.35464106001759, -1.22719591140301,
             -3.13270623939232, -3.10588248552621, -2.80077727800692, -3.43016135173596,
             -1.03751137562884, 2.37815800893698, -1.34999431179513, -1.7974281778018,
             1.59907711836505, 0.174793210256249, -1.46055921060779, 2.05910929993713,
             1.4351570265233, -3.11315282125059, 2.34641112770724, -1.8728471411523,
             3.05233463450168, 0.502365433634216, -2.09566276465018, 1.38609035755964,
             1.9867658756299, 1.04385596772587, 2.24166412912973, 0.48301016971414,
             2.95924449985295, 0.815309389798107, 0.547408556116637, 0.402658803643701,
             0.347904243958205, 0.293148591169946, 0.33381302448024, 0.9872125457884,
             0.580921887380074, 0.448280970418907, 0.417370784017859, 0.259186993195619,
             0.832639703107403, 0.214921383832783, 0.296985946606762, 1.02714850549377,
             1.41377744833991, 1.46971420917116, 1.00282382763453, 1.21398240258295,
             1.04693838787649, 1.22122857201799, -0.278142239285745, -1.66163204016444,
             -2.97764891938217, -2.72835223764781, -2.69519338494753, -3.3786609568656,
             -2.27101263544741, -0.128952134210897, 0.096927945318385, -0.229987080991259,
             0.00665070752832922, 0.301143994569227, -0.316074474269264, -0.380465207679694,
             0.302912353732166, -0.065321989262304, -0.607572288802094, 0.392367084365822,
             -0.127209264023268, -0.206946350354197, 0.808881708209294, -0.438317888033248,
             -0.0312989853432367, -1.12867627606029, 0.830779146814718, 0.0641115345868785,
             -0.22363758201639, -0.144068772066544, 0.00459242742023205, -0.30467454720269,
             -0.172618848409753, 0.1321666634488, -0.0398515900419236, 0.0803532071321392,
             -0.0198171373665182, 0.135419065526702, 0.00527843003313361,
             -0.30467454720269, -0.172618848409752, 0.1321666634488, -0.0398515900419236,
             0.0803532071321392, -0.0198171373665182, 0.151431610867278, 0.00148948368736681,
             -0.30467454720269, -0.172618848409753, 0.1321666634488, -0.0398515900419236,
             0.0803532071321391, -0.0198171373665182, -0.213154567577446,
             0.00108763281515294, -0.30467454720269, -0.172618848409752, 0.1321666634488,
             -0.0398515900419236, 0.0803532071321392, -0.0198171373665182,
             -0.0444868686616684, -0.00229461030160591, -0.30467454720269,
             -0.172618848409753, 0.1321666634488, -0.0398515900419236, 0.0803532071321391,
             -0.0198171373665182, 0.0805546426839155, -0.00366057469586981,
             -0.30467454720269, -0.172618848409753, 0.1321666634488, -0.0398515900419236,
             0.0803532071321392, -0.0198171373665182, -0.0243231015211255,
             -0.00312574627200779, -0.30467454720269, -0.172618848409752,
             0.1321666634488, -0.0398515900419236, 0.0803532071321391, -0.0198171373665182,
             0.0530880909433676)

##  model checking
resid.mle <- sigex.resid(psi.mle,mdl,data.ts)[[1]]
resid.mle <- sigex.load(t(Re(resid.mle)),start(data.ts),frequency(data.ts),colnames(data.ts),TRUE)
resid.acf <- acf(resid.mle,lag.max=2*53,plot=FALSE)$acf

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

# HERE


## manage output
#psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
#hess <- fit.mle[[1]]$hessian
#par.mle <- fit.mle[[2]]












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

