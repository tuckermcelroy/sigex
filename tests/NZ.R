#########################################
###  Script for Daily Immigration Data
#########################################

## wipe
rm(list=ls())

library(devtools)

setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex")
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
psi.mle <- c(0.42339322295077, 0.251689511813033, 0.314297263134772, 0.292884125504594,
             0.257913164651892, 0.148602188503815, 0.712452057692297, 0.419275740694886,
             0.297429319854366, 0.207634209891519, 0.186822165568141, 0.472031179535268,
             0.45988914772637, 0.413178135931064, 0.408962630726071, 0.624491760357869,
             0.562091978394646, 0.357407169620627, 0.556706146637411, 0.476584809876487,
             0.583332175453021, -3.19866501190753, -3.85775135156783, -4.11467196059505,
             -4.12374509082116, -4.08224616367834, -4.2642872390669, -3.91329379347115,
             1.103090080827, 0.516527181388701, 0.497603007560039, -0.117397825595269,
             0.186268760242241, -0.538319607809195, 0.451918434348235, 0.49598863791487,
             0.0957095518741189, 0.442821257794782, 0.530601282475884, 0.394910156898995,
             0.0991632309965594, 0.446239448570845, 0.0610477578408184, 1.73426975902616,
             -0.385486868923238, 0.628493278474925, 1.99713719328107, -0.492822137460735,
             1.32399206628853, -1.70593014011961, -0.911006847269805, -1.12499863618325,
             -2.7454219597856, -2.8810577831584, -2.27085248435103, -3.07945084889349,
             -1.12534658785179, 2.40256784511735, -0.544987506790035, -1.87317007554894,
             1.94366938375088, 0.300541874228645, -1.24530804950616, 1.41447655896811,
             1.80549340435235, -2.89712955567044, 2.15116699963435, -0.0402193973964761,
             2.88009047625405, 1.61030864683225, -1.59176431992338, 2.23767655819279,
             1.11396524023326, 0.29378384031964, 2.05515869450879, 0.431432852547541,
             2.33761218055479, 0.788314781946021, 0.469604631936633, 0.274031979259612,
             0.129919526414844, 0.118649752196389, 0.238409353664854, 0.974974354108441,
             0.455929002064592, 0.301534046312906, 0.240394857099853, -0.0129878570036447,
             0.886549609460141, 0.0987115319322978, 0.130721944514428, 0.679177199349905,
             1.19306232513805, 1.2198314526085, 0.892904754774129, 1.0104564572229,
             0.710051415679332, 0.801417772593156, -0.567073640100783, -1.76719424658498,
             -3.0285806739591, -2.65767609836693, -2.76534313228085, -3.31951953854884,
             -2.47233158131481, -0.101256371548154, 0.130020620914774, -0.200241833762783,
             0.0309157519905209, 0.27890358279936, -0.328181439944388, -0.417844236422465,
             0.267461634096443, -0.0661338967900519, -0.530123153322867, 0.403562948918957,
             -0.0913230565664687, -0.163930981177709, 0.682399968554364, -0.425258343078195,
             -0.117124250962673, -1.05842180324038, 0.756027603009854, -0.012999473505258,
             -0.193099353094616, -0.219810519935073, 0.00513551502962987,
             -0.340205300755847, -0.187977048051066, 0.136650273760756, 0.00386244713360125,
             0.0527682412851075, -0.00711720278859726, 0.137999609126442,
             -0.00115686493669871, -0.340205300755846, -0.187977048051066,
             0.136650273760756, 0.00386244713360125, 0.0527682412851075, -0.00711720278859722,
             0.136551158042287, -0.00309756937625004, -0.340205300755847,
             -0.187977048051067, 0.136650273760756, 0.00386244713360124, 0.0527682412851075,
             -0.00711720278859731, -0.185872432356347, 0.00263140416117883,
             -0.340205300755847, -0.187977048051066, 0.136650273760756, 0.00386244713360123,
             0.0527682412851075, -0.00711720278859726, -0.0361260089399556,
             0.000374322462658291, -0.340205300755847, -0.187977048051066,
             0.136650273760756, 0.00386244713360129, 0.0527682412851075, -0.00711720278859725,
             0.0939340431145078, 0.000592795432065807, -0.340205300755847,
             -0.187977048051066, 0.136650273760756, 0.00386244713360126, 0.0527682412851075,
             -0.00711720278859725, -0.0619191509115667, -0.0011743954320999,
             -0.340205300755847, -0.187977048051066, 0.136650273760756, 0.00386244713360125,
             0.0527682412851075, -0.00711720278859726, 0.041665222148067)


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

