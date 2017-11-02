sigex.graph <- function(extract,reg=NULL,start,period,series,displace,color,fade)
{
	
	###############################
	#   sigex.graph
	#	by Tucker McElroy	
	#
	#	Adds signal extraction estimates, with shaded bands, 
	#		to an existing time series plot
	#		extract is the output of sigex.extract
	#		reg is TxN matrix of fixed effects, to be added to extraction,
	#			only column "series" is utilized.
	#			note: for multiple reg effects, add these all up before hand
	#		start is starting year of the time series
	#		period is number of seasons per year
	#		series is which of the N series to plot
	#		displace gives a vertical shift to the plot
	#		color given as a number in the range of colors()
	#		fade gives shading proportion for uncertainty
	#
	#################################

	T <- dim(extract[[1]])[1]
	rgbCol <- col2rgb(color)
	hexCol <- rgb(rgbCol[1,1],rgbCol[2,1],rgbCol[3,1],maxColorValue=255)
	if(length(reg)==0) { reg.shift <- matrix(0,nrow=T,ncol=1) } else {
		reg.shift <- as.matrix(reg[,series]) }

	time <- start[1] + (start[2]-1)/period + (seq(1,T)-1)/period
# which is right????
#	time <- start[1] + start[2]/period + (seq(1,T)-1)/period

	lines(ts(extract[[1]][,series]+displace+reg.shift,
		start=start,frequency=period),col=hexCol,lwd=2)
	polygon(c(time,rev(time)),c(extract[[3]][,series]+rep(displace,T)+reg.shift,
		rev(extract[[2]][,series]+rep(displace,T)+reg.shift)),
		col=paste(hexCol,fade,sep=""),border=NA)
}	
