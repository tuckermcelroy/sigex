sigex.graph <- function(extract,reg=NULL,start.date,period,series,displace,color,fade)
{
	
	##########################################################################
	#
	#	sigex.graph
	# 	    Copyright (C) 2017  Tucker McElroy
	#
	#    This program is free software: you can redistribute it and/or modify
	#    it under the terms of the GNU General Public License as published by
	#    the Free Software Foundation, either version 3 of the License, or
	#    (at your option) any later version.
	#
	#    This program is distributed in the hope that it will be useful,
	#    but WITHOUT ANY WARRANTY; without even the implied warranty of
	#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	#    GNU General Public License for more details.
	#
	#    You should have received a copy of the GNU General Public License
	#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
	#
	############################################################################

	################# Documentation ############################################
	#
	#	Purpose: adds signal extraction estimates, with shaded bands, 
	#		to an existing time series plot
	#	Background:	
	#		A sigex model consists of process x = sum y, for 
	#		stochastic components y.  Each component process y_t
	#		is either stationary or is reduced to stationarity by
	#		application of a differencing polynomial delta(B), i.e.
	#			w_t = delta(B) y_t   is stationary.
	#		We have a model for each w_t process, and can compute its
	#		autocovariance function (acf), and denote its autocovariance
	#		generating function (acgf) via gamma_w (B).
	#		The signal extraction filter for y_t is determined from
	#		this acgf and delta.
	#	Inputs:
	#		extract: T x N matrix of the signal estimates, e.g. output of sigex.extract
	#		reg: a T x N matrix of fixed effects, to be added to extraction.
	#			Note: only the column of reg corresponding to "series" is utilized.
	#			(To incorporate multiple reg effects, add these all up before hand.)
	#		start.date: date of first time obersvation; the
	#			 format is c(year,season)
	#		period: number of seasons per year
	#		series: index of the particular series under consideration
	#		displace: gives a vertical shift to the plot
	#		color: given as a number in the range of colors()
	#		fade: gives shading proportion for uncertainty (NULL if none)
	#	Outputs: none
	#
	####################################################################
 
	T <- dim(extract[[1]])[1]
	rgbCol <- col2rgb(color)
	hexCol <- rgb(rgbCol[1,1],rgbCol[2,1],rgbCol[3,1],maxColorValue=255)
	if(length(reg)==0) { reg.shift <- matrix(0,nrow=T,ncol=1) } else {
		reg.shift <- as.matrix(reg[,series]) }

	time <- start.date[1] + (start.date[2]-1)/period + (seq(1,T)-1)/period
	lines(ts(extract[[1]][,series]+displace+reg.shift,
		start=start.date,frequency=period),col=hexCol,lwd=1)
	if(length(fade)>0) {
	polygon(c(time,rev(time)),c(extract[[3]][,series]+rep(displace,T)+reg.shift,
		rev(extract[[2]][,series]+rep(displace,T)+reg.shift)),
		col=paste(hexCol,fade,sep=""),border=NA) }
}	
