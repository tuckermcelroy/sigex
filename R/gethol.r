#' generates holiday regressors from holiday dates
#'
#' @param hol.dates list of dates in month-day-year format;
#'			this is a matrix with one row for each holiday date,
#'			and 3 columns for month, day in month, and year
#' @param hol.fore number of days before holiday, where effect begins
#' @param hol.aft number of days after holiday, where effect ends
#' @param start.date start date of time series in month-day-year format,
#'  		a 3-element vector
#' @param end.date end date of time series in month-day-year format,
#'			a 3-element vector
#'
#' @return subset.hol.reg: regressors corresponding to the holiday effect,
#'			with an activity window determined by hol.fore and hol.aft,
#'			and at dates between start.date and end.date
#' @export
#'

gethol <- function(hol.dates,hol.fore,hol.aft,start.date,end.date)
{

	##########################################################################
	#
	#	gethol
	# 	    Copyright (C) 2018  Tucker McElroy
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

	################# Documentation #####################################
	#
	#	Purpose: generates holiday regressors from holiday dates
	#
	#	Inputs:
	#		hol.dates: list of dates in month-day-year format;
	#			this is a matrix with one row for each holiday date,
	#			and 3 columns for month, day in month, and year
	#		hol.fore: number of days before holiday, where effect begins
	#		hol.aft: number of days after holiday, where effect ends
	#		start.date: start date of time series in month-day-year format,
	#			a 3-element vector
	#		end.date: end date of time series in month-day-year format,
	#			a 3-element vector
	#	Outputs:
	#		subset.hol.reg: regressors corresponding to the holiday effect,
	#			with an activity window determined by hol.fore and hol.aft,
	#			and at dates between start.date and end.date
	#	Requires: date2day
	#	Notes:  this version presumes start.date and end.date are in separate years
	#
	#####################################################################

	hol.day <- NULL
	nHol <- dim(hol.dates)[1]
	for(year in 1:nHol)
	{
		hol.day <- c(hol.day,date2day(hol.dates[year,1],
			hol.dates[year,2],hol.dates[year,3]))
	}

	hol.reg <- NULL
	for(year in hol.dates[1,3]:hol.dates[nHol,3])
	{
		hol.index <- which(hol.dates[,3]==year)
		hol.new <- rep(0,date2day(12,31,year))
		hol.new[hol.day[hol.index]] <- 1
		hol.reg <- c(hol.reg,hol.new)
	}
	hol.index <- which(hol.reg==1)

	ny <- 0
	leap.index <- date2day(2,28,1600) + 1	# equals 60
	hol.means <- rep(0.0,366)
	for(i in 1:length(hol.index))
	{
		year.len <- date2day(12,31,hol.dates[i,3])
		if((hol.index[i]-hol.fore > 0) && (hol.index[i]+hol.aft <= length(hol.reg)))
		{ hol.reg[(hol.index[i]-hol.fore):(hol.index[i]+hol.aft)] <- 1.0 }
		hol.year.reg <- hol.reg[(ny+1):(ny+year.len)]
		if(year.len==365) {
			hol.means[-leap.index] <- hol.means[-leap.index] + hol.year.reg } else
		{	hol.means <- hol.means + hol.year.reg }
		ny <- ny + year.len
	}
	hol.means <- hol.means/length(hol.index)

	start.day <- date2day(start.date[1],start.date[2],start.date[3])
	end.day <- date2day(end.date[1],end.date[2],end.date[3])

	subset.hol.reg <- NULL
	ny <- 0
	for(year in hol.dates[1,3]:end.date[3])
	{
		if(year >= start.date[3]) {
			i1 <- date2day(1,1,year) + ny
			i2 <- date2day(12,31,year) + ny
			hol.new <- hol.reg[i1:i2]
			if(date2day(12,31,year)==365) {
				hol.new <- hol.new - hol.means[-leap.index]
			} else { hol.new <- hol.new - hol.means }
			if(year==start.date[3]) { hol.new <- hol.new[start.day:length(hol.new)] }
			if(year==end.date[3]) { hol.new <- hol.new[1:end.day] }
			subset.hol.reg <- c(subset.hol.reg,hol.new)
		}
		ny <- ny + date2day(12,31,year)
	}

	return(subset.hol.reg)
}

