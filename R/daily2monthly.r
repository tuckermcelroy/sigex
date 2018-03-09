daily2monthly <- function(daily.series,start.date,dayofweek)
{
	# Takes a daily.series from start.date to end.date, and computes
	#	a monthly time series by aggregating all daily values of
	#	a given dayofweek type that fall in a given month
	#
  #
	# uses day2date, date2day, day2week

	my.dates <- t(apply(as.matrix(seq(0,length(daily.series)-1,1)),1,day2date,start.date))
	my.days <- apply(my.dates,1,day2week)
	dated.daily <- cbind(daily.series,my.days,my.dates)
	sub.series <- dated.daily[dated.daily[,2]==dayofweek,]
	monthly.series <- unique(sub.series[,c(3,5)])
	monthly.sums <- NULL
	for(i in 1:dim(monthly.series)[1])
	{
		indices <- (sub.series[,3]==monthly.series[i,1]) +
			(sub.series[,5]==monthly.series[i,2])
		new.sum <- sum(sub.series[indices==2,1])
		monthly.sums <- c(monthly.sums,new.sum)
	}
	monthly.series <- cbind(monthly.sums,monthly.series)

	return(monthly.series)
}


