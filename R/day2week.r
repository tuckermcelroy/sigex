day2week <- function(date)
{
	# Takes month day year format date, and returns day of week:
	#	1 is Sunday, 2 is Monday, etc.
	#  Calibrated to Jan 1, 1600 is a Saturday

	month <- date[1]
	day <- date[2]
	year <- date[3]
	total.days <- 0
	if(year > 1600) {
	for(years in 1600:(year-1))
	{
		total.days <- total.days + date2day(12,31,years)
	} }
	total.days <- total.days + date2day(month,day,year) -1
	dayofweek <- total.days %% 7
	if(dayofweek==0) dayofweek <- 7
	return(dayofweek)
}

