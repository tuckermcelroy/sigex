day2date <- function(day,start.date)
{
	# related to date2day, takes a day count past start.date, and 
	#	returns resulting date in month-day-year format.
	#	if day is negative, yields date before start.date

	year.counter <- start.date[3]
	
	day.index <- day + date2day(start.date[1],start.date[2],start.date[3])
	day.counting <- day.index

	if(day.index > 0) {
	while(day.index > 0)
	{
		leap.flag <- 0
		if(((year.counter %% 4 == 0) && (year.counter %% 100 != 0)) || 
			(year.counter %% 400 == 0)){ leap.flag <- 1 }
		day.counting <- day.index
		day.index <- day.index - 365
		if(leap.flag) day.index <- day.index -1
		year.counter <- year.counter + 1
	}
	year.counter <- year.counter - 1
	day.index <- day.counting  }

	if(day.index <= 0) {
	year.counter <- year.counter - 1
 	while(day.index <= 0)
	{
		leap.flag <- 0
		if(((year.counter %% 4 == 0) && (year.counter %% 100 != 0)) || 
			(year.counter %% 400 == 0)){ leap.flag <- 1 }
		day.counting <- day.index
		day.index <- day.index + 365
		if(leap.flag) day.index <- day.index +1
		year.counter <- year.counter - 1
	}
	year.counter <- year.counter + 1  }
	
	if((1 <= day.index) && (day.index <= 31)) month.counter <- 1
	if((32 <= day.index) && (day.index <= (59+leap.flag))) month.counter <- 2
	if(((60+leap.flag) <= day.index) && (day.index <= (90+leap.flag))) month.counter <- 3
	if(((91+leap.flag) <= day.index) && (day.index <= (120+leap.flag))) month.counter <- 4
	if(((121+leap.flag) <= day.index) && (day.index <= (151+leap.flag))) month.counter <- 5
	if(((152+leap.flag) <= day.index) && (day.index <= (181+leap.flag))) month.counter <- 6
	if(((182+leap.flag) <= day.index) && (day.index <= (212+leap.flag))) month.counter <- 7
	if(((213+leap.flag) <= day.index) && (day.index <= (243+leap.flag))) month.counter <- 8
	if(((244+leap.flag) <= day.index) && (day.index <= (273+leap.flag))) month.counter <- 9
	if(((274+leap.flag) <= day.index) && (day.index <= (304+leap.flag))) month.counter <- 10
	if(((305+leap.flag) <= day.index) && (day.index <= (334+leap.flag))) month.counter <- 11
	if(((335+leap.flag) <= day.index) && (day.index <= (365+leap.flag))) month.counter <- 12
	day.counter <- day.index - date2day(month.counter,1,year.counter) + 1

	return(c(month.counter,day.counter,year.counter))
}
