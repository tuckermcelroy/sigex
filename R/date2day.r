date2day <- function(month,day,year)
{
	# Takes month day year date, and returns day index from 1 to 366

	leap.flag <- 0
	if(((year %% 4 == 0) && (year %% 100 != 0)) || (year %% 400 == 0)){ 
		leap.flag <- 1}
	
	day.count <- 0
	if(month==1) { day.count <- day } else {
	for(month.count in 1:(month-1))
	{
		if(month.count==1) day.add <- 31
		if(month.count==2){ 
			day.add <- 28
			if(leap.flag) day.add <- 29
		}
		if(month.count==3) day.add <- 31
		if(month.count==4) day.add <- 30
		if(month.count==5) day.add <- 31
		if(month.count==6) day.add <- 30
		if(month.count==7) day.add <- 31
		if(month.count==8) day.add <- 31
		if(month.count==9) day.add <- 30
		if(month.count==10) day.add <- 31
		if(month.count==11) day.add <- 30
		day.count <- day.count + day.add
	}
	day.count <- day.count + day
	}

	return(day.count)	
}
