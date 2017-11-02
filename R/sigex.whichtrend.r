sigex.whichtrend <- function(mdl)
{
	#############################
	#	sigex.whichtrend
	#		by Tucker McElroy
	#
	#	Determines which component (index) corresponds to trend
	#
	##########################

numcomp <- length(mdl[[3]])
trendcomp <- NULL
for(i in 1:numcomp)
{
	if(sum(mdl[[3]][[i]])==0) trendcomp <- i
}
return(trendcomp)
}
