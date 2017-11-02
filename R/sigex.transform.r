sigex.transform <- function(data,transform,aggregate=FALSE)
{
	
	#############################
	#	SigEx transform
	#		by Tucker McElroy
	#
	#	Transforms data according to specified function:
	#		"log", "logistic", "none"
	#		aggregate is Boolean, TRUE indicates to sum series
	#
	##############################3
	
	if(aggregate) data <- rowSums(data)
	if(transform=="log") data <- log(data)
	if(transform=="logistic") data <- log(data) - log(1-data)
	if(transform=="none") data <- data
	
	return(data)
}


