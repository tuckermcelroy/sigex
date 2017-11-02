sigex.zetalen <- function(mdlType)
{

	#################################
	#   sigex.zetalen
	#	by Tucker McElroy	
	#
	#	Takes   mdlType and produces length of corresponding zeta
	#		psi = [xi,zeta,beta]
	#			xi ~ all pre-parameters for covariance matrices
	#			zeta ~ all pre-parameters for t.s. models
	#			beta ~ all regression parameters
	#
	#################################


	if(mdlType == "wn") { zetalen <- 0 }
	if(mdlType == "canonWN") { zetalen <- 0 }
	if(mdlType == "AR1") { zetalen <- 1 }
	if(mdlType %in% c("cycleBW1","cycleBW2","cycleBW3","cycleBW4","cycleBW5",
		"cycleBW6","cycleBW7","cycleBW8","cycleBW9","cycleBW10"))
		{ zetalen <- 2 }
	if(mdlType %in% c("canonCycleBW1","canonCycleBW2","canonCycleBW3",
		"canonCycleBW4","canonCycleBW5","canonCycleBW6","canonCycleBW7",
		"canonCycleBW8","canonCycleBW9","canonCycleBW10"))
		{ zetalen <- 2 }
	if(mdlType %in% c("cycleBAL1","cycleBAL2","cycleBAL3","cycleBAL4","cycleBAL5",
		"cycleBAL6","cycleBAL7","cycleBAL8","cycleBAL9","cycleBAL10"))
		{ zetalen <- 2 }
	if(mdlType %in% c("canonCycleBAL1","canonCycleBAL2","canonCycleBAL3",
		"canonCycleBAL4","canonCycleBAL5","canonCycleBAL6","canonCycleBAL7",
		"canonCycleBAL8","canonCycleBAL9","canonCycleBAL10"))
		{ zetalen <- 2 }
#	if(mdlType == "VAR1") { zetalen <- N^2 }
	if(mdlType == "ARMA22") { zetalen <- 4 }
	if(mdlType == "damped") { zetalen <- 1 + as.integer(N*(N+1)/2) }

	return(zetalen)
}
