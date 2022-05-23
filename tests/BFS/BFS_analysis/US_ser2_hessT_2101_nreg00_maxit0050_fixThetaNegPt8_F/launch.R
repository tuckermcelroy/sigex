
rm(list=ls())
devtools::load_all(".")

# Returns path to this file. This is where results will be saved
modelPath <- getSrcDirectory(function(x) {x})

# Data processing
begin <- c(2006,1)
end <- c(2020,27)
period <- 52
transform <- "log"
aggregate <- FALSE
subseries <- 2
range <- NULL

# Delta operator
delta.s <- ubgenerator(365.25/7,NULL,1000,1)
delta.full <- polymult(c(1,-1), delta.s)
delta <- delta.full

# Mdl
class <- "sarma"
order <- c(0, 1, 0, 1, 365.25/7)
constraint <- matrix(c(-2.197, 0, 0, 1, 0), nrow = 1) # add constraint
regVars <- NULL

# Should fit.mle be run?
runMLE <- TRUEm

source(file.path(modelPath, "run.R"))


