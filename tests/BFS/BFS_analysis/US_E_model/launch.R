
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
subseries <- 1
range <- NULL

# Delta operator
delta.s <- ubgenerator(365.25/7,NULL,1000,1)
delta.full <- polymult(c(1,-1), delta.s)
delta <- polymult(c(1, -1), delta.full) # add additional 1st diff

# Mdl
class <- "sarma"
order <- c(1, 2, 0, 1, 365.25/7)
constraint <- NULL
regVars <- NULL

# Should fit.mle be run?
runMLE <- FALSE

source(file.path(modelPath, "run.R"))


