
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
s.period <- 365.25/7
rho.s <- 1
s.div <- floor(s.period)
s.frac <- s.period - s.div
sar.op <- c(1, rep(0, s.div - 1), (s.frac - 1) * rho.s, -1 * s.frac * rho.s)
delta <- sar.op

# Mdl
class <- "sarmaf"
order <- c(0,1,0,1,365.25/7)
# Constraint making SMA = -.8
constraint <- constraint <- matrix(c(-2.197, 0, 0, 1, 0), nrow = 1)
regVars <- NULL

# Should fit.mle be run?
runMLE <- TRUE

source(file.path(modelPath, "run.R"))


