
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
s.period <- 365.25/7
rho.s <- 1
s.div <- floor(s.period)
s.frac <- s.period - s.div
sar.op <- c(1, rep(0, s.div - 1), (s.frac - 1) * rho.s, -1 * s.frac * rho.s)
delta <- sar.op

# Mdl
class <- "sarmaf"
order <- c(2, 0, 0, 1, 365.25/7)
constraint <- NULL
regVars <- c("nyd.reg", "mlk.reg", "labor.reg")

# Should fit.mle be run?
runMLE <- FALSE

source(file.path(modelPath, "run.R"))


