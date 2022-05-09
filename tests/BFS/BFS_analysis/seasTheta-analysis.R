# Purpose: Want to understand relationship between seasonal-theta and
#   the diagnostics in a sigex fit.


# ---- Do analysis for regular SARIMA ----

x <- AirPassengers

fit <- arima(x, order = c(0, 1, 1), seasonal = c(0, 1, 1))
plot(resid(fit))
acf(resid(fit))
fit$coef



fit2 <- arima(log(x), order = c(0, 1, 1), seasonal = c(0, 1, 1))


plot(AirPassengers)
plot(log(AirPassengers))



dx <- diff(diff((AirPassengers)), 12)
dlogx <- diff(diff(log(AirPassengers)), 12)
plot(dx)
plot(dlogx)


acf(dx)



fit3 <- arima(log(x), order = c(0, 1, 1), seasonal = c(0, 1, 1), fixed = c(NA, -.8))
acf(resid(fit3), plot = FALSE)$acf[13]
fit3$coef


seasTheta_seq <- seq(-.9, 0, .05)
N <- length(seasTheta_seq)
lag12acf <- numeric(N)
thetaEst <- numeric(N)

for(i in seq_along(seasTheta_seq)){
  seasTheta <- seasTheta_seq[i]

  fit3 <- arima(log(x), order = c(0, 1, 1),
                seasonal = c(0, 1, 1),
                fixed = c(NA, seasTheta))
  lag12acf[i] <- acf(resid(fit3), plot = FALSE)$acf[13]
  thetaEst[i] <- fit3$coef[1]
}

plot(lag12acf ~ seasTheta_seq)
abline(h = 1.96*1/sqrt(length(AirPassengers)))
abline(h = -1.96*1/sqrt(length(AirPassengers)))
abline(h = 0 , lty = 'dotted')




# ---- replicate study in sigex ----
# Data processing
begin <- c(2006,1)
end <- c(2020,27)
period <- 52
transform <- "none"
aggregate <- FALSE
subseries <- 1
range <- NULL

fakeData <- matrix(rnorm(757*4), 757, 4)

dataALL.ts <- sigex.load(fakeData,
                         begin,
                         period,
                         c("bfs-ba","bfs-hba","bfs-wba","bfs-cba"),
                         FALSE)


data.ts <- sigex.prep(dataALL.ts,
                      transform,
                      aggregate,
                      subseries,
                      range,
                      TRUE)

N <- ncol(data.ts)
T <- nrow(data.ts)

# simple function to write delta and class into global.env
defDelta <- function(x = c("JD", "US")){
  if(x == "JD"){
    # Delta operator (JD+)
    s.period <- 365.25/7
    rho.s <- 1
    s.div <- floor(s.period)
    s.frac <- s.period - s.div
    sar.op <- c(1, rep(0, s.div - 1), (s.frac - 1) * rho.s, -1 * s.frac * rho.s)
    delta <<- sar.op
    CLASS <<- "sarmaf"
    print("delta = JD+ diff operator, and CLASS = \"sarmaf\"")
  } else if (x == "US") {
    # Delta operator (US)
    delta.s <- ubgenerator(365.25/7,NULL,1000,1)
    delta.full <- polymult(c(1,-1), delta.s)
    delta <<- delta.full
    CLASS <<- "sarma"
    print("delta = our diff operator, and CLASS = \"sarma\"")
  } else {
    warning("Not a valid input for defDelta")
  }
}

# model construction
mdl <- NULL
defDelta("US")
mdl <- sigex.add(mdl = mdl,
                 vrank = 1:N,
                 class = CLASS,
                 order = c(0, 1, 0, 1, 365.25/7),
                 bounds = NULL,
                 name = "process",
                 delta = delta)
mdl <- sigex.meaninit(mdl,data.ts,0)
par_default <- sigex.default(mdl, data.ts, constraint = NULL)

# Fit
# mle fit
fit <- sigex.mlefit(data.ts = data.ts,
                    param = par_default,
                    constraint = NULL,
                    mdl = mdl,
                    method = "bfgs",
                    debug = TRUE)
# Extract par from fit
psi.mle <- fit[[1]]$par
par.mle <- sigex.psi2par(psi.mle, mdl, data.ts)

# Find residuals of fit and store the lag 50-54 ACF values of residuals
resid.mle <- sigex.resid(psi.mle, mdl, data.ts)[[1]]
resid.mle <- sigex.load(t(resid.mle), start(data.ts), frequency(data.ts), colnames(data.ts), FALSE)
resid.acf <- acf(resid.mle, lag.max = 4 * 53, plot = TRUE)$acf
lag52acf[i, ] <- resid.acf[51:55]

# Sequence of seasonal thetas to loop over
n_seq <- 10
seasTheta_seq <- seq(0, .95, length.out = n_seq)

# storage for acf lags 50-54 at each iteration of loop
lag52acf <- matrix(nrow = n_seq, ncol = 5)
colnames(lag52acf) <- paste0('lag', 50:54)

# storage for each fit param list at each iteration of loop
ALL_fit_fixed <- vector(mode = 'list', length = n_seq)

# Main loop
for(i in seq_along(seasTheta_seq)){

  seasTheta <- seasTheta_seq[i]

  print(i)
  print(seasTheta)

  # Find value of psi that maps to desired seasTheta
  par_temp <- par_default
  par_temp[[3]][[1]][2] <- seasTheta
  psi_temp <- sigex.par2psi(par_temp, mdl)
  psi_fixed <- psi_temp[3]

  # Set constraint to correct value of psi
  constraint <- matrix(c(psi_fixed, 0, 0, 1, 0), nrow = 1)

  # initialized par with correct constraint
  par_seas <- sigex.default(mdl, data.ts, constraint)

  # mle fit
  fit_fixed <- sigex.mlefit(data.ts = data.ts,
                            param = par_seas,
                            constraint = constraint,
                            mdl = mdl,
                            method = "bfgs",
                            debug=FALSE)

  # Store output of mlefit
  ALL_fit_fixed[[i]] <- fit_fixed

  # Print and beep to alert that 1 iteration has finished
  print("finished a fit")
  beepr::beep(2)

  # Extract par from fit (extra complexity due to constraint)
  eta.mle <- fit_fixed[[1]]$par
  psi.mle <- sigex.eta2psi(eta = eta.mle, constraint = constraint)
  par.mle <- sigex.psi2par(psi.mle, mdl, data.ts)

  # Find residuals of fit and store the lag 50-54 ACF values of residuals
  resid.mle <- sigex.resid(psi.mle, mdl, data.ts)[[1]]
  resid.mle <- sigex.load(t(resid.mle), start(data.ts), frequency(data.ts), colnames(data.ts), FALSE)
  resid.acf <- acf(resid.mle, lag.max = 4 * 53, plot = TRUE)$acf
  lag52acf[i, ] <- resid.acf[51:55]

}

# save(lag52acf, file = "~/GitHub/sigex/tests/BFS/BFS_analysis/Ser2_JD_lag52acf.RData")
load(file = "~/GitHub/sigex/tests/BFS/BFS_analysis/Ser1_US_lag52acf.RData")

# Create empty plot ACF plot
n <- nrow(data.ts)
par(mfrow = c(1, 1))
plot(range(0, 1.1), range(lag52acf), type = "n", xlab = 'seasonal Theta', ylab = 'ACF')
abline(h = c(1, -1)*1.96*1/sqrt(n)) # bounds for WN
abline(h = 0 , lty = 'dotted')

# Get correct x-axis values
n_seq <- nrow(lag52acf)
seasTheta_seq <- seq(0, .95, length.out = n_seq)

# Make line for each lag of ACF 50, 51, 52, 53, 54
legend("topright", legend = colnames(lag52acf), pch = 19, col = 1:5)
for(i in 1:5){
  points(lag52acf[, i] ~ seasTheta_seq, col = i, pch = 19, type = 'b')
}





































































