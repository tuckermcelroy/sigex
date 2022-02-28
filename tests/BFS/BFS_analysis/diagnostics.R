print(modelPath)

# ---- Check convergence of optim
fit.mle[[1]]$convergence == 0

# ---- * manage output ----
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par, constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]

# Likelihood value at optim value
sigex.lik(psi = psi.mle, mdl = mdl, data.ts = data.ts)

# Print first seasonal parameter afer non-sesaonal AR/MA params
# Seasonal Theta or seasonal Phi
PAR <- sigex.psi2par(psi = fit.mle[[1]]$par, mdl = mdl, data.ts = data.ts)
seasTheta <- PAR[[3]][[1]][sum(order[1:2]) + 1]
seasTheta <- round(seasTheta, 3)

## check on standard errors and get t statistics
eigen(hess)$values

# residual analysis
resid.mle <- sigex.resid(psi.mle, mdl, data.ts)[[1]]
resid.mle <- sigex.load(t(resid.mle), start(data.ts), frequency(data.ts), colnames(data.ts), TRUE)
resid.acf <- acf(resid.mle, lag.max = 4 * 53, plot = TRUE)$acf

# t statistics for parameters
sigex.tstats(mdl,psi.mle,hess,constraint)

## model checking
sigex.portmanteau(resid.mle, 4*period, length(psi.mle))
sigex.gausscheck(resid.mle)
