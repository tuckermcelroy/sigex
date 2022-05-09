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


# ---- Signal extraction from fit mdl ----
# bundle
analysis.mle <- sigex.bundle(data.ts, transform, mdl, psi.mle)

## load up the fitted model for signal extraction
dataNA.ts <- analysis.mle[[1]]
mdl <- analysis.mle[[3]]
psi <- analysis.mle[[4]]
param <- sigex.psi2par(psi, mdl, dataNA.ts)
param2 <- with(analysis.mle, sigex.psi2par(psi, model, data))

## get fixed effects
reg.trend <- sigex.fixed(data.ts,mdl,1,param,"Trend")
# reg.nyd <- sigex.fixed(data.ts,mdl,1,param,"NewYearDay")
# reg.mlk <- sigex.fixed(data.ts,mdl,1,param,"MLK")
dataLIN.ts <- data.ts - ts(reg.trend,
                           start=start(data.ts),frequency=period)

## define trend and SA weekly filters
week.period <- 365.25/7
half.len <- floor(week.period/2)
x11.filters <- x11filters(week.period,1)
trend.filter <- x11.filters[[1]]
seas.filter <- x11.filters[[2]]
sa.filter <- x11.filters[[3]]
shift <- (dim(sa.filter)[3]-1)/2

## compute extractions
trend.comp <- sigex.adhocextract(psi,mdl,dataNA.ts,trend.filter,half.len,0,TRUE)
sa.comp <- sigex.adhocextract(psi,mdl,dataNA.ts,sa.filter,shift,0,TRUE)
# AO.errs <- dataLIN.ts[AO.times] - data.casts[[1]]
# sa.comp[[1]][AO.times] <- sa.comp[[1]][AO.times] + AO.errs

## plotting
trendcol <- "tomato"
cyccol <- "orchid"
seascol <- "seagreen"
sacol <- "navyblue"
fade <- 60

#pdf(file="bfs-signal-ba.pdf",height=8,width=10)
plot(data.ts)
sigex.graph(extract = trend.comp,
            reg = reg.trend,
            start.date = start(data.ts),
            period = period,
            series = 1,
            displace = 0,
            color = trendcol,
            fade = fade)
sigex.graph(sa.comp,reg.trend,start(data.ts),
            period,1,0,sacol,fade)
# dev.off()

## spectral diagnostics: seasonal adjustment
sigex.specar(sa.comp[[1]],FALSE,1,period)
dev.off()
