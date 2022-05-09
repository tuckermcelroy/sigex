
a <-
sigex.acf(L.par = fit_fixed[[2]][[1]][[1]],
          D.par = fit_fixed[[2]][[2]][[1]],
          mdl = mdl,
          comp = 1,
          mdlPar = fit_fixed[[2]][[3]][[1]],
          delta = delta,
          maxlag = 70,
          freqdom = FALSE)


plot(a, type = 'h')
abline(v = 52, lty = 'dotted')
