
acf(diff(data.ts), lag.max = 4 * 53)$acf
pacf(diff(data.ts), lag.max = 4 * 53)


spec.ar(diff(data.ts), order = 200)

spec.ar(data.ts)

spec.ar(AirPassengers)


TSA::periodogram(diff(data.ts))

abline(v = (1:52)/52, lty = 'dotted', col = 'red')
abline(v = (1:53)/53, lty = 'dotted', col = 'blue')
p <- 365.25/7
abline(v = (1:p)/p, lty = 'dotted', col = 'blue')



plot(diff(data.ts))

abline(v = 2006:2020, lty = 'dotted', col = 'red')
