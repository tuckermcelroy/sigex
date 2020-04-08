#' The purpose of this script is to mimic the workflow we expect from a user
#' of our software.
#' The current steps are to:
#'     1. load data
#'     2. prep data to mts object
#'     3. envoke sigex
#'     4. look at diagnositc plots


# ---- 1. Load data ----

head(starts)
class(starts)

# ---- 2. prep as MTS object ----

start.date <- c(1964,1)
period <- 12
x <- ts(starts[, c(1, 3)], start = start.date, frequency = period)


# ---- 3. Envoke sigex ----

fit <- sigex.fit(x = x, model_components = c("trend", "seasonal"))
fit2 <- sigex.fit(x = x, model_components = c("trend"))



# ---- 4. Diagnostics


# get signal filters from output
signal.trend <- fit[[1]]
signal.seas <- fit[[2]]
signal.sa <- fit[[3]]

# get extractions from output
extract.trend <- fit[[4]]
extract.seas <- fit[[5]]
extract.sa <- fit[[6]]

trendcol <- "tomato"
cyccol <- "orchid"
seascol <- "seagreen"
sacol <- "navyblue"
fade <- 60
par(mfrow=c(1,2))
for(i in 1:2)
{
  plot(x[,i],xlab="Year",ylab="",ylim=c(min(x[,i])-20,max(x[,i])),lwd=1)
  sigex.graph(extract.sa,reg.trend,begin.date,period,i,0,sacol,fade)
  sigex.graph(extract.trend,reg.trend,begin.date,period,i,0,trendcol,fade)
  sigex.graph(extract.seas,NULL,begin.date,period,i,-5,seascol,fade)
}
