sigex_fit <- function(x, model_components){

  # Verify that x is mts
  if(!('mts' %in% class(x))) stop('x must be of class mts')

  # Sigex.prep defaults
  transform <- "none"
  aggregate <- FALSE
  subseries <- 1:dim(x)[2]
  begin.date <- start(x)
  end.date <- end(x)
  range <- list(begin.date,end.date)

  # prep data
  data.ts <- sigex.prep(data.ts = dataALL.ts,
                        transform = transform,
                        aggregate = aggregate,
                        subseries = subseries,
                        range = range,
                        plot = FALSE)

  # dimensionality
  N <- dim(data.ts)[2]
  T <- dim(data.ts)[1]

  # setup mdl
  mdl <- NULL
  if("trend" %in% model_components){
    mdl <- sigex.add(mdl, seq(1,N), "arma", c(0,0), 0, "trend", c(1,-2,1))
  }
  if("seasonal" %in% model_components){
    mdl <- sigex.add(mdl, seq(1,N), "arma", c(0,0), 0, "seasonal", rep(1, 12))
  }
  mdl <- sigex.add(mdl, seq(1,N), "arma", c(0,0), 0, "irregular", 1)
  # regressors
  mdl <- sigex.meaninit(mdl,data.ts,0)

  # Model fitting
  constraint <- NULL
  par.mle <- sigex.default(mdl,data.ts,constraint)
  psi.mle <- sigex.par2psi(par.mle,mdl)

  # bundle for default span
  analysis.mle <- sigex.bundle(data.ts = data.ts,
                               transform = transform,
                               mdl = mdl,
                               psi = psi.mle)

  # Signal Extraction
  data.ts <- analysis.mle[[1]]
  mdl <- analysis.mle[[3]]
  psi <- analysis.mle[[4]]
  param <- sigex.psi2par(psi,mdl,data.ts)

  # get signal filters
  signal.trend <- sigex.signal(data.ts,param,mdl,1)
  signal.seas <- sigex.signal(data.ts,param,mdl,2)
  signal.sa <- sigex.signal(data.ts,param,mdl,c(1,3))

  # get extractions
  extract.trend <- sigex.extract(data.ts,signal.trend,mdl,param)
  extract.seas <- sigex.extract(data.ts,signal.seas,mdl,param)
  extract.sa <- sigex.extract(data.ts,signal.sa,mdl,param)

}
