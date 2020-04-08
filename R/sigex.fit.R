sigex.fit <- function(x, model_components){

  # Verify that x is mts
  if(!('mts' %in% class(x))) stop('x must be of class mts')

  # dimensionality
  N <- dim(x)[2]
  T <- dim(x)[1]

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
  mdl <- sigex.meaninit(mdl,x,0)

  # Model fitting
  constraint <- NULL
  par.mle <- sigex.default(mdl,x,constraint)
  psi.mle <- sigex.par2psi(par.mle,mdl)

  # bundle for default span
  analysis.mle <- sigex.bundle(data.ts = x,
                               transform = transform,
                               mdl = mdl,
                               psi = psi.mle)

  # Signal Extraction
  data.ts <- analysis.mle[[1]]
  mdl <- analysis.mle[[3]]
  psi <- analysis.mle[[4]]
  param <- sigex.psi2par(psi,mdl,data.ts)

  # get signal filters - OLD
  # signal.trend <- sigex.signal(data.ts,param,mdl,1)
  # signal.seas <- sigex.signal(data.ts,param,mdl,2)
  # signal.sa <- sigex.signal(data.ts,param,mdl,c(1,3))

  # get signal filters
  signal.filters <- list()
  numComp <- length(mdl[[1]])
  for(i in 1:(numComp - 1)){
    signal.filters[[i]] <- sigex.signal(data.ts, param, mdl, i)
  }

  # get extractions - OLD
  # extract.trend <- sigex.extract(data.ts,signal.trend,mdl,param)
  # extract.seas <- sigex.extract(data.ts,signal.seas,mdl,param)
  # extract.sa <- sigex.extract(data.ts,signal.sa,mdl,param)

  # get extractions
  signal.extract <- list()
  for(i in 1:(numComp - 1)){
    signal.extract[[i]] <- sigex.extract(data.ts, signal.filters[[i]], mdl, param)
  }

  # return(list(
  #   signal.trend,
  #   signal.seas,
  #   signal.sa,
  #   extract.trend,
  #   extract.seas,
  #   extract.sa
  # ))

  return(list(
    signal.filters = signal.filters,
    signal.extract = signal.extract
  ))
}
