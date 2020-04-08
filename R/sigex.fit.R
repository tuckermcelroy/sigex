sigex.fit <- function(x,
                      model_components = "trend",
                      return_extractions = NULL){

  # Verify that x is mts
  if(!('mts' %in% class(x))) stop('x must be of class mts')

  # dimensionality
  N <- dim(x)[2]
  T <- dim(x)[1]
  period <- frequency(x)

  # setup mdl if character vector is specified by user
  if(is.character(model_components))
  {
    mdl <- NULL
    if("trend" %in% model_components){
      mdl <- sigex.add(mdl, seq(1,N), "arma", c(0,0), 0, "trend", c(1,-2,1))
    }
    if("seasonal" %in% model_components){
      mdl <- sigex.add(mdl, seq(1,N), "arma", c(0,0), 0, "seasonal", rep(1, period))
    }
    mdl <- sigex.add(mdl, seq(1,N), "arma", c(0,0), 0, "irregular", 1)
    # regressors
    mdl <- sigex.meaninit(mdl, x, 0)
  }

  # Define what extraction components to return
  # Null means user did not supply values
  if(is.null(return_extractions)){
    return_extractions <- list()
    for(i in 1:(length(mdl[[1]])-1)) return_extractions[[i]] = i
  }


  # Model fitting
  constraint <- NULL
  par.mle <- sigex.default(mdl,x,constraint)
  psi.mle <- sigex.par2psi(par.mle,mdl)

  ## run fitting: commented out, this takes a while
  # fit.mle <- sigex.mlefit(x, par.mle, constraint, mdl, "bfgs", debug=FALSE)

  ## manage output
  # psi.mle <- sigex.eta2psi(fit.mle[[1]]$par, constraint)
  # hess <- fit.mle[[1]]$hessian
  # par.mle <- fit.mle[[2]]

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

  # get signal filters
  signal.filters <- list()
  numComp <- length(return_extractions)
  for(i in 1:(numComp)){
    signal.filters[[i]] <- sigex.signal(data.ts, param, mdl, return_extractions[[i]])
  }

  # get extractions
  signal.extract <- list()
  for(i in 1:(numComp)){
    signal.extract[[i]] <- sigex.extract(data.ts, signal.filters[[i]], mdl, param)
  }

  return(list(
    signal.extract = signal.extract
  ))
}
