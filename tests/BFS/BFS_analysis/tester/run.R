# This script fits sigex model assuming the following are loaded into the
#    global environment:
# begin
# end
# period
# transform
# aggredgate
# subseries
# range
# delta
# class
# order
# constraint
# regVars [a character vector of regressors]


# * ts object and plot ----
dataALL.ts <- sigex.load(bfs[,3:6],
                         begin,
                         period,
                         c("bfs-ba","bfs-hba","bfs-wba","bfs-cba"),
                         FALSE)

data.ts <- sigex.prep(dataALL.ts,
                      transform,
                      aggregate,
                      subseries,
                      range,
                      FALSE)


# #### Model Declaration ############################################

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]

source(file = "tests/BFS/BFS_analysis/regressors.R")


# model construction
mdl <- NULL
mdl <- sigex.add(mdl    = mdl,
                 vrank  = seq(1,N),
                 class  = class,
                 order  = order,
                 bounds = NULL,
                 name   = "process",
                 delta  = delta )
mdl <- sigex.meaninit(mdl,data.ts,0)

# Add regressors from named character vector regVars
for(v in regVars){
  varName <- strsplit(v, split = "\\.")[[1]][1]
  mdl <- sigex.reg(mdl = mdl,
                   series = 1,
                   reg = ts(as.matrix(get(v)),
                            start = start(get(v)),
                            frequency = period,
                            names = varName)
  )
}

# Set default params
par.default <- sigex.default(mdl, data.ts, constraint)
psi.default <- sigex.par2psi(par.default, mdl)

if(runMLE){
  ## run fitting:
  message_file <- file(file.path(modelPath, 'mle_log.txt'), open = 'w')
  sink(message_file, type = 'message')
  message(Sys.time())
  fit.mle <- sigex.mlefit(data.ts = data.ts,
                          param = par.default,
                          constraint = constraint,
                          mdl = mdl,
                          method = "bfgs",
                          debug = TRUE,
                          hess = FALSE)
  message(Sys.time())
  sink(type = 'message')
  close(message_file)

  #params
  psi.mle <-           fit.mle[[1]]$par
  par.mle <- sigex.psi2par(psi.mle, mdl, data.ts)
  # Save everything
  save.image(file = file.path(modelPath, "results.RData"))
}

