# Top-level directory path (using relative sigex project pathing)
topDir <- 'tests/BFS/BFS_analysis/'

# List off all directories in top-level directory
list.dirs(topDir, full.names = "FALSE")

# List of models (directories) to load results.RData for
model_directories <-
  c(
  "JD_A_model",
  "JD_A2_model",
  "JD_A3_model",
  "JD_A4_model",
  "JD_A5_model",
  "JD_A6_model",
  "US_A_model",
  "US_A3_model",
  "US_A3_model",
  "US_A4_model",
  "US_A5_model",
  "US_A6_model"
  )


par(mfrow = c(2, 3), mar = c(5, 4, 4, 2))
for(d in model_directories){
  # Load results file and create short model name
  load(file.path(topDir, d, 'results.RData'))
  shortMdlName <- substr(d, 1, 5)
  print(shortMdlName)

  # Check convergence of optim
  optimConverge <- fit.mle[[1]]$convergence == 0
  print(optimConverge)
  print(fit.mle[[1]]$convergence)

  # Seasonal Theta (or seasonal AR)
  PAR <- sigex.psi2par(psi = fit.mle[[1]]$par, mdl = mdl, data.ts = data.ts)
  seasTheta <- PAR[[3]][[1]][sum(order[1:2]) + 1]
  seasTheta <- round(seasTheta, 3)

  # Make acf plot of residuals
  plotTitle <- paste(shortMdlName, seasTheta)
  # plot
  resid.mle <- sigex.resid(psi.mle, mdl, data.ts)[[1]]
  resid.mle <- sigex.load(t(resid.mle), start(data.ts), frequency(data.ts), colnames(data.ts), FALSE)
  resid.acf <- acf(resid.mle, lag.max = 4 * 53, plot = TRUE, main = plotTitle)
}
