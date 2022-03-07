# Top-level directory path (using relative sigex project pathing)
topDir <- 'tests/BFS/BFS_analysis/'

# List off all directories in top-level directory
list.dirs(topDir, full.names = "FALSE")

# List of models (directories) to load results.RData for
model_directories <-
  c(
    "JD_ser1_2101_nreg00_maxit2500_initmdl_A2",
    "JD_ser1_2101_nreg00_maxit2500_initmdl_A2",
    "JD_ser1_2201_nreg00_maxit2500_initmdl_A3",
    "JD_ser1_3201_nreg00_maxit2500_initmdl_A4",
    "JD_ser1_2210_nreg00_maxit2500_initmdl_A5",
    "JD_ser1_4001_nreg00_maxit2500_initmdl_A6",
    "US_ser1_2101_nreg00_maxit2500_initmdl_A2",
    "US_ser1_2201_nreg00_maxit2500_initmdl_A3",
    "US_ser1_3201_nreg00_maxit2500_initmdl_A4",
    "US_ser1_2210_nreg00_maxit2500_initmdl_A5",
    "US_ser1_4001_nreg00_maxit2500_initmdl_A6"
  )

par(mfrow = c(2, 3), mar = c(5, 4, 4, 2))
for(d in model_directories){
  # Load results file and create short model name
  load(file.path(topDir, d, 'results.RData'))
  print(d)

  # Check convergence of optim
  optimConverge <- fit.mle[[1]]$convergence == 0
  print(optimConverge)
  print(fit.mle[[1]]$convergence)

  # Seasonal Theta (or seasonal AR)
  PAR <- sigex.psi2par(psi = fit.mle[[1]]$par, mdl = mdl, data.ts = data.ts)
  seasTheta <- PAR[[3]][[1]][sum(order[1:2]) + 1]
  seasTheta <- round(seasTheta, 3)

  # Make acf plot of residuals
  if(TRUE){

    substrRight <- function(x, n){
      substr(x, nchar(x)-n+1, nchar(x))
    }
    shortName <- paste(substr(d, 1, 2), substrRight(d, 2), sep = "-")
    plotTitle <- paste(shortName, seasTheta)
    # plot
    resid.mle <- sigex.resid(psi.mle, mdl, data.ts)[[1]]
    resid.mle <- sigex.load(t(resid.mle), start(data.ts), frequency(data.ts), colnames(data.ts), FALSE)
    resid.acf <- acf(resid.mle, lag.max = 4 * 53, plot = TRUE, main = plotTitle)
  }
}
