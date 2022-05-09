# Top-level directory path (using relative sigex project pathing)
topDir <- 'tests/BFS/BFS_analysis/'

# List off all directories in top-level directory
list.dirs(topDir, full.names = "FALSE")

grep("JD_ser1_2101_nreg00_maxit[0-9]{4}",
     x = list.dirs(topDir, full.names = "FALSE"),
     value = TRUE)


# List of models (directories) to load results.RData for
model_directories <-
  c(
    "JD_ser2_2101_nreg00_maxit0050_initmdl_A",
    "JD_ser2_2101_nreg12_maxit0050_allReg_B",
    "JD_ser2_2101_nreg03_maxit0050_signifReg_C",
    "JD_ser2_0101_nreg00_maxit0050_visual_D"
  )

par(mfrow = c(2, 2), mar = c(5, 4, 4, 2))
for(d in model_directories){
  # Load results file and create short model name
  load(file.path(topDir, d, 'results.RData'))
  print(d)

  # Check convergence of optim
  print(paste("Optim convergence code = ", fit.mle[[1]]$convergence))
  # optimConverge <- fit.mle[[1]]$convergence == 0

  # Seasonal Theta (or seasonal AR)
  PAR <- sigex.psi2par(psi = fit.mle[[1]]$par, mdl = mdl, data.ts = data.ts)
  seasTheta <- PAR[[3]][[1]][sum(order[1:2]) + 1]
  # seasTheta <- round(seasTheta, 3)
  message("Seasonal theta = ", seasTheta)
  print(PAR)
  print("-------------------------------------------")

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








