# Top-level directory path (using relative sigex project pathing)
topDir <- 'tests/BFS/BFS_analysis/'

# List off all directories in top-level directory
list.dirs(topDir, full.names = "FALSE")

g <-
grep("JD_ser1_2101_nreg00_maxit[0-9]{4}",
     x = list.dirs(topDir, full.names = "FALSE"),
     value = TRUE)

g <-
grep("hessF",
     x = list.dirs(topDir, full.names = "FALSE"),
     value = TRUE)



# List of models (directories) to load results.RData for
model_directories <-
  c("JD_ser1_hessF_0101_nreg00_maxit0500_visual_D",
    "JD_ser1_hessF_1201_nreg00_maxit0500_autoArima_E",
    "JD_ser1_hessF_2101_nreg00_maxit0500_initmdl_A2",
    "JD_ser1_hessF_2101_nreg12_maxit0500_allreg_B",
    "JD_ser1_hessF_2201_nreg00_maxit0500_initmdl_A3",
    "JD_ser1_hessF_2210_nreg00_maxit0500_initmdl_A5",
    "JD_ser1_hessF_3201_nreg00_maxit0500_initmdl_A4",
    "JD_ser1_hessF_4001_nreg00_maxit0500_initmdl_A6",
    "JD_ser1_hessT_0101_nreg00_maxit0050_visual_D",
    "JD_ser1_hessT_1201_nreg00_maxit0050_autoArima_E",
    "JD_ser1_hessT_2101_nreg00_maxit0050_fixThetaNegPt8_F",
    "JD_ser1_hessT_2101_nreg00_maxit0500_initmdl_A2",
    "JD_ser1_hessT_2101_nreg00_maxit2500_initmdl_A",
    "JD_ser1_hessT_2101_nreg00_maxit2500_initmdl_A2",
    "JD_ser1_hessT_2101_nreg03_maxit0050_signifReg_C",
    "JD_ser1_hessT_2101_nreg12_maxit0050_allreg_B",
    "JD_ser1_hessT_2201_nreg00_maxit0500_initmdl_A3",
    "JD_ser1_hessT_2201_nreg00_maxit2500_initmdl_A3",
    "JD_ser1_hessT_2210_nreg00_maxit0500_initmdl_A5",
    "JD_ser1_hessT_2210_nreg00_maxit2500_initmdl_A5",
    "JD_ser1_hessT_3201_nreg00_maxit0500_initmdl_A4",
    "JD_ser1_hessT_3201_nreg00_maxit2500_initmdl_A4",
    "JD_ser1_hessT_4001_nreg00_maxit0500_initmdl_A6",
    "JD_ser1_hessT_4001_nreg00_maxit2500_initmdl_A6",
    "JD_ser2_hessF_0101_nreg00_maxit0500_visual_D",
    "JD_ser2_hessF_1201_nreg00_maxit0500_autoArima_E",
    "JD_ser2_hessF_2101_nreg00_maxit0500_fixThetaNegPt8_F",
    "JD_ser2_hessF_2101_nreg00_maxit0500_initmdl_A",
    "JD_ser2_hessF_2101_nreg03_maxit0500_signifReg_C",
    "JD_ser2_hessF_2101_nreg12_maxit0500_allreg_B",
    "JD_ser2_hessT_0101_nreg00_maxit0050_visual_D",
    "JD_ser2_hessT_1201_nreg00_maxit0050_autoArima_E",
    "JD_ser2_hessT_2101_nreg00_maxit0050_fixThetaNegPt8_F",
    "JD_ser2_hessT_2101_nreg00_maxit0050_initmdl_A",
    "JD_ser2_hessT_2101_nreg03_maxit0050_signifReg_C",
    "JD_ser2_hessT_2101_nreg12_maxit0050_allreg_B",
    "US_ser1_hessF_0101_nreg00_maxit0500_visual_D",
    "US_ser1_hessF_1201_nreg00_maxit0500_autoArima_E",
    "US_ser1_hessF_2101_nreg00_maxit0500_fixThetaNegPt8_F",
    "US_ser1_hessF_2101_nreg00_maxit0500_initmdl_A",
    "US_ser1_hessF_2101_nreg00_maxit0500_initmdl_A2",
    "US_ser1_hessF_2101_nreg03_maxit0500_signifReg_C",
    "US_ser1_hessF_2101_nreg12_maxit0500_allreg_B",
    "US_ser1_hessF_2201_nreg00_maxit0500_initmdl_A3",
    "US_ser1_hessF_2210_nreg00_maxit0500_initmdl_A5",
    "US_ser1_hessF_3201_nreg00_maxit0500_initmdl_A4",
    "US_ser1_hessF_4001_nreg00_maxit0500_initmdl_A6",
    "US_ser1_hessT_0101_nreg00_maxit0050_visual_D",
    "US_ser1_hessT_1201_nreg00_maxit0050_autoArima_E",
    "US_ser1_hessT_2101_nreg00_maxit0050_fixThetaNegPt8_F",
    "US_ser1_hessT_2101_nreg00_maxit0050_initmdl_A",
    "US_ser1_hessT_2101_nreg00_maxit2500_initmdl_A",
    "US_ser1_hessT_2101_nreg00_maxit2500_initmdl_A2",
    "US_ser1_hessT_2101_nreg03_maxit0050_signifReg_C",
    "US_ser1_hessT_2101_nreg12_maxit0050_allreg_B",
    "US_ser1_hessT_2201_nreg00_maxit0500_initmdl_A3",
    "US_ser1_hessT_2201_nreg00_maxit2500_initmdl_A3",
    "US_ser1_hessT_2210_nreg00_maxit0500_initmdl_A5",
    "US_ser1_hessT_2210_nreg00_maxit2500_initmdl_A5",
    "US_ser1_hessT_3201_nreg00_maxit0500_initmdl_A4",
    "US_ser1_hessT_3201_nreg00_maxit2500_initmdl_A4",
    "US_ser1_hessT_4001_nreg00_maxit0500_initmdl_A6",
    "US_ser1_hessT_4001_nreg00_maxit2500_initmdl_A6",
    "US_ser2_hessF_0101_nreg00_maxit0500_visual_D",
    "US_ser2_hessF_1201_nreg00_maxit0500_autoArima_E",
    "US_ser2_hessF_2101_nreg00_maxit0500_fixThetaNegPt8_F",
    "US_ser2_hessF_2101_nreg00_maxit0500_initmdl_A",
    "US_ser2_hessF_2101_nreg03_maxit0500_signifReg_C",
    "US_ser2_hessF_2101_nreg12_maxit0500_allreg_B",
    "US_ser2_hessT_0101_nreg00_maxit0050_visual_D",
    "US_ser2_hessT_1201_nreg00_maxit0050_autoArima_E",
    "US_ser2_hessT_2101_nreg00_maxit0050_initmdl_A",
    "US_ser2_hessT_2101_nreg03_maxit0050_signifReg_C",
    "US_ser2_hessT_2101_nreg12_maxit0050_allreg_B"
  )

colNames <- c("shortName",	"Series",	"Hess",	"Delta",	"nReg",	"order",
              "parEst",	"lag51acf",	"lag52acf",	"lag53acf")
df <- data.frame(matrix(NA, length(model_directories), length(colNames)))
colnames(df) <- colNames

# par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
for(i in seq_along(model_directories)){

  d <- model_directories[i]

  ss <- strsplit(d, "_")[[1]]

  df$shortName[i] <- ss[7]
  df$Series[i]    <- ss[2]
  df$Hess[i]      <- ss[3]
  df$Delta[i]     <- ss[1]
  df$nReg[i]      <- ss[5]
  df$order[i]     <- ss[4]

  # Load results file and create short model name
  load(file.path(topDir, d, 'results.RData'))
  print(d)

  # Check convergence of optim
  # print(paste("Optim convergence code = ", fit.mle[[1]]$convergence))
  # optimConverge <- fit.mle[[1]]$convergence == 0

  # Extract par from fit (extra complexity due to constraint)
  eta.mle <- fit.mle[[1]]$par
  psi.mle <- sigex.eta2psi(eta = eta.mle, constraint = constraint)
  par.mle <- sigex.psi2par(psi.mle, mdl, data.ts)
  tsParams <- par.mle[[3]][[1]]

  # Round to 4 decimal places and put in df
  tsParams <- round(tsParams, 4)
  tsParams <- as.character(tsParams)
  tsParams <- paste(tsParams, collapse = ", ")

  df$parEst[i] <- tsParams

  # Seasonal Theta (or seasonal AR)
  # seasTheta <- par.mle[[3]][[1]][sum(order[1:2]) + 1]

  # Find residuals of fit and store the lag 50-54 ACF values of residuals
  resid.mle <- sigex.resid(psi.mle, mdl, data.ts)[[1]]
  resid.mle <- sigex.load(t(resid.mle), start(data.ts), frequency(data.ts), colnames(data.ts), FALSE)
  resid.acf <- acf(resid.mle, lag.max = 4 * 53, plot = FALSE)$acf
  # lag52acf[i, ] <- resid.acf[51:55]

  df$lag51acf[i] <- round(resid.acf[52], 4)
  df$lag52acf[i] <- round(resid.acf[53], 4)
  df$lag53acf[i] <- round(resid.acf[54], 4)

  # Make acf plot of residuals
  if(FALSE){
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









