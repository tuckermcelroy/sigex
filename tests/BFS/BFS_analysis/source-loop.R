
# Top-level directory path (using relative sigex project pathing)
topDir <- 'tests/BFS/BFS_analysis/'

# Directories containing launch.R files to be run
sourceDirs <- c(
    "JD_ser1_hessF_0101_nreg00_maxit0500_visual_D",
    "JD_ser1_hessF_1201_nreg00_maxit0500_autoArima_E",
    "JD_ser1_hessF_2101_nreg00_maxit0500_fixThetaNegPt8_F",
    "JD_ser1_hessF_2101_nreg00_maxit0500_initmdl_A2",
    "JD_ser1_hessF_2101_nreg12_maxit0500_allreg_B",
    "JD_ser1_hessF_2201_nreg00_maxit0500_initmdl_A3",
    "JD_ser1_hessF_2210_nreg00_maxit0500_initmdl_A5",
    "JD_ser1_hessF_3201_nreg00_maxit0500_initmdl_A4",
    "JD_ser1_hessF_4001_nreg00_maxit0500_initmdl_A6",
    "JD_ser2_hessF_0101_nreg00_maxit0500_visual_D",
    "JD_ser2_hessF_1201_nreg00_maxit0500_autoArima_E",
    "JD_ser2_hessF_2101_nreg00_maxit0500_fixThetaNegPt8_F",
    "JD_ser2_hessF_2101_nreg00_maxit0500_initmdl_A",
    "JD_ser2_hessF_2101_nreg03_maxit0500_signifReg_C",
    "JD_ser2_hessF_2101_nreg12_maxit0500_allreg_B",
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
    "US_ser2_hessF_0101_nreg00_maxit0500_visual_D",
    "US_ser2_hessF_1201_nreg00_maxit0500_autoArima_E",
    "US_ser2_hessF_2101_nreg00_maxit0500_fixThetaNegPt8_F",
    "US_ser2_hessF_2101_nreg00_maxit0500_initmdl_A",
    "US_ser2_hessF_2101_nreg03_maxit0500_signifReg_C",
    "US_ser2_hessF_2101_nreg12_maxit0500_allreg_B"
)

# Check that all files have runMLE flag = TRUE
launchFiles <- file.path(topDir, sourceDirs, 'launch.R')
out <- data.frame(model = sourceDirs, runMLE = NA)
for(i in seq_along(launchFiles)){
  f <- launchFiles[i]
  g <- grepl(pattern = "runMLE <- TRUE", x = readLines(f))
  out$runMLE[i] <- any(g)
}
all(out$runMLE) # If TRUE then all have correct runMLE flag
View(out)



# If any out$runMLE are FALSE we set the runMLE flag TRUE
if(FALSE){
  launchFiles <- file.path(topDir, sourceDirs, 'launch.R')
  for(i in seq_along(launchFiles)){
    f <- launchFiles[i]
    x <- readLines(f)

    # If runMLE already TRUE then skip
    g <- grepl(pattern = "runMLE <- TRUE", x = x)
    if(any(g)) next

    # line with runMLE flag
    mleLine <- grep(pattern = "runMLE <- FALSE", x = x)

    # Print that we are changing flag
    print("Changing runMLE flag in:")
    print(sourceDirs[i])

    # Change line
    x[mleLine] <- "runMLE <- TRUE"

    # Overwrite old run.R file
    writeLines(text = x, f)
  }
}


# Set Debug flag
if(FALSE){
  runFiles <- file.path(topDir, sourceDirs, 'run.R')
  for(i in seq_along(runFiles)){
    f <- runFiles[i]
    x <- readLines(f)

    # Find and replace with desired flag
    x_new <- gsub(pattern     = "debug = TRUE",
                  replacement = "debug = FALSE",
                  x = x)

    # Overwrite old run.R file
    writeLines(text = x_new, f)
  }
}



# !!! CAREFUL: Add hess = FALSE to sigex.mlefit( ) in run.R
if(FALSE){
  runFiles <- file.path(topDir, sourceDirs, 'run.R')
  for(i in seq_along(runFiles)){
    f <- runFiles[i]
    x <- readLines(f)

    # Does hess=FALSE already exist?
    g <- grepl(pattern = "hess = FALSE", x = x)
    if(any(g)) next

    # first line for sigex.mlefit
    firstLine <- grep(pattern = "sigex.mlefit[(]", x = x)

    # last line ending sigex.mlefit
    endLine <- grep(pattern = "[)]$", x = x)
    endLine <- min(endLine[endLine > firstLine])

    # Print that we are inserting hess = FALSE
    print("Inserting hess = FALSE in:")
    print(sourceDirs[i])

    # insert hess=FALSE
    insertLine <- endLine - 1
    if(insertLine <= firstLine) stop("issue with firstLine/lastLine")
    x_new <- c(x[1:(insertLine - 1)],
               "hess = FALSE, ",
               x[(insertLine):length(x)])

    # Overwrite old run.R file
    writeLines(text = x_new, f)
  }
}


# !!!! CAREFUL: Delete results.RData from
if(FALSE){
    for(d in sourceDirs){
      launchPath <- paste0(topDir, d)
      launchFile <- file.path(launchPath, 'launch.R')
      resultsFile <- file.path(launchPath, "results.RData")
      if(file.exists(resultsFile)){
        print("Deleting results.RData from:")
        print(d)
        file.remove(resultsFile)
      }
    }
}


# ---- Main loop sourcing all launch.R files in specified directories ----
for(d in sourceDirs){
  topDir <- 'tests/BFS/BFS_analysis/'
  launchPath <- paste0(topDir, d)
  launchFile <- file.path(launchPath, 'launch.R')
  print(d)
  source(launchFile)
}



# loop over all model directories and check for results.RData
resultExists <- rep(FALSE, length(sourceDirs))
for(i in seq_along(sourceDirs)){
  d <- sourceDirs[i]
  resultsFile <- file.path(topDir, d, "results.RData")
  resultExists[i] <- file.exists(resultsFile)
}
data.frame(sourceDirs, resultExists) |> View()


# update
idx <- which(!resultExists)
sourceDirs <- sourceDirs[idx]

