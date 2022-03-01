
# Top-level directory path (using relative sigex project pathing)
topDir <- 'tests/BFS/BFS_analysis/'

# Directories containing launch.R files to be run
sourceDirs <- c(
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

# Check that all files have runMLE flag = TRUE
countFlag <- 0
for(d in sourceDirs){
    if(sum(grepl("runMLE <- TRUE", readLines(launchFile))) > 0){
      countFlag <- countFlag + 1
    }
}
length(sourceDirs)
countFlag

# !!!! CAREFUL: Delete results.RData from
if(FALSE){
    for(d in sourceDirs){
      launchPath <- paste0(topDir, d)
      launchFile <- file.path(launchPath, 'launch.R')
      print("Deleting results.RData from:")
      print(d)
      resultsFile <- file.path(launchPath, "results.RData")
      if(file.exists(resultsFile)) file.remove(resultsFile)
    }
}


# ---- Main loop sourcing all launch.R files in specified directories ----
for(d in sourceDirs){
  launchPath <- paste0(topDir, d)
  launchFile <- file.path(launchPath, 'launch.R')
  source(launchFile)
}


### Check for results.RData files in all subdirectories

# create char vector of all subdirs (ignore all .R files and README)
modelDirsIdx <- grep(pattern = "*_", x = list.files(topDir))
modelDirs <- list.files(topDir)[modelDirsIdx]

# loop over all model directories and check for results.RData
modelDir_withResults <- numeric()
for(d in modelDirs){
  print(d)
  print(file.exists(file.path(topDir, d, "results.RData")))
}

# create vector of all model dirs with results files
idx <- file.exists(file.path(topDir, modelDirs, "results.RData"))
modelDir_withResults <- modelDirs[idx]

