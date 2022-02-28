
# Top-level directory path (using relative sigex project pathing)
topDir <- 'tests/BFS/BFS_analysis/'

# Directories containing launch.R files to be run
sourceDirs <- c(
  "JD_A2_model",
  "JD_A3_model",
  "JD_A4_model",
  "JD_A5_model",
  "JD_A6_model",
  "US_A2_model",
  "US_A3_model",
  "US_A4_model",
  "US_A5_model",
  "US_A6_model"
)

# Main loop sourcing all launch.R files in specified directories
for(d in sourceDirs){
  launchPath <- paste0(topDir, d)
  launchFile <- file.path(launchPath, 'launch.R')
  source(launchFile)
}

