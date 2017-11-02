install.packages("~/sigex/", repos = NULL, type = "source")
library(sigex)
ls("package:sigex")
lsf.str("package:sigex")

# Check which file has an error
for (f in list.files("~/sigex/R/", full.names=TRUE)){
  print(f)
  parse(f)
}
