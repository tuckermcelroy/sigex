install.packages("~/GitHub/sigex/", repos = NULL, type = "source")
library(sigex)
# list all functions from loading sigex package
ls("package:sigex")
lsf.str("package:sigex")

# Check which file has an error
for (f in list.files("~/GitHub/sigex/R/", full.names=TRUE)){
  print(f)
  parse(f)
}

