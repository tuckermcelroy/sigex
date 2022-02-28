zeroReg <- function(psi, mdl){
  outPsi <- psi
  nARMAparams <- sum(mdl$type[[1]][[2]][1:4])
  outPsi[(nARMAparams + 1):length(psi)] <- 0
  return(outPsi)
}
