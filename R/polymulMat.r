polymulMat <- function(amat,bmat)
{
        p <- dim(amat)[3]-1
        q <- dim(bmat)[3]-1
        N <- dim(amat)[2]

        r <- p+q
        bmat.pad <- array(0,c(N,N,r+1))
        for(i in 1:(q+1)) { bmat.pad[,,i] <- bmat[,,i] }
        cmat <- array(0,c(N,N,r+1))
        cmat[,,1] <- amat[,,1] %*% bmat.pad[,,1]
	      if(r > 0) {
          for(j in 2:(r+1))
          {
                cmat[,,j] <- amat[,,1] %*% bmat.pad[,,j]
                if(p > 0) {
                  for(k in 1:min(p,j-1))
                  { cmat[,,j] <- cmat[,,j] + amat[,,k+1] %*% bmat.pad[,,j-k] } }
        } }

        return(cmat)
}

