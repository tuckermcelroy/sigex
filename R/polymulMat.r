polymulMat <- function(amat,bmat)
{
	p <- dim(amat)[3]
	q <- dim(bmat)[3]
	m <- dim(amat)[2]
	amatd <- amat[,,p:1]
	if(q > 1) amatd <- array(c(matrix(0,m,m*(q-1)),amatd),c(m,m,p+q-1))
	bigmat <- NULL
	for(i in 1:(p+q-1)) 
	{
		nextmat <- matrix(amatd[,,1:(p+q-1)],m,m*(p+q-1))
		bigmat <- rbind(nextmat,bigmat)
		amatd <- amatd[,,-1]
		amatd <- array(c(amatd,matrix(0,m,m)),c(m,m,p+q-1))
	}
	bigmat <- bigmat[,1:(m*q)]
	out <- bigmat %*% t(matrix(bmat[,,q:1],m,m*q))
	out <- array(out,c(m,p+q-1,m))
	return(out)
}

