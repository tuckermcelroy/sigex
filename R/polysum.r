polysum <- function(a,b)
{
	# a centered summation, assuming both polys are symmetric power series
        n <- length(a)
        m <- length(b)
        if (m > n) out <- b + c(rep(0,(m-n)/2),a,rep(0,(m-n)/2))
        if (n > m) out <- a + c(rep(0,(n-m)/2),b,rep(0,(n-m)/2))
        if (n==m) out <- a + b
        return(out)
}
