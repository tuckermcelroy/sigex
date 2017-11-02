polymul <- function(a,b)
{
	u <- convolve(Re(a),rev(Re(b)),type="open") - convolve(Im(a),rev(Im(b)),type="open")
	v <- convolve(Im(a),rev(Re(b)),type="open") + convolve(Re(a),rev(Im(b)),type="open")
	return(Re(u + 1i*v))
}
