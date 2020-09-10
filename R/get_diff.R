
#' Helper function to allow user to specify integer for trend differencing
#'
#' @param delta integer specifying the number of differences for the trend
#'
#' @return coefficients of polynomial differencing operator

get_trend_delta <- function(delta){

  if(length(delta) == 1){
    delta <-
      switch( delta,
              "1" = c(1, -1),
              "2" = c(1, -2, 1),
              "3" = c(1, -3, 3, -1)
            )
  }
  if(is.null(delta)) stop('delta for trend must 1, 2, 3 or polynomial coefs ')

  return(delta)
}

