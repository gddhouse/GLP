#' Calculate BIC index for LOESS model
#'
#' @param loess_obj An object of class "loess"
#'
#' @return A BIC value
#' @export
#'
#' @examples
#'
#' cars.lo <- loess(dist ~ speed, cars)
#' calculate_bic(cars.lo)
#'
calculate_bic <- function(loess_obj) {
  n <- length(loess_obj$fitted)
  rss <- sum(loess_obj$residuals^2)
  k <- loess_obj$trace.hat

  bic <- n * log(rss / n) + log(n) * k
  return(c(BIC = bic))
}
