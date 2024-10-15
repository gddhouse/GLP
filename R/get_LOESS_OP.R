#' Compute optimal span width for LOESS regression
#'
#' @param x a vector with Predictor Variable
#' @param y a vector with Response Variable
#' @param span_values Default span values
#'
#' @return the optimal span value
#' @export
#'
#'
get_loess_op <- function(x, y,
                         span_values = seq(from = 0.01, to = 0.1, by = 0.01)) {
  bic_values <- numeric(length(span_values))
  for (j in seq_along(span_values)) {
    span <- span_values[j]
    loess_fit <- suppressWarnings({
      loess(y ~ x, span = span)
    })
    bic <- calculate_bic(loess_fit)
    bic_values[j] <- bic["BIC"]
  }
  optimal_span_bic <- span_values[which.min(bic_values)]
  return(optimal_span_bic)
}
