#' LOESS with Tukey’s biweight robust statistical method
#'
#' @param x a vector with Predictor Variable
#' @param y a vector with Response Variable
#' @param span the parameter α which controls the degree of smoothing.
#' @param iterations Number of iterations for outlier detection.
#' @param percentage the parameter controlling the window size for outlier detection.
#'
#' @return A list containing three elements:
#'         1.weight: gene weight for LOESS regression(outliers are 0).
#'         2.local_mads: gene's local median absolute deviation.
#'         3.model: An final robust LOESS regression model.
#' @export
#'
#'
robust_loess <- function(x, y, span = span, iterations = 1, percentage = 0.2) {
  n <- length(y)
  weights <- rep(1, n) # 初始化权重
  loess_fit <- list()

  x_range <- diff(range(x)) * percentage # 计算全局范围的百分比

  for (i in 1:iterations) {
    # 使用当前权重拟合loess模型
    fit <- loess(y ~ x, weights = weights, span = span)
    residuals <- y - predict(fit)

    # 向量化局部范围计算
    x_matrix <- outer(x, x, "-") # 计算所有x之间的差值
    local_matrix <- abs(x_matrix) <= x_range / 2 # 定义局部范围的布尔矩阵

    # 计算局部残差的MAD
    local_mads <- apply(local_matrix, 1, function(row) mad(residuals[row]))

    # 更新权重，只调整正残差的极端值
    normalized_residuals <- residuals / (6 * local_mads)
    weights[normalized_residuals > 1] <- 0
    weights[x > 0.99] <- 0
  }

  # 最终拟合loess模型
  loess_fit$weight <- weights
  loess_fit$local_mads <- local_mads
  loess_fit$model <- loess(y ~ x, weights = weights, span = span)

  return(loess_fit)
}
