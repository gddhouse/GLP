#' GLP selection of high variable genes
#'
#' @param df a data.table where each row represents a gene,
#'           with four columns: average expression level, positive ratio,
#'           positive cell count, and gene name.
#' @param npc Minimum number of cells in which a gene must be expressed.
#' @param nfeatures the number of high variable genes.
#'
#' @return A list containing three elements:
#'         1.Genes.summary: The input data.table with summary statistics for each gene.
#'         2.Genes.regression: A data.table with the main results of the robust LOESS regression.
#'         3.HVG: A vector of highly variable genes identified through the analysis."
#' @importFrom stats loess
#' @importFrom stats mad
#' @importFrom stats pnorm
#' @importFrom stats predict
#' @importFrom stats sd
#' @export
#'
#' @examples
#' library(GLP)
#' counts <- matrix(sample(0:10,3000*1000, replace = TRUE), nrow = 3000, ncol = 1000)
#' rownames(counts) <- paste0("gene",1:3000)
#' gene.info <- calculate_gpr(counts)
#' result <- glp(gene.info, 3, 1000)
#'
glp <- suppressWarnings(function(df, npc = 3, nfeatures = 1000) {
  result <- list()
  # 过滤 pcount >= 3 的数据
  result[["Genes.summary"]] <- df
  df <- df[pcount >= npc]
  #------------------------------------------------- loess
  optimal_span_bic <- get_loess_op(df$positive_rate, df$exp.mean)
  if (length(optimal_span_bic) == 0) {
    optimal_span_bic <- get_loess_op(df$positive_rate,
                                     df$exp.mean,
                                     span_values = seq(from = 0.1, to = 0.9, by = 0.1))
  }
  loess_fit <- robust_loess(df$positive_rate,
                            df$exp.mean,
                            span = optimal_span_bic)
  df[, fit := predict(loess_fit$model)]
  df[, weight := loess_fit$weight]
  df[, local_mads := loess_fit$local_mads]
  df[, span := optimal_span_bic]
  df[, se := exp.mean - fit] # fit.se
  #------------------------------------------------- calculate confidence interval
  n <- 100
  df[, `:=`(se.sd = NA_real_, se.mean = NA_real_)] # 初始化新列
  df[, c("se.sd", "se.mean") := {
    start_idx <- pmax(1, .I - n)
    end_idx <- pmin(nrow(df), .I + n)
    sample <- df[start_idx:end_idx]
    list(sd(sample$se), mean(sample$se))
  }, by = seq_len(nrow(df))]
  df[, adjust.se := (se - se.mean) / se.sd]
  df[, pv := 1 - pnorm(adjust.se)]
  #------------------------------------------------- define hvg
  df <- df[order(pv)] # 按pv列从小到大排序
  df[, hvg := ifelse(.I <= nfeatures, "yes", "no")] # 根据排序后的行号取前1000个
  #------------------------------------------------- merge all genes
  hvg <- df[hvg == "yes", gene]
  result[["HVG"]] <- hvg
  result[["Genes.regression"]] <- df
  return(result)
})
