#' Compute relevant expression statistics for each gene.
#'
#' @param counts single-cell expression matrix where rows represent genes
#'               and columns represent cells.
#'
#' @return Return a data.table where each row represents a gene,
#'         with four columns: average expression level,
#'         positive ratio, positive cell count, and gene name.
#' @importFrom data.table as.data.table
#' @importFrom data.table data.table
#' @importFrom data.table setnames
#' @importFrom data.table :=
#' @export
#'
#'
#' @examples
#' library(GLP)
#' counts <- matrix(sample(0:10,3000*1000, replace = TRUE), nrow = 3000, ncol = 1000)
#' rownames(counts) <- paste0("gene",1:3000)
#' gene.info <- calculate_gpr(counts)
#'
calculate_gpr <- function(counts) {
  # 转换数据框为 data.table
  data <- as.data.table(counts, keep.rownames = TRUE)
  # 过滤掉全为零的基因
  data <- data[rowSums(data > 0) >= 1, , drop = FALSE]
  gene <- data$rn
  data <- data[, -1]
  # 计算细胞数量
  cellnumber <- ncol(data)
  # 计算每个基因的总表达量、均值和标准差
  sum_exp <- rowSums(data)
  exp.mean <- sum_exp / cellnumber
  # 计算每个基因的正值细胞数量和比例
  counts_matrix_pos <- data > 0
  positive_cell_number <- rowSums(counts_matrix_pos)
  positive_cell_rate <- positive_cell_number / cellnumber
  # 创建数据框
  df <- data.table(
    exp.mean = exp.mean,
    positive_rate = positive_cell_rate,
    pcount = positive_cell_number,
    gene = gene
  )
  # 按照 positive_rate 排序
  df <- df[order(positive_cell_rate)]
  # 将数据框的列名调整为原始的
  setnames(df, c("exp.mean", "positive_rate", "pcount", "gene"))
  # 返回数据框
  return(df)
}
