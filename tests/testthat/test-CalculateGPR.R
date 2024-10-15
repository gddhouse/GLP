test_that("calculate_gpr", {
  counts <- matrix(sample(1:100,20), nrow = 5, ncol = 4)
  rownames(counts) <- c("geneA","geneB","geneC","geneD","geneE")
  gene.info <- calculate_gpr(counts)
})
