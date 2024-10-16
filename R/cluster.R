#' Cluster based on input HVG
#'
#' @param counts an single cell expression matrix
#' @param gene_list input HVG
#' @param res Value of the resolution parameter(default is 0.8), bigger value obtain a larger number of clusters
#'
#' @return A seurat object
#' @export
#'
#' @examples
Cluster=function(counts, gene_list,res=0.8){
  SeuObj = CreateSeuratObject(counts=data)
  SeuObj <- NormalizeData(SeuObj,verbose = FALSE)
  all.gene = rownames(SeuObj)
  SeuObj <- ScaleData(SeuObj, features = all.gene, verbose = FALSE)
  if(ncol(data)<50){
    SeuObj <- RunPCA(SeuObj, features = gene_list,verbose = FALSE,npcs=20)
  }else{
    SeuObj <- RunPCA(SeuObj, features = gene_list,verbose = FALSE)
  }
  if(length(gene_list)<10){
    SeuObj <- FindNeighbors(SeuObj, verbose = FALSE, dims = 1:(length(gene_list)-1))
  }else{
    SeuObj <- FindNeighbors(SeuObj, verbose = FALSE)
  }
  SeuObj <- FindClusters(SeuObj, verbose = FALSE, res=res)
  if(ncol(SeuObj@reductions$pca@cell.embeddings)<20){
    SeuObj <- RunUMAP(SeuObj, reduction = "pca", dims = 1:ncol(SeuObj@reductions$pca@cell.embeddings), verbose=F)
  }else{
    SeuObj <- RunUMAP(SeuObj, reduction = "pca", dims = 1:20, verbose=F)
  }
  return(SeuObj)
}
