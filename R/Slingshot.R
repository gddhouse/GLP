setwd("/DATA/TLS/Users/xieyue/01LPR/01figure/10Final_result_V2/08trajectory/annotation/glp/")

obj <- readRDS("/DATA/TLS/Users/xieyue/01LPR/01figure/00data/00rds/Martins.rds")

counts <- obj@assays$RNA@counts

load("/DATA/TLS/Users/xieyue/01LPR/01figure/10Final_result_V2/01scRNA/01GLP_v7_1000_adotped/36Martins/cluster_result.RData")

load("/DATA/TLS/Users/xieyue/01LPR/01figure/10Final_result_V2/01scRNA/01GLP_v6_0.05/36Martins/GLPReg.RData")

sce <- SingleCellExperiment(assays = List(counts = counts))

sce <- sce[glp$HVG, ]

FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)

reducedDims(sce) <- SimpleList(UMAP = result$lpr$umap)
colData(sce)$annotation <- obj$annotation

df <- data.frame(row.names = unique(sce$annotation), seq = 1:7)

plot(reducedDims(sce)$UMAP, col = brewer.pal(9,'Set1')[df[sce$annotation, "seq"]], pch=16, asp = 1)

sce <- slingshot(sce, clusterLabels = 'annotation', reducedDim = 'UMAP', start.clus= 'RGs_DividingPr')

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
plotcol2 <- colors[cut(sce$slingPseudotime_2, breaks=100)]
#plotcol3 <- colors[cut(sce$slingPseudotime_3, breaks=100)]

pdf("slingshot_pseudotime.pdf")
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

plot(reducedDims(sce)$UMAP, col = plotcol2, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

#plot(reducedDims(sce)$UMAP, col = plotcol3, pch=16, asp = 1)
#lines(SlingshotDataSet(sce), lwd=2, col='black')
dev.off()

pdf("slingshot_traject.pdf")
plot(reducedDims(sce)$UMAP, col = brewer.pal(9,'Set3')[df[sce$annotation, "seq"]], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
dev.off()

save(sce, file = "slingshot.Rdata")
