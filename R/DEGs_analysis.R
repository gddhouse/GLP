### calculate DEGs of unsupervised clusters derived from each method
obj <- readRDS("/DATA/TLS/Users/xieyue/01LPR/01figure/00data/00rds/Martins.rds")

load("/DATA/TLS/Users/xieyue/01LPR/01figure/10Final_result_V2/01scRNA/01GLP_v6_0.05/36Martins/cluster_result.RData")

load("/DATA/TLS/Users/xieyue/01LPR/01figure/10Final_result_V2/01scRNA/02OtherMethods/36Martins/cluster_result.RData")

obj@active.ident <- as.factor(obj$annotation)

table(obj@active.ident)

markers <- FindAllMarkers(obj, only.pos = TRUE)

write.csv(markers, "martins_deg.csv")


### merge all DEGs and compare the log2FC
deg <- list.files("./DEGs/")
deg

me <- sapply(deg, function(x) strsplit(x, "_")[[1]][1])

df <- data.frame()
top <- data.frame()
for(i in 1:length(deg)){
    deg_tmp <- read.csv(paste0("./DEGs/", deg[i]))
    deg_tmp$method <- me[i]
    deg_tmp$cluster <- as.character(deg_tmp$cluster)
    df <- rbind(df, deg_tmp)
    top_tmp <- deg_tmp %>% group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% slice_head(n=50)
    top <- rbind(top, top_tmp)
}

df$method <- factor(df$method, levels = c("LPM", "FEAST", 'genebasisR', "HVG", "SCT", "M3Drop", 'NBDrop', "HRG", "SCMarker"))
top$method <- factor(top$method, levels = c("LPM", "FEAST", 'genebasisR', "HVG", "SCT", "M3Drop", 'NBDrop', "HRG", "SCMarker"))

pdf("deg_lf_top50.pdf")
ggplot(top, aes(x=method, y=avg_log2FC, fill=method)) + geom_boxplot() + 
    scale_fill_manual(values = col.df[levels(top$method),'color'])+  
    theme(axis.title.x = element_text(size=20),
         axis.text.x = element_text(size=15, angle = 45, hjust = 0.9),
         axis.title.y = element_text(size=20),
            axis.text.y = element_text(size=15),
          legend.position = "top",
          panel.grid = element_blank(),
          panel.background = element_rect(fill="white", color="black"))+ 
          stat_compare_means(ref.group = "LPM", paired = FALSE,  label = "p.signif", method = "wilcox.test", 
                             method.args = list(alternative = "less"),size=8)
dev.off()
