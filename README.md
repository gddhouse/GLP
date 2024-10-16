# GLP: Redefining the high variable genes by LOESS regression with positive ratio

GLP is an R package designed for efficient feature selection in single-cell RNA sequencing (scRNA-seq) analyses. It identifies highly variable genes (HVGs) by analyzing the relationship between the positive ratio and average expression level. The package uses optimized LOESS regression with automatic bandwidth selection based on the Bayesian Information Criterion (BIC) and applies Tukey’s biweight method to reduce the impact of outliers. GLP offers improved accuracy in gene selection, helping to preserve essential biological information for downstream analyses.

## Installation

You can install the development version of GLP from [GitHub](https://github.com/) with:

``` r
install.packages("devtools")
library(devtools)
install_github("gddhouse/GLP")
library(GLP)
```

## Example

This is a basic example which shows you how to solve a common problem:

First，we need to read in data:
```{r example}
library(Seurat)
library(ggplot2)
data <- readRDS("./tests/PBMC.rds")
```


Then, the second step is to calculate expression informations of each gene. Here we offer you `calculate_gpr`:

```{r cars}
counts <- data@assays$RNA@counts
gene.info <- calculate_gpr(counts)
head(gene.info)
```
| exp.mean | positive_rate | pcount |     gene |
|---------:|---------------|-------:|----------|
|         0|              0|       0|MIR1302-10|
|         0|              0|       0|   FAM138A|
|         0|              0|       0|     OR4F5|

The genes with zero expression level will be removed in next step.

The second step is the core step, we will perform robust LOESS regression to identify potential highly variable genes：

```{r pressure, echo = FALSE}
hvg <- glp(df)
reg <- hvg$Genes.regression
ggplot(reg, aes(x=positive_rate, y=exp.mean, color=hvg)) + geom_point() +
    geom_line(data=reg, aes(x=positive_rate, y=fit), color="green",size=1.5)+ 
    coord_cartesian(xlim = c(0, 0.5),ylim=c(0,7))+
    theme(axis.title.x = element_text(size=20),
         axis.text.x = element_text(size=15, angle = 0, hjust = 0.9),
         axis.title.y = element_text(size=20),
            axis.text.y = element_text(size=15),
          legend.position = "top",
          panel.grid = element_blank(),
          panel.background = element_rect(fill="white", color="black"))
```
<picture>
 <source media="(prefers-color-scheme: dark)" srcset="https://github.com/gddhouse/GLP/blob/main/tests/regression.png">
 <source media="(prefers-color-scheme: light)" srcset="https://github.com/gddhouse/GLP/blob/main/tests/regression.png">
 <img alt="YOUR-ALT-TEXT" src="https://github.com/gddhouse/GLP/blob/main/tests/regression.png" width="420" height="400">
</picture>


The potential HVGs selected by GLP are stored in hvg$HVG.
