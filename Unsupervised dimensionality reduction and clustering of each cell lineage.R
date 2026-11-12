library(dplyr)
library(ggplot2)
library(cowplot)
library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(ggsci)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
options(stringsAsFactors = F)

Myeloid[["percent.mito"]] <- PercentageFeatureSet(Myeloid, pattern = "^MT-")
Myeloid[["percent.ribo"]] <- PercentageFeatureSet(Myeloid, pattern = "^RPS|^RPL")
table(Myeloid$percent.mito<30)
Myeloid <- subset(Myeloid, subset = percent.mito < 30)
table(Myeloid$nGene>500)
Myeloid <- subset(Myeloid, subset = nGene > 500)
table(Myeloid$nUMI>1000)
Myeloid <- subset(Myeloid, subset = nUMI > 1000)

Myeloid <- NormalizeData(Myeloid[, Myeloid$sample %in% sample])
Myeloid <- FindVariableFeatures(Myeloid)
VariableFeaturePlot(Myeloid)
Myeloid <- ScaleData(Myeloid, verbose = T)
Myeloid <- RunPCA(Myeloid, npcs = 50, verbose = FALSE)
Myeloid <- RunHarmony(Myeloid, c("patient"), max_iter = 20)

Myeloid <- FindNeighbors(Myeloid, reduction = "harmony", dims = 1:30, do.plot = T)
Myeloid <- FindClusters(Myeloid, resolution = 1)
Myeloid <- RunUMAP(Myeloid, reduction = "harmony", dims = 1:30)

p1 <- DimPlot(Myeloid, reduction = "umap", group.by = "patient", label = TRUE, cols = celltype[-1], pt.size = 1, raster = T)+ NoLegend()
p2 <- DimPlot(Myeloid, reduction = "umap", group.by = "group", label = TRUE, cols = celltype[-1], pt.size = 1, raster = T) + NoLegend()
plot_grid(p1,p2)

DimPlot(Myeloid,
        reduction = "umap",
        group.by = "celltype2",
        label = F,
        repel = TRUE
)




