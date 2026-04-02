
# title: "clustering_global.R"
# output: Apply QC metrics to unfiltered object, normalize, scale, and cluster.
# author: Haewon Shin

library(Seurat)
library(dplyr)
library(harmony)

# Normalize
# Set memory limit first
mem.maxVSize()
new_max_vsize <- 1024 * 1024 * 1024 * 1024  
mem.maxVSize(new_max_vsize)
options(future.globals.maxSize = 8000*1024^2)

# read in 10X
x <- readRDS("./10x/objects/preQC.rds")

# filtering
x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 35000 & percent.mito < 5)
# SCTransform, regressing by mito %
x <- SCTransform(x, vars.to.regress = "percent.mito", verbose = FALSE, conserve.memory = TRUE)
saveRDS(x, file = "./10x/objects/post_sct.rds")
x <- RunPCA(x, features = VariableFeatures(object = x), npcs=100, verbose=TRUE)
x <- RunHarmony(x, "sample", dims.use = 1:100, reduction = "pca", reduction.save = "harmony", assay.use = "SCT")
x <- RunUMAP(x, reduction = "harmony", dims = 1:100, return.model=TRUE)
x <- FindNeighbors(x, reduction = "harmony", dims = 1:100)
x <- FindClusters(x, graph.name = "SCT_snn", algorithm = 3, resolution = 0.5, verbose = TRUE)

saveRDS(x, file = "./10x/objects/global_0.5.rds")

Idents(x) <- "SCT_snn_res.0.5"
DefaultAssay(x) <- 'SCT'
x <- PrepSCTFindMarkers(x, assay = "SCT", verbose = TRUE)
rna.markers <- FindAllMarkers(x, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers, file ="./10x/marker_tables/global_0.5.csv", quote = FALSE)
saveRDS(x, file = "./10x/objects/global_0.5.rds")

# read in icell
x <- readRDS("./icell/objects/preQC.rds")

x <- subset(x, subset = nCount_RNA > 10000 & nCount_RNA < 1e6 & percent.mito < 30)
# SCTransform, regressing by mito %
x <- SCTransform(x, vars.to.regress = "percent.mito", verbose = FALSE, conserve.memory = TRUE)
saveRDS(x, file = "./icell/objects/post_sct.rds")
x <- RunPCA(x, features = VariableFeatures(object = x), npcs=100, verbose=TRUE)
x <- RunHarmony(x, "sample", dims.use = 1:100, reduction = "pca", reduction.save = "harmony", assay.use = "SCT")
x <- RunUMAP(x, reduction = "harmony", dims = 1:100, return.model=TRUE)
x <- FindNeighbors(x, reduction = "harmony", dims = 1:100)
x <- FindClusters(x, graph.name = "SCT_snn", algorithm = 3, resolution = 0.9, verbose = TRUE)

Idents(x) <- "SCT_snn_res.0.9"
DefaultAssay(x) <- 'SCT'
x <- PrepSCTFindMarkers(x, assay = "SCT", verbose = TRUE)
rna.markers <- FindAllMarkers(x, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers, file ="./icell/marker_tables/global_0.9.csv", quote = FALSE)
saveRDS(x, file = "./icell/objects/global_0.9.rds")
