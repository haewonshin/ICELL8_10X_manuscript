
# title: "annotation.R"
# output: Read global clustered objects and annotate based on marker genes, recluster after removing junk
# author: Haewon Shin
  
library(Seurat)
library(dplyr)
library(harmony)

seurat <- readRDS("./icell/objects/global_0.9.rds")

print("annotating")
# annotate
fun <- function(x) {
  if (x == "0") {"fibroblast"}
  else if (x == "1") {"endothelial_1"}
  else if (x == "2") {"junk"}
  else if (x == "3") {"macrophage"}
  else if (x == "4") {"pericyte"}
  else if (x == "5") {"cm"}
  else if (x == "6") {"endocardial"}
  else if (x == "7") {"endothelial_2"}
  else if (x == "8") {"t_cell"}
  else if (x == "9") {"junk"}
  else if (x == "10") {"junk"}
  else if (x == "11") {"monocyte"}
  else if (x == "12") {"proliferating"}
  else if (x == "13") {"lymphatic"}
  else if (x == "14") {"mast_cells"}
  else if (x == "15") {"junk"}
}
seurat$cell_type <- mapply(fun, seurat$SCT_snn_res.0.9)

Idents(seurat) <- "cell_type"

# rescale and recluster
seurat <- subset(seurat, idents = "junk", invert = TRUE)

# Normalize
# Set memory limit first
mem.maxVSize()
new_max_vsize <- 1024 * 1024 * 1024 * 1024
mem.maxVSize(new_max_vsize)
options(future.globals.maxSize = 8000*1024^2)

# SCTransform
seurat <- SCTransform(seurat, vars.to.regress = "percent.mito", verbose = FALSE, conserve.memory = TRUE)
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), npcs=100, verbose=TRUE)
seurat <- RunHarmony(seurat, "sample", dims.use = 1:100, reduction = "pca", reduction.save = "harmony", assay.use = "SCT")
seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:100, return.model=TRUE)
seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:100)
seurat <- FindClusters(seurat, graph.name = "SCT_snn", algorithm = 3, resolution = 0.9, verbose = FALSE)

print("saving icell seurat")
saveRDS(seurat, "./icell/objects/global_reclustered.rds")


seurat <- readRDS("./10x/objects/global_0.5.rds")

print("annotating")
# annotate
fun <- function(x) {
  if (x == "0") {"fibroblast_1"}
  else if (x == "1") {"endothelial_1"}
  else if (x == "2") {"macrophage"}
  else if (x == "3") {"pericyte"}
  else if (x == "4") {"cm_1"}
  else if (x == "5") {"endothelial_2"}
  else if (x == "6") {"endocardial"}
  else if (x == "7") {"cm_2"}
  else if (x == "8") {"t_cell"}
  else if (x == "9") {"fibroblast_2"}
  else if (x == "10") {"monocyte"}
  else if (x == "11") {"smc"}
  else if (x == "12") {"junk"}
  else if (x == "13") {"junk"}
  else if (x == "14") {"glia"}
  else if (x == "15") {"mast_cell"}
  else if (x == "16") {"neutrophil"}
  else if (x == "17") {"lymphatic"}
  else if (x == "18") {"fibroblast_3"}
  else if (x == "19") {"junk"}
  else if (x == "20") {"fibroblast_4"}
  else if (x == "21") {"proliferating"}
  else if (x == "22") {"junk"}
  else if (x == "23") {"fibroblast_5"}
  else if (x == "24") {"junk"}
  else if (x == "25") {"epicardial"}
  else if (x == "26") {"junk"}
  
}
seurat$cell_type <- mapply(fun, seurat$SCT_snn_res.0.5)

Idents(seurat) <- "cell_type"

# rescale and recluster
seurat <- subset(seurat, idents = "junk", invert = TRUE)

# SCTransform
seurat <- SCTransform(seurat, vars.to.regress = "percent.mito", verbose = FALSE, conserve.memory = TRUE)
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), npcs=100, verbose=TRUE)
seurat <- RunHarmony(seurat, "sample", dims.use = 1:100, reduction = "pca", reduction.save = "harmony", assay.use = "SCT")
seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:100, return.model=TRUE)
seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:100)
seurat <- FindClusters(seurat, graph.name = "SCT_snn", algorithm = 3, resolution = 0.5, verbose = FALSE)

print("saving 10X seurat")
saveRDS(seurat, "./10x/objects/global_reclustered.rds")