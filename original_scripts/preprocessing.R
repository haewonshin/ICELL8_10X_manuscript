
# title: "preprocessing.R"
# output: Create Seurat objects from count matrices.
# author: Haewon Shin


library(Seurat)
library(dplyr)
library(harmony)

samples <- c("D1118","D1192","D1340","T1227","T1342","T1346")
conditions <- c("donor","donor","donor","nicm","nicm","nicm")

# 10X object
print("10X")

list <- list()

for (i in seq_along(samples)) {
  print(samples[[i]])
  directory <- paste0("./samples_10x/", samples[[i]], "/filtered_feature_bc_matrix/")
  print(directory)
  counts <- Read10X(data.dir = directory)

  list[[i]] <- CreateSeuratObject(counts = counts)
  list[[i]]$sample <- samples[[i]]
  list[[i]]$condition <- conditions[[i]]
}

print("merging 10X")
x <- merge(list[[1]], y = c(list[[2]],list[[3]],list[[4]],list[[5]],list[[6]]), add.cell.ids = c("D1118","D1192","D1340","T1227","T1342","T1346"))

# add percent mitochondria + other genes
x <- PercentageFeatureSet(x, pattern = "^MT-", col.name = "percent.mito")
x <- PercentageFeatureSet(x, pattern = "^RPS|^RPL", col.name = "percent.ribo")
x <- PercentageFeatureSet(x, pattern = "^HB[AB]", col.name = "percent.hemo")

print("saving 10X")
saveRDS(x, "./10x/objects/preQC.rds")

# ICELL object
print("ICELL")

gene_info <- read.table(file = "./samples_icell/D1118/analyze/gene_info.csv", sep = ",", header = TRUE)

list <- list()

for (i in seq_along(samples)) {
  print(samples[[i]])
  counts <- read.table(file = paste0("./samples_icell/", samples[[i]], "/analyze/analyze_genematrix.csv"), sep = ",", header = TRUE)
  
  # Now just replace with the gene name so that counts can be summed up
  gene_names <- gene_info$Gene_Name
  counts$GeneID <- gene_names
  
  counts_df <- as.data.frame(counts)
  print("summing")
  counts_df <- counts_df %>%
    group_by(GeneID) %>%
    summarise(across(everything(), sum), .groups = "drop")
  
  counts <- as.data.frame(counts_df)
  
  # Format gene names for Seurat in next step and append counts to list
  rownames(counts) <- counts$GeneID
  counts$GeneID <- NULL
  
  print("making seurat object")
  list[[i]] <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 300)
  list[[i]]$sample <- samples[[i]]
  list[[i]]$condition <- conditions[[i]]
}

print("merging ICELL")
icell <- merge(list[[1]], y = c(list[[2]],list[[3]],list[[4]],list[[5]],list[[6]]), add.cell.ids = c("D1118","D1192","D1340","T1227","T1342","T1346"))

# add percent mitochondria + other genes
icell <- PercentageFeatureSet(icell, pattern = "^MT-", col.name = "percent.mito")
icell <- PercentageFeatureSet(icell, pattern = "^RPS|^RPL", col.name = "percent.ribo")
icell <- PercentageFeatureSet(icell, pattern = "^HB[AB]", col.name = "percent.hemo")

print("saving ICELL")
saveRDS(icell, "./icell/objects/preQC.rds")