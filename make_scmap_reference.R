metadata <- read.table("Cell_clusterID_cycle.txt",  header=T)
clust_labs <- read.delim("Cluster_names.txt", sep="\t", header=F)


mat <- read.table("FilteredCounts.csv", sep=",", header=T)
rownames(mat) <- mat[,1]
mat <- mat[,-1]
require("Matrix")
mat <- as(mat, "Matrix")

mat <- mat[!grepl("^MT-", rownames(mat)),] #remove mt genes

rownames(metadata) <- metadata[,1];

require("SingleCellExperiment")
sce <- SingleCellExperiment(assays=list(counts=mat), colData=metadata)

require("scater")
sce <-calculateQCMetrics(sce);
sce <- normalize(sce)

assays(sce)[["logcounts"]] <- as.matrix(assays(sce)[["logcounts"]])
assays(sce)[["counts"]] <- as.matrix(assays(sce)[["counts"]])

rowData(sce)$feature_symbol <- rownames(sce);
colData(sce)$cell_type1 <- clust_labs[sce$Cluster.,3]
require("scmap")
sce <- selectFeatures(sce, suppress_plot = FALSE)
sce <- indexCluster(sce)
set.seed(1)
sce <- indexCell(sce)

assays(sce)[["logcounts"]] <- as(assays(sce)[["logcounts"]], "Matrix")
assays(sce)[["counts"]] <- as(assays(sce)[["counts"]], "Matrix")

saveRDS(sce, "scmap_reference.rds")

