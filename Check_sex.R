Ygenes <- read.delim("/cluster/projects/macparland/TA/ExternalData/NCBI_Gene_Ychr_NOT_Xchr_Human_10Feb2021.txt")

Ygenes <- Ygenes[,"Symbol"]

Fgenes <- c("XIST")


require(Seurat)
obj <- readRDS("Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno_dimreduce.rds")

expr_mat_M <- obj@assays$RNA@data[rownames(obj) %in% c(Ygenes),]
expr_mat_F <- obj@assays$RNA@data[rownames(obj) %in% c(Fgenes),]

M_score <- Matrix::colMeans(expr_mat_M);

