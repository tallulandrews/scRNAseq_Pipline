require(Seurat)
require(ggplot2)
require(M3Drop)
source("~/scripts/My_R_Scripts.R")

# Pearson Residuals

Slices <- c(
				"BA1_1", "BA1_2", "BA1_4",
				"BA2_1", "BA2_2", "BA2_3", "BA2_4",
				"NL_1", "NL_2", "NL_3", "NL_4",
				"NonBA_1", "NonBA_2", "NonBA_3", "NonBA_4"
				)
Mt_gene_pattern <- "^MT-";
norm_method="Pearson"
integration_method="RPCA"
#class_name <- "BA1"
class_name <- commandArgs(trailingOnly=TRUE)
outname <- paste(class_name, norm_method, integration_method, sep="_")

Slices <- Slices[grepl(class_name, Slices)]

print_summary_stats <- function(slice_obj) {
	slice_obj[["percent.mt"]] <- PercentageFeatureSet(slice_obj, pattern = "^MT-")
	slice_obj[["percent.ribo"]] <- PercentageFeatureSet(slice_obj, pattern = "^RP[SL]")

	ng = nrow(slice_obj)
	print(paste("n Genes =", ng))
	nc = ncol(slice_obj)
	print(paste("n Spots =", nc))
	tc = sum(slice_obj@meta.data$nCount_Spatial)
	print(paste("Total UMIs =", tc))
	umi_pC = median(slice_obj@meta.data$nCount_Spatial)
	print(paste("UMI/spot =", umi_pC))
	gene_pc = median(slice_obj@meta.data$nFeature_Spatial)
	print(paste("Gene/spot =", gene_pc))
	pct_mt = sum(slice_obj@meta.data$nCount_Spatial*slice_obj@meta.data$percent.mt/100)/tc*100
	print(paste("% Mito =", pct_mt))
	pct_ribo = sum(slice_obj@meta.data$nCount_Spatial*slice_obj@meta.data$percent.ribo/100)/tc*100
	print(paste("% Ribo =", pct_ribo))
}

simpson_index <- function(samples, clusters) {
	tab <- table(samples, clusters);
	tab <- t(t(tab)/colSums(tab))
	return(colSums(tab^2))
}

simpson_index_optimal <- function(samples) {
	tab <- table(samples)
	tab <- tab/sum(tab)
	return(sum(tab^2))
}

simpson_plot <- function(object, samples, clusters, verbose=FALSE, sample.colours=NULL) {
	if (length(samples) == 1) {
		if (verbose) {print("Collected samples from metadata.")}
		samples <- object@meta.data[,samples]
	}
	if (length(clusters) == 1) {
		if (verbose) {print("Collected clusters from metadata.")}
		clusters <- object@meta.data[,clusters]
	}
	toplot <- table(samples, clusters)

	if (is.null(sample.colours)) {
		sample.colours <- RColorBrewer::brewer.pal(nrow(toplot))
	}
	simpson.cluster <- round(simpson_index(samples, clusters), digits=2)
	simpson.optimal <- round(simpson_index_optimal(samples), digits=2)
	b_loc <- barplot(toplot/ncol(object), col=sample.colours)
	text(x=b_loc, y=colSums(toplot/ncol(object)), labels=simpson, pos=3)
	legend("topright", c(
		paste("Exp. Simpson =", simpson.optimal, digits=2),
		paste("Avg. Simpson =", round(mean(simpson), digits=2)),
		paste("Max. Simpson =", max(simpson))
	), bty="n")
}

default_clustering <- function(object) {
	object <- RunPCA(object, assay="SCT", verbose = FALSE, features=consensus_features)
	object <- FindNeighbors(object, assay="SCT", reduction = "pca", dims = 1:20)
	object <- RunUMAP(object, assay="SCT", reduction = "pca", dims = 1:20)
	object <- FindClusters(object)
	return(object)
}


set.seed(101)

obj_list <- list();
merged_obj1 <- NULL;
for(slice in Slices) {
	class_name <- unlist(strsplit(slice, "_"))[1]
	obj <- Load10X_Spatial(slice)
	obj <- obj[rowSums(obj@assays$Spatial@counts) > 5,]
	print(slice);
#	print_summary_stats(obj);
	
	obj@meta.data$orig.ident <- rep(slice, ncol(obj));
	obj@meta.data$cell_barcode <- colnames(obj);
	obj@meta.data$cell_id <- paste(obj@meta.data$orig.ident, obj@meta.data$cell_barcode, sep="_")
	obj@meta.data$sample_set <- rep(class_name, nrow(obj@meta.data))
	obj <- FindVariableFeatures(obj)
	
	pearson <- NBumiPearsonResiduals(obj@assays$Spatial@counts)
	obj@assays$Spatial@data <- as(pearson,"dgCMatrix")

	obj_list[[slice]] <- obj

	if (is.null(merged_obj1)) {
		merged_obj1 <- obj;
	} else {
		merged_obj1 <- merge(merged_obj1, obj, add.cell.ids = c("", slice));
	}
}

consensus_features <- my_consensus_features(obj_list)

# Version 1: Individual Pearson
set.seed(897)
local_pearson <- merged_obj1@assays$Spatial@data
merged_obj1@assays$Spatial@scale.data <- as.matrix(merged_obj1@assays$Spatial@data)
merged_obj1 <- RunPCA(merged_obj1, assay="Spatial", verbose = FALSE, features=consensus_features)
merged_obj1 <- FindNeighbors(merged_obj1, assay="Spatial", reduction = "pca", dims = 1:20)
merged_obj1 <- RunUMAP(merged_obj1, assay="Spatial", reduction = "pca", dims = 1:20)
merged_obj1 <- FindClusters(merged_obj1)
cluster_comp <- table(merged_obj1@meta.data$orig.ident, merged_obj1@meta.data$seurat_clusters)
simpson <- round(simpson_index(merged_obj1@meta.data$orig.ident, merged_obj1@meta.data$seurat_clusters), digits=2)
pdf(paste(outname, "UMAP_LocalPearson.pdf", sep="_"), width=5, height=6.5)
b_loc <- barplot(cluster_comp/ncol(merged_obj1), col=RColorBrewer::brewer.pal(nrow(cluster_comp), "Set3"))
text(x=b_loc, y=colSums(cluster_comp/ncol(merged_obj1)), labels=simpson, pos=3)
legend("topright", c(
		paste("Avg. Simpson =", round(mean(simpson), digits=2)),
		paste("Max. Simpson =", max(simpson))
	), bty="n")
print(DimPlot(merged_obj1, reduction="umap", group.by="orig.ident"))
dev.off()


# Version 2: Global Pearson
set.seed(897)
global_pearson <- NBumiPearsonResiduals(merged_obj1@assays$Spatial@counts)
merged_obj1@assays$Spatial@scale.data <- as.matrix(global_pearson)
merged_obj1 <- RunPCA(merged_obj1, assay="Spatial", verbose = FALSE, features=consensus_features)
merged_obj1 <- FindNeighbors(merged_obj1, assay="Spatial", reduction = "pca", dims = 1:20)
merged_obj1 <- RunUMAP(merged_obj1, assay="Spatial", reduction = "pca", dims = 1:20)
merged_obj1 <- FindClusters(merged_obj1)
cluster_comp <- table(merged_obj1@meta.data$orig.ident, merged_obj1@meta.data$seurat_clusters)
pdf(paste(outname, "UMAP_GlobalPearson.pdf", sep="_"), width=5, height=6.5)
b_loc <- barplot(cluster_comp/ncol(merged_obj1), col=RColorBrewer::brewer.pal(nrow(cluster_comp), "Set3"))
simpson <- round(simpson_index(merged_obj1@meta.data$orig.ident, merged_obj1@meta.data$seurat_clusters), digits=2)
text(x=b_loc, y=colSums(cluster_comp/ncol(merged_obj1)), labels=simpson, pos=3)
legend("topright", c(
		paste("Avg. Simpson =", round(mean(simpson), digits=2)),
		paste("Max. Simpson =", max(simpson))
	), bty="n")
print(DimPlot(merged_obj1, reduction="umap", group.by="orig.ident"))
dev.off()

# Version 3: SCT
set.seed(897)
merged_obj1 <- SCTransform(merged_obj1, assay="Spatial", vars.to.regress=c("orig.ident", "nFeature_Spatial"))
merged_obj1 <- ScaleData(merged_obj1, assay="SCT")
merged_obj1 <- RunPCA(merged_obj1, assay="SCT", verbose = FALSE, features=consensus_features)
merged_obj1 <- FindNeighbors(merged_obj1, assay="SCT", reduction = "pca", dims = 1:20)
merged_obj1 <- RunUMAP(merged_obj1, assay="SCT", reduction = "pca", dims = 1:20)
merged_obj1 <- FindClusters(merged_obj1)
cluster_comp <- table(merged_obj1@meta.data$orig.ident, merged_obj1@meta.data$seurat_clusters)
pdf(paste(outname, "UMAP_SCT.pdf", sep="_"), width=5, height=6.5)
b_loc <- barplot(cluster_comp/ncol(merged_obj1), col=RColorBrewer::brewer.pal(nrow(cluster_comp), "Set3"))
simpson <- round(simpson_index(merged_obj1@meta.data$orig.ident, merged_obj1@meta.data$seurat_clusters), digits=2)
text(x=b_loc, y=colSums(cluster_comp/ncol(merged_obj1)), labels=simpson, pos=3)
legend("topright", c(
		paste("Avg. Simpson =", round(mean(simpson), digits=2)),
		paste("Max. Simpson =", max(simpson))
	), bty="n")
print(DimPlot(merged_obj1, reduction="umap", group.by="orig.ident"))
dev.off()

# Version 4: LogNormalize
set.seed(897)
merged_obj1 <- NormalizeData(merged_obj1, normalization.method = "LogNormalize", assay="Spatial")
merged_obj1@assays$Spatial@scale.data <- matrix(0);
merged_obj1 <- ScaleData(merged_obj1, assay="Spatial")
merged_obj1 <- RunPCA(merged_obj1, assay="Spatial", verbose = FALSE, features=consensus_features)
merged_obj1 <- FindNeighbors(merged_obj1, assay="Spatial", reduction = "pca", dims = 1:20)
merged_obj1 <- RunUMAP(merged_obj1, assay="Spatial", reduction = "pca", dims = 1:20)
merged_obj1 <- FindClusters(merged_obj1, graph.name="Spatial_snn")
cluster_col <- "seurat_clusters"
cluster_comp <- table(merged_obj1@meta.data$orig.ident, merged_obj1@meta.data[,cluster_col])
simpson <- round(simpson_index(merged_obj1@meta.data$orig.ident, merged_obj1@meta.data[,cluster_col]), digits=2)
pdf(paste(outname, "UMAP_LogNormalize.pdf", sep="_"), width=5, height=6.5)
b_loc <- barplot(cluster_comp/ncol(merged_obj1), col=RColorBrewer::brewer.pal(nrow(cluster_comp), "Set3"))
text(x=b_loc, y=colSums(cluster_comp/ncol(merged_obj1)), labels=simpson, pos=3)
legend("topright", c(
		paste("Avg. Simpson =", round(mean(simpson), digits=2)),
		paste("Max. Simpson =", max(simpson))
	), bty="n")
print(DimPlot(merged_obj1, reduction="umap", group.by="orig.ident"))
print(DimPlot(merged_obj1, reduction="umap", group.by=cluster_col))
dev.off()



### Integration ###
obj_list2 <- SplitObject(merged_obj1, split.by="orig.ident") 

# Pipeline
set.seed(291) 

## SCT ##
for (integration_method in c("cca", "rpca")) {
	set.seed(897)
	if (integration_method == "rpca") {
		for (i in 1:length(obj_list2)) {
			obj_list2[[i]] <- RunPCA(obj_list2[[i]], assay="SCT", features=consensus_features)
		}
	}
	anchors <- FindIntegrationAnchors(obj_list2, assay=rep("SCT", length(obj_list)), 
			anchor.features = consensus_features, reduction = integration_method)

	integrated_obj <- IntegrateData(anchorset=anchors, features=consensus_features, normalization.method="SCT")
	DefaultAssay(integrated_obj) <- "integrated"
	integrated_obj <- RunPCA(integrated_obj, assay="integrated", verbose = FALSE, features=consensus_features)
	integrated_obj <- FindNeighbors(integrated_obj, assay="integrated", reduction = "pca", dims = 1:20)
	integrated_obj <- RunUMAP(integrated_obj, assay="integrated", reduction = "pca", dims = 1:20)
	integrated_obj <- FindClusters(integrated_obj, graph.name="integrated_snn")
	cluster_comp <- table(integrated_obj@meta.data$orig.ident, integrated_obj@meta.data$seurat_clusters)
	simpson <- round(simpson_index(integrated_obj@meta.data$orig.ident, integrated_obj@meta.data$seurat_clusters), digits=2)
	pdf(paste(outname, integration_method, "UMAP_SCT_Integrated.pdf", sep="_"), width=5, height=6.5)
	b_loc <- barplot(cluster_comp/ncol(integrated_obj), col=RColorBrewer::brewer.pal(nrow(cluster_comp), "Set3"))
	text(x=b_loc, y=colSums(cluster_comp/ncol(integrated_obj)), labels=simpson, pos=3)
	legend("topright", c(
			paste("Avg. Simpson =", round(mean(simpson), digits=2)),
			paste("Max. Simpson =", max(simpson))
		), bty="n")
	print(DimPlot(integrated_obj, reduction="umap", group.by="orig.ident"))
	if (integration_method == "cca") {
		saveRDS(integrated_obj, paste(outname, integration_method, "UMAP_SCT_Integrated.rds", sep="_"))
	}
	dev.off()
}

## LogNormalize ##
for (integration_method in c("cca", "rpca")) {
	set.seed(897)
	if (integration_method == "rpca") {
		for (i in 1:length(obj_list2)) {
			obj_list2[[i]] <- ScaleData(obj_list2[[i]], assay="Spatial", features=consensus_features)
			obj_list2[[i]] <- RunPCA(obj_list2[[i]], assay="Spatial", features=consensus_features)
		}
	}
	anchors <- FindIntegrationAnchors(obj_list2, assay=rep("Spatial", length(obj_list)), 
			anchor.features = consensus_features, reduction = integration_method)

	integrated_obj <- IntegrateData(anchorset=anchors, features=consensus_features, normalization.method="LogNormalize")
	DefaultAssay(integrated_obj) <- "integrated"
	integrated_obj <- ScaleData(integrated_obj, assay="integrated", verbose = FALSE, features=consensus_features)
	integrated_obj <- RunPCA(integrated_obj, assay="integrated", verbose = FALSE, features=consensus_features)
	integrated_obj <- FindNeighbors(integrated_obj, assay="integrated", reduction = "pca", dims = 1:20)
	integrated_obj <- RunUMAP(integrated_obj, assay="integrated", reduction = "pca", dims = 1:20)
	integrated_obj <- FindClusters(integrated_obj, graph.name="integrated_snn")
	cluster_comp <- table(integrated_obj@meta.data$orig.ident, integrated_obj@meta.data$seurat_clusters)
	pdf(paste(outname, integration_method, "UMAP_LogNormalize_Integrated.pdf", sep="_"), width=5, height=6.5) 
	b_loc <- barplot(cluster_comp/ncol(integrated_obj), col=RColorBrewer::brewer.pal(nrow(cluster_comp), "Set3"))
	simpson <- round(simpson_index(integrated_obj@meta.data$orig.ident, integrated_obj@meta.data$seurat_clusters), digits=2)
	text(x=b_loc, y=colSums(cluster_comp/ncol(integrated_obj)), labels=simpson, pos=3)
	legend("topright", c(
			paste("Avg. Simpson =", round(mean(simpson), digits=2)),
			paste("Max. Simpson =", max(simpson))
		), bty="n")
	print(DimPlot(integrated_obj, reduction="umap", group.by="orig.ident"))
	dev.off()
}


## Local Pearson ##
merged_obj1@assays$SCT@counts <- merged_obj1@assays$Spatial@counts
merged_obj1@assays$SCT@data <- local_pearson
merged_obj1@assays$SCT@scale.data <- as.matrix(local_pearson)
obj_list2 <- SplitObject(merged_obj1, split.by="orig.ident") 
for (integration_method in c("cca", "rpca")) {
	set.seed(897)
	if (integration_method == "rpca") {
		for (i in 1:length(obj_list2)) {
			obj_list2[[i]] <- RunPCA(obj_list2[[i]], assay="SCT", features=consensus_features)
		}
	}
	anchors <- FindIntegrationAnchors(obj_list2, assay=rep("SCT", length(obj_list2)), 
			anchor.features = consensus_features, reduction = integration_method)

	integrated_obj <- IntegrateData(anchorset=anchors, features=consensus_features, normalization.method="SCT")
	DefaultAssay(integrated_obj) <- "integrated"
	integrated_obj <- RunPCA(integrated_obj, verbose = FALSE, features=consensus_features)
	integrated_obj <- FindNeighbors(integrated_obj, reduction = "pca", dims = 1:20)
	integrated_obj <- RunUMAP(integrated_obj, reduction = "pca", dims = 1:20)
	integrated_obj <- FindClusters(integrated_obj)
	cluster_comp <- table(integrated_obj@meta.data$orig.ident, integrated_obj@meta.data$seurat_clusters)
	pdf(paste(outname, integration_method, "UMAP_LocalPearson_Integrated.pdf", sep="_"), width=5, height=6.5)
	b_loc <- barplot(cluster_comp/ncol(integrated_obj), col=RColorBrewer::brewer.pal(nrow(cluster_comp), "Set3"))
	simpson <- round(simpson_index(integrated_obj@meta.data$orig.ident, integrated_obj@meta.data$seurat_clusters), digits=2)
	text(x=b_loc, y=colSums(cluster_comp/ncol(integrated_obj)), labels=simpson, pos=3)
	legend("topright", c(
			paste("Avg. Simpson =", round(mean(simpson), digits=2)),
			paste("Max. Simpson =", max(simpson))
		), bty="n")
	print(DimPlot(integrated_obj, reduction="umap", group.by="orig.ident"))
	dev.off()
	if (integration_method == "cca") {
		saveRDS(integrated_obj, paste(outname, integration_method, "UMAP_LocalPearson_Integrated.rds", sep="_"))
	}
}

## Global Pearson ##
merged_obj1@assays$SCT@data <- as(global_pearson, "dgCMatrix")
merged_obj1@assays$SCT@scale.data <- as.matrix(global_pearson)
obj_list2 <- SplitObject(merged_obj1, split.by="orig.ident") 
for (integration_method in c("cca", "rpca")) {
	set.seed(897)
	if (integration_method == "rpca") {
		for (i in 1:length(obj_list2)) {
			obj_list2[[i]] <- RunPCA(obj_list2[[i]], assay="SCT", features=consensus_features)
		}
	}
	anchors <- FindIntegrationAnchors(obj_list2, assay=rep("SCT", length(obj_list)), 
			anchor.features = consensus_features, reduction = integration_method)

	integrated_obj <- IntegrateData(anchorset=anchors, features=consensus_features, normalization.method="SCT")
	DefaultAssay(integrated_obj) <- "integrated"
	integrated_obj <- RunPCA(integrated_obj, verbose = FALSE, features=consensus_features)
	integrated_obj <- FindNeighbors(integrated_obj, reduction = "pca", dims = 1:20)
	integrated_obj <- RunUMAP(integrated_obj, reduction = "pca", dims = 1:20)
	integrated_obj <- FindClusters(integrated_obj)
	cluster_comp <- table(integrated_obj@meta.data$orig.ident, integrated_obj@meta.data$seurat_clusters)
	pdf(paste(outname, integration_method, "UMAP_GlobalPearson_Integrated.pdf", sep="_"), width=5, height=6.5)
	b_loc <- barplot(cluster_comp/ncol(integrated_obj), col=RColorBrewer::brewer.pal(nrow(cluster_comp), "Set3"))
	simpson <- round(simpson_index(integrated_obj@meta.data$orig.ident, integrated_obj@meta.data$seurat_clusters), digits=2)
	text(x=b_loc, y=colSums(cluster_comp/ncol(integrated_obj)), labels=simpson, pos=3)
	legend("topright", c(
			paste("Avg. Simpson =", round(mean(simpson), digits=2)),
			paste("Max. Simpson =", max(simpson))
		), bty="n")
	print(DimPlot(integrated_obj, reduction="umap", group.by="orig.ident"))
	dev.off()
	if (integration_method == "cca") {
		saveRDS(integrated_obj, paste(outname, integration_method, "UMAP_GlobalPearson_Integrated.rds", sep="_"))
	}
}
