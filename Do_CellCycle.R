my_metadata_table <- read.table("~/scripts/LiverMap2.0/Metadata20LiverMapPlusParams.csv", sep=",", header=T, stringsAsFactors=FALSE);

all_cc <- list();


for(dataset_row in 1:25) {

set.seed(my_metadata_table$Seed[dataset_row])
name <- my_metadata_table$Name[dataset_row]
npcs <- my_metadata_table$nPCs[dataset_row]
nkNN <- my_metadata_table$kNN[dataset_row]
res <- 5


print(name);
require(dplyr)
require(Seurat)
require(Matrix)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#scmap based
if (!file.exists(paste(name, "Anno_SeurObj.rds", sep="_"))) {next;}
myseur <- readRDS(paste(name,"Anno_SeurObj.rds", sep="_"));
myseur <- NormalizeData(myseur)
myseur <- FindVariableFeatures(myseur, selection.method = "vst")
myseur <- ScaleData(myseur, features = rownames(myseur))

myseur <- CellCycleScoring(myseur, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#RidgePlot(myseur, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
#myseur <- RunPCA(myseur, features = c(s.genes, g2m.genes))
#DimPlot(myseur)

all_cc[[name]] <- mysur@meta.data

}

saveRDS(all_cc, "All20_cellcycle.rds");





