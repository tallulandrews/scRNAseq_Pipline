library(Seurat)
library(CellChat)
source("C:/Users/tandrews/Documents/UHNSonya/scripts/LiverMap2.0/My_R_Scripts.R")

source("Cell_cell_interactions_colourscheme.R")

## CellChat Analysis

meta <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/NKT_fullmetadata.rds")
NKT_obj <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/AllGenes/NKT_harmony_Subcluster_Allgenes.rds")
NKT_obj@meta.data <- meta
NKT_obj@meta.data$Subset <- "NKT"

meta <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/Macrophage_fullmetadata.rds")
Mac_obj <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/AllGenes/Macrophage_harmony_Subcluster_Allgenes.rds")
Mac_obj@meta.data <- meta
Mac_obj@meta.data$Subset <- "Macrophage"

meta <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/Endo_fullmetadata.rds")
Endo_obj <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/AllGenes/Endo_harmony_Subcluster_Allgenes.rds")
Endo_obj@meta.data <- meta
Endo_obj@meta.data$Subset <- "Endothelial"

meta <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/Stellate_fullmetadata.rds")
Stellate_obj <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/AllGenes/Stellate_harmony_Subcluster_Allgenes.rds")
Stellate_obj@meta.data <- meta
Stellate_obj@meta.data$Subset <- "Stellate"

## Create CellChat Object

MergedSeuratObj <- merge(Mac_obj, Endo_obj, add.cell.ids=c("Mac", "Endo"))
MergedSeuratObj <- merge(MergedSeuratObj, Stellate_obj, add.cell.ids=c("", "Ste"))
MergedSeuratObj <- merge(MergedSeuratObj, NKT_obj, add.cell.ids=c("", "NKT"))

rownames(MergedSeuratObj@meta.data) <- colnames(MergedSeuratObj@assays$RNA@data)

# Clean it up - remove genes that are highest in contamination, remove contamination cells.
contam_types <- c("Debris", "HepContam", "Contam", "Flush", "Doublet")
mean_type_expr <- group_rowmeans(MergedSeuratObj@assays$RNA@data, MergedSeuratObj@meta.data$Subcluster_Manual)

hep_contam <- which(colnames(mean_type_expr) %in% c("Debris", "HepContam", "Contam") )

exclude <- apply(mean_type_expr, 1, function(x){
						out<-which(x==max(x)); 
						if (length(out) > 1) {return(hep_contam[1])}
						else {return(out)}})
CleanedMergedSeurat <- MergedSeuratObj[!exclude, ! MergedSeuratObj@meta.data$Subcluster_Manual %in% contam_types]

cellchat_obj <- createCellChat(object=CleanedMergedSeurat, meta=CleanedMergedSeurat@meta.data, group.by="Subcluster_Manual")
cellchat_obj <- setIdent(cellchat_obj, ident.use = "Subcluster_Manual") 
groupSize <- as.numeric(table(cellchat_obj@idents)) 

# Load Interactors
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
#CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor")
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
cellchat_obj@DB <- CellChatDB.use

# DE
cellchat <- subsetData(cellchat_obj)
future::plan("multiprocess", workers = 2) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

# Cell-cell communication
cellchat <- computeCommunProb(cellchat, type="truncatedMean", trim=0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat) # Pathway
cellchat <- aggregateNet(cellchat) # Aggregate by cell-type

# Understand results
df.net <- subsetCommunication(cellchat)


# Visualization
groupSize <- as.numeric(table(cellchat@idents))
#par(mfrow = c(1,2), xpd=TRUE)
#netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
#netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
heatmap(cellchat@net$count, scale="none")
png("CellChat_Weight_Cell-Cell_heatmap.png", width=8, height=8, units="in", res=300)
heatmap(cellchat@net$weight, scale="none")
dev.off()

# Show interactions
table(cellchat@meta$Subcluster_Manual)
png("CellChat_Endo_vs_Mac_bubble_plot.png", width=6, height=6, units="in", res=300)
netVisual_bubble(cellchat, sources.use = c(7,8, 22, 26), targets.use = c(12, 20), remove.isolate = FALSE)
dev.off()
png("CellChat_Endo_vs_NKT_bubble_plot.png", width=6, height=6, units="in", res=300)
netVisual_bubble(cellchat, sources.use =  c(7,8, 22, 26), targets.use = c(6,16), remove.isolate = FALSE)
dev.off()

png("CellChat_Endo_vs_Mac_chord_plot.png", width=6, height=6, units="in", res=300)
netVisual_chord_gene(cellchat, sources.use = c(7,8, 22, 26), targets.use = c(12, 20), lab.cex = 0.5,legend.pos.y = 30)
dev.off()
png("CellChat_Endo_vs_NKT_chord_plot.png", width=6, height=6, units="in", res=300)
netVisual_chord_gene(cellchat, sources.use = c(7,8, 22, 26), targets.use = c(6, 16), lab.cex = 0.5,legend.pos.y = 30)
dev.off()


pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout="hierarchy")

