library(nichenetr)
library(Seurat)
library(dplyr)
library(tidyverse)
source("C:/Users/tandrews/Documents/UHNSonya/scripts/LiverMap2.0/My_R_Scripts.R")

source("Cell_cell_interactions_colourscheme.R")

#obj <- readRDS("C37_full_importedClusters_forNicheNet.rds")

#ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#ligand_target_matrix[1:5,1:5]
#saveRDS(ligand_target_matrix, "NicheNet_ligand_target_matrix.rds")

#lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#head(lr_network)

#saveRDS(lr_network, "NicheNet_ligand_receptor_matrix.rds")

#weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
#weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
#head(weighted_networks$lr_sig)
#head(weighted_networks$gr)

#saveRDS(weighted_networks, "NicheNet_prior_wieghts.rds")



## Jan 2022 - Proper NicheNet Analysis
ligand_target_matrix <- readRDS("NicheNet_ligand_target_matrix.rds")
lr_network <- readRDS("NicheNet_ligand_receptor_matrix.rds")
weighted_networks <- readRDS("NicheNet_prior_wieghts.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

meta <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/NKT_fullmetadata.rds")
Mac_obj <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/AllGenes/NKT_harmony_Subcluster_Allgenes.rds")
Mac_obj@meta.data <- meta
Mac_obj@meta.data$Subset <- "NKT"
Reciever_Cell_Type = "NKT"

meta <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/Endo_fullmetadata.rds")
Endo_obj <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/AllGenes/Endo_harmony_Subcluster_Allgenes.rds")
Endo_obj@meta.data <- meta
Endo_obj@meta.data$Subset <- "Endothelial"

meta <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/Stellate_fullmetadata.rds")
Stellate_obj <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/AllGenes/Stellate_harmony_Subcluster_Allgenes.rds")
Stellate_obj@meta.data <- meta
Stellate_obj@meta.data$Subset <- "Stellate"

MergedSeuratObj <- merge(Mac_obj, Endo_obj, add.cell.ids=c("Mac", "Endo"))
MergedSeuratObj <- merge(MergedSeuratObj, Stellate_obj, add.cell.ids=c("", "Ste"))
MergedSeuratObj@meta.data$phony_condition <- MergedSeuratObj@meta.data$Subcluster_Manual
#MergedSeuratObj@meta.data$phony_condition <- "None"
#MergedSeuratObj@meta.data$phony_condition[MergedSeuratObj@meta.data$Subcluster_Manual %in% c("Inflam", "InflaSynap")] <- "Inflammatory"
#MergedSeuratObj@meta.data$phony_condition[MergedSeuratObj@meta.data$Subcluster_Manual %in% c("NonInf", "PhagoNonInf", "ResNonInf")] <- "NonInflammatory"
#MergedSeuratObj@meta.data$phony_condition[MergedSeuratObj@meta.data$Subcluster_Manual %in% c("InfActiv", "Activated")] <- "Activated"

MergedSeuratObj@active.ident <- factor(MergedSeuratObj@meta.data$Subcluster_Manual, levels=c(Reciever_Cell_Type, sort(unique(MergedSeuratObj@meta.data$Subcluster_Manual))))
MergedSeuratObj@active.ident[MergedSeuratObj@meta.data$Subset == Reciever_Cell_Type] <- Reciever_Cell_Type
rownames(MergedSeuratObj@meta.data) <- colnames(MergedSeuratObj@assays$RNA@data)
# To do: make a "best of" hepatocyte object

## receiver
receiver = Reciever_Cell_Type
condition_io ="cNKcell"
condition_ref="lrNKcell"

expressed_genes_receiver = rownames(MergedSeuratObj)[Matrix::rowMeans(MergedSeuratObj@assays$RNA@counts[,MergedSeuratObj@meta.data$Subset == Reciever_Cell_Type] > 0) > 0.1]
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
#sender_celltypes = c("Fibro-aHSC", "Fibro-qHSC", "cvEndo", "cvLSEC", "ppLSEC", "VasEndo")
#sender_colours <- c(cell_cell_colours[["Stellate"]], cell_cell_colours[["Stellate"]], 
#				cell_cell_colours[["Endo"]], cell_cell_colours[["cvLSEC"]], 
#				cell_cell_colours[["ppLSEC"]],cell_cell_colours[["Endo"]])

sender_celltypes = c("cvEndo", "cvLSEC", "ppLSEC", "VasEndo")
sender_colours <- c( cell_cell_colours[["Endo"]], cell_cell_colours[["cvLSEC"]], 
			   cell_cell_colours[["ppLSEC"]],"yellow")

names(sender_colours) <- sender_celltypes
 
list_expressed_genes_sender = unique(sapply(sender_celltypes, function(x) {rownames(MergedSeuratObj)[Matrix::rowMeans(MergedSeuratObj@assays$RNA@counts[,MergedSeuratObj@meta.data$Subcluster_Manual == x] > 0) > 0.1]} )) 
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

# genesets


DE_table_receiver <- FindMarkers(object = MergedSeuratObj, ident.1 = condition_io, ident.2 = condition_ref, group.by="phony_condition", min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

# Active Ligands
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
				background_expressed_genes = background_expressed_genes, 
				ligand_target_matrix = ligand_target_matrix, 
				potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities

best_upstream_ligands = ligand_activities$test_ligand[ligand_activities$pearson > 0.03] 


#best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
#MergedSeuratObj<-RunUMAP(MergedSeuratObj, features=unique(c(VariableFeatures(Mac_obj), VariableFeatures(Endo_obj), VariableFeatures(Stellate_obj))))

DotPlot(Endo_obj, features = rev(best_upstream_ligands), group.by="Subcluster_Manual", cols = "RdYlBu") + RotatedAxis()
DotPlot(Stellate_obj, features = rev(best_upstream_ligands), group.by="Subcluster_Manual", cols = "RdYlBu") + RotatedAxis()


# Active Ligand Targets
#active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na() # macrophage thresholds
#active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33) # macrophage thresholds

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 5000) %>% bind_rows() %>% drop_na() # NKT thresholds
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33) # NKT thresholds


order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

# Receptors of the ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
    
dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

# Constrained to proper interactions
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict

# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
#DE_table_all = lapply(sender_celltypes, get_lfc_celltype, seurat_obj = MergedSeuratObj, condition_colname = "phony_condition", 
#								condition_oi = condition_io, condition_reference = condition_ref,
#								expression_pct = 0.10, celltype_col = "Subcluster_Manual") % reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
#DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information


ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_receiver %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

#nichenet_output <- nichenet_seuratobj_aggregate(seurat_obj=MergedSeuratObj, receiver="Macrophage", sender=sender,
#	condition_colname="phony_condition", condition_oi="NonInflammatory", condition_reference="Inflammatory",
#	ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, 
#	organism = "human")


###### Summary Figures ######
#Dotplot of expression in senders

ligand_2_target <- vis_ligand_target

anno_targets<-function(ligand, n=5, l2t=ligand_2_target, de=DE_table_receiver) {
	if(! ligand %in% rownames(l2t)) {warning(paste(ligand, "not in ligand2target matrix")); return()}
	targets <- tail(sort(l2t[rownames(l2t) == ligand,]), n)
	targets <- names(targets[targets > 0]) # account for those with fewer than n targets.
	if (length(targets) == 0) {warning(paste(ligand, "has no targets")); return()}
	dir_targets = de[match(targets, de$gene),];
	dir_targets = sign(dir_targets$avg_log2FC);
	dir_targets[is.na(dir_targets)] <- 0;
	return(data.frame(gene=targets, direction=dir_targets)) # +ve Condition_io, -ve Condition_ref
}


# heat map of ligand 2 receptors.
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

these_ligands <- rownames(ligand_2_target)
these_receptors <- rownames(vis_ligand_receptor_network)

ligand_2_receptors <- lr_network %>% filter(from %in% these_ligands & to %in% these_receptors) %>% distinct(from,to)
heatmap(  (vis_ligand_receptor_network>0.3)+1, scale="none")

tmp <- merge(Stellate_obj, Endo_obj)
tmp <- tmp[,tmp@meta.data$Subcluster_Manual %in% sender_celltypes]
DotPlot(tmp, features = these_ligands, group.by="Subcluster_Manual", cols = "RdYlBu")  + coord_flip() +scale_y_discrete(position = "right") 
DotPlot(Mac_obj, features = these_receptors, group.by="Subcluster_Manual", cols = "RdYlBu")  + coord_flip() +scale_y_discrete(position = "right") 



#heatmap of targets genes - colour bar for DE in Inflam / NonInflam
hist(vis_ligand_target, breaks=20)
vis_ligand_target_clean <- vis_ligand_target
#vis_ligand_target_clean[vis_ligand_target_clean < 0.003] <- 0 # Macrophage Threshold
vis_ligand_target_clean[vis_ligand_target_clean < 0.00135] <- 0 # NKT Threshold
p_ligand_target_network = vis_ligand_target_clean %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

# Ligand 2 target
ligand_2_targets <- vis_ligand_target_clean[rowSums(vis_ligand_target_clean) > 0, colSums(vis_ligand_target_clean) > 0]
DotPlot(Mac_obj, features = sub("\\.","-", colnames(ligand_2_targets)), group.by="Subcluster_Manual", cols = "RdYlBu")  + coord_flip() +scale_y_discrete(position = "right") 
mac_means <- group_rowmeans(Mac_obj@assays$RNA@data, Mac_obj@meta.data$Subcluster_Manual)



# get rid of hepatocyte / ambient stuff
hep_contam <- which(colnames(mac_means) %in% c("Debris", "HepContam") )
exclude <- apply(mac_means, 1, function(x){
						out<-which(x==max(x)); 
						if (length(out) > 1) {return(hep_contam)}
						else {return(out)}}
			)
exclude <- rownames(mac_means)[exclude == hep_contam]
DotPlot(Mac_obj, features = colnames(ligand_2_targets)[! colnames(ligand_2_targets) %in% exclude], 
		group.by="Subcluster_Manual", cols = "RdYlBu")  + coord_flip() +scale_y_discrete(position = "right") 

ligand_2_targets<-ligand_2_targets[,! colnames(ligand_2_targets) %in% exclude]
ligand_2_targets <- ligand_2_targets[rowSums(ligand_2_targets) > 0,]

Links <- which(ligand_2_targets > 0, arr.ind=T)
Links <- data.frame(ligand=rownames(ligand_2_targets)[Links[,1]],
				target=paste(" ", colnames(ligand_2_targets)[Links[,2]], sep=""),
				weight=ligand_2_targets[which(ligand_2_targets>0)])

#target_colour_scheme <- RColorBrewer::brewer.pal(8, "RdBu")
target_colour_scheme <- PurpleAndYellow(12); target_colour_scheme<-target_colour_scheme[-1*seq(length(target_colour_scheme)/2-1, length(target_colour_scheme)/2+2)]
target_colour <- DE_table_receiver[match(sub("\\.", "-", colnames(ligand_2_targets)), DE_table_receiver$gene), "avg_log2FC"]
limit <- max(abs(target_colour))*1.0001
target_colour <- target_colour_scheme[cut(target_colour, breaks=seq(from=-1*limit, to=limit, length=length(target_colour_scheme)+1))]
names(target_colour) <- colnames(ligand_2_targets)

# Assign each ligand to the cell-type with highest mean expression
sender_means <- group_rowmeans(MergedSeuratObj@assays$RNA@data, MergedSeuratObj@meta.data$Subcluster_Manual)
sender_detect <- group_rowmeans(MergedSeuratObj@assays$RNA@data>0, MergedSeuratObj@meta.data$Subcluster_Manual)
sender_means <- sender_means[,colnames(sender_means) %in% sender_celltypes]

max_exp_type <- apply(sender_means[match(sub("\\.", "-", rownames(ligand_2_targets)), rownames(sender_means)),], 1, function(x) {
							colnames(sender_means)[which(x == max(x))]})
max_exp_val <-  apply(sender_detect[match(sub("\\.", "-", rownames(ligand_2_targets)), rownames(sender_means)),], 1, function(x) {
							x[which(x == max(x))]})

		
ligand_2_celltype <- data.frame(names(max_exp_type), max_exp_type, colour=sender_colours[max_exp_type])
ligand_2_celltype <- ligand_2_celltype[order(ligand_2_celltype[,3], ligand_2_celltype[,2]),]
ligand_colour <- ligand_2_celltype$colour; names(ligand_colour) <- rownames(ligand_2_celltype)

# Circos visualization:
# Sender ligand genes = colour by source cell-type
# Reciever target genes = colour by DE 

sorted_target_colour <- target_colour[order(as.numeric(factor(target_colour, levels=target_colour_scheme)))]

all_genes_in_order <- c(paste(" ",names(sorted_target_colour), sep=""), ligand_2_celltype[,1])
all_genes_in_order_col <- c(sorted_target_colour, ligand_2_celltype[,3]);
names(all_genes_in_order_col) <- all_genes_in_order

width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

Links$ligand_type <- ligand_2_celltype[match(sub("\\.", "-", Links[,1]), sub("\\.", "-", ligand_2_celltype[,1])),2]#cell type expressing the ligand for each link
Links$colour <- ligand_2_celltype[match(sub("\\.", "-", Links[,1]), sub("\\.", "-", ligand_2_celltype[,1])),3]
Links <- Links[Links[,1] %in% all_genes_in_order & Links[,2] %in% all_genes_in_order,]


gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = sum(all_genes_in_order_col %in% sorted_target_colour)-1),
  width_different_cell
)

for (type in unique(ligand_2_celltype[,2])) {
  gaps <- c(gaps, 
       rep(width_same_cell_same_ligand_type, times = sum(ligand_2_celltype[,2] == type)-1),
       width_different_cell
  )
}

png("NKT_nichenetr_circos2.png", width=8, height=8, units="in", res=300)
require(circlize)

circos.clear()
#circos.par(gap.degree = gaps[1:(length(gaps)-1)])
circlize::chordDiagram(Links, directional=1, direction.type = c("diffHeight", "arrows"),
				order=names(all_genes_in_order_col), link.sort = TRUE, link.decreasing = FALSE,
				transparency = 0, grid.col=all_genes_in_order_col,
				diffHeight = 0.005, 
				annotationTrack="grid", preAllocateTracks = list(track.height = 0.075),
				link.arr.type = "big.arrow", big.gap=0.5, small.gap=0.5				
				)
circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) 
dev.off()

# Legends
png("NKT_nichenetr_circos2_legend.png", width=5, height=5, units="in", res=300)
source("C:/Users/tandrews/Documents/UHNSonya/scripts/LiverMap2.0/My_R_Scripts.R")

plot(1,1, col="white", xaxt="none", yaxt="none", bty="none", xlab="", ylab="", main="")
legend("topleft", bty="n", fill=sender_colours, names(sender_colours))

colour_gradiant_scale <- signif(seq(from=-1*limit, to=limit, length=length(target_colour_scheme)+1), digits=1)
legend("bottomleft",
       legend = colour_gradiant_scale[-1*ceiling(length(colour_gradiant_scale)/2)],
       col = target_colour_scheme,
       border = NA, y.intersp = 0.75, pt.cex=2,
       cex = 1, bty="n", pch=15, title="L2FC")
dev.off()

# Igraph
Links <- which(ligand_2_targets > 0, arr.ind=T)
Links <- data.frame(ligand=rownames(ligand_2_targets)[Links[,1]],
				target=colnames(ligand_2_targets)[Links[,2]],
				weight=ligand_2_targets[which(ligand_2_targets>0)])


node_colours <- c(target_colour, ligand_colour)
require(igraph)
G <- graph_from_edgelist(cbind(Links[,1], Links[,2]))
E(G)$Weight <- Links[,3]
V(G)$Colour <- node_colours[match(sub("\\.", "-", names(V(G))), sub("\\.", "-", names(node_colours)))]
plot(G, vertex.color=V(G)$Colour, edge.width=E(G)$Weight*500, edge.arrow.size=E(G)$Weight*200)


require(ComplexHeatmap)

col_fun <- function(x){circlize::colorRamp2(c(0, max(vis_ligand_target)), c("white", "purple"))}
Heatmap(vis_ligand_target, color_space=col_fun)













# My Analysis
ligand_target_matrix <- readRDS("NicheNet_ligand_target_matrix.rds")


meta <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/Macrophage_fullmetadata.rds")
obj <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/AllGenes/Macrophage_harmony_Subcluster_Allgenes.rds")
obj@meta.data <- meta

obj <- obj[,! (obj@meta.data$Subcluster_Manual %in% c("Debris", "Doublet"))]

require(Matrix)
is.3pr = obj@meta.data$assay_type == "3pr"
is.5pr = obj@meta.data$assay_type == "5pr"
pseudobulks_3pr <- group_rowmeans(obj@assays$RNA@scale.data[,is.3pr], obj@meta.data[is.3pr,"Subcluster_Manual"])
pseudobulks_5pr <- group_rowmeans(obj@assays$RNA@scale.data[,is.5pr], obj@meta.data[is.5pr,"Subcluster_Manual"])
pseudobulks <- group_rowmeans(obj@assays$RNA@scale.data, obj@meta.data[,"Subcluster_Manual"])

identical(rownames(pseudobulks_3pr), rownames(pseudobulks_5pr))

# Harmonize Genes
ligand_target_matrix <- ligand_target_matrix[rownames(ligand_target_matrix) %in% rownames(pseudobulks),]
ligand_target_matrix <- ligand_target_matrix[,colSums(ligand_target_matrix)>0]
pseudobulks_3pr <- pseudobulks_3pr[match(rownames(ligand_target_matrix), rownames(pseudobulks_3pr)),]
pseudobulks_5pr <- pseudobulks_5pr[match(rownames(ligand_target_matrix), rownames(pseudobulks_5pr)),]
pseudobulks <- pseudobulks[match(rownames(ligand_target_matrix), rownames(pseudobulks)),]


#pseudobulks_3pr <- t(t(pseudobulks_3pr)/colSums(pseudobulks_3pr)*5000)
#pseudobulks_5pr <- t(t(pseudobulks_5pr)/colSums(pseudobulks_5pr)*5000)
#pseudobulks <- t(t(pseudobulks)/colSums(pseudobulks)*5000)

# Pathway Weights - by rescaling the ligand-target-matrix the scores become weighted means where the sum of weights = 1
ligand_target_matrix <- t(t(ligand_target_matrix)/colSums(ligand_target_matrix)) 

scores_3pr <- t(ligand_target_matrix) %*% pseudobulks_3pr
scores_5pr <- t(ligand_target_matrix) %*% pseudobulks_5pr
scores <- t(ligand_target_matrix) %*% pseudobulks


# save result table
write.table(scores_3pr, file="Macrophage_NicheNet_3pr_scores.csv", sep=",")
write.table(scores_5pr, file="Macrophage_NicheNet_5pr_scores.csv", sep=",")
write.table(scores, file="Macrophage_NicheNet_scores.csv", sep=",")

#Norm_Scores <- t( t(scores)/colSums(scores))
#size_factors <- colSums(abs(pseudobulks))
#Norm_Scores <- t( t(scores)/(size_factors/median(size_factors)) )


# Bootstrap confidence? - This is way too slow!
bootstrap_weighted_confidence <- function(vals, weights, n = 10000) {
	score = sum(vals*weights)
	return(sum(rnd_scores > score)/n)
	#return(rnd_score)
}

# Use correlations for significance.
p_matrix <- matrix(-1, nrow=nrow(Norm_Scores), ncol=ncol(Norm_Scores))
r_matrix <- matrix(-1, nrow=nrow(Norm_Scores), ncol=ncol(Norm_Scores))

time = Sys.time()
for (pathway in 1:ncol(ligand_target_matrix)) {
	for (type in 1:ncol(pseudobulks)) {
		out <- cor.test(pseudobulks[,type], ligand_target_matrix[,pathway])
		p_matrix[pathway, type] <- out$p.value
		r_matrix[pathway, type] <- out$estimate

	}
}
colnames(r_matrix) <- colnames(scores)
rownames(r_matrix) <- rownames(scores)

Bonferroni_threshold <- 0.05/prod(dim(p_matrix))

tidied_scores <- scores
tidied_scores[p_matrix > Bonferroni_threshold] <- 0
r_matrix[p_matrix > Bonferroni_threshold] <- 0

#tidied_scores <- t( t(tidied_scores)/size_factors )


hvp <- apply(tidied_scores, 1, var)
ntop = 50
heatmap.2(r_matrix[hvp > quantile(hvp, prob=1-ntop/nrow(tidied_scores)),], scale="none", 
		hclustfun=function(x){hclust(x, method="ward.D")}, trace="none")


tmp <- Norm_Scores
tmp[p_matrix > Bonferroni_threshold] <- 0

# Keep max score from either 3pr or 5pr
tmp_3pr <- scores_3pr
tmp_3pr[scores_3pr < scores_5pr] <- 0
tmp_5pr <- scores_5pr
tmp_5pr[scores_5pr < scores_3pr] <- 0

overall_score <- tmp_3pr+tmp_5pr

head(scores_5pr[order(scores_5pr[,4], decreasing=TRUE),], 20)
head(scores_3pr[order(scores_3pr[,4], decreasing=TRUE),], 20)
head(overall_score[order(overall_score[,4], decreasing=TRUE),], 20)

write.table(overall_score, "NKT_NichNet_Scores.tsv", sep="\t", row.names=T, col.names=T)

