do_fgsea <- function(scored_genes, pathways=MSigAll, fdr=0.05, nmax=20, seed=101){
	set.seed(101)
	res <- fgsea(pathways, scored_genes, minSize=15, maxSize=1000, eps=0, nPermSimple=20000)
	res <- res[!is.na(res$pval) & res$padj < fdr,]
	res <- res[order(res$NES),]
	res_full <- res;
	if (nrow(res) > nmax) {
		res_pos <- data.frame(res[unlist(res$NES) >0,])
		res_pos <- res_pos[!is.na(unlist(res_pos[,1])),]
		res_neg <- data.frame(res[unlist(res$NES) <0,])
		res_neg <- res_neg[!is.na(unlist(res_neg[,1])),]
		res_pos <- res_pos[order(abs(unlist(res_pos$NES)), decreasing=TRUE),]
		res_neg <- res_neg[order(abs(unlist(res_neg$NES)), decreasing=TRUE),]
		res <- rbind(res_pos[1:min(nrow(res_pos), nmax),], res_neg[1:min(nrow(res_neg), nmax),])
		res <- res[order(res$NES),]
	}

	size <- abs(res$NES)
	colour <- sign(res$NES)
	col_palette <- c("dodgerblue", "grey50", "firebrick")
	gene_lists <- res[,"leadingEdge"]
	if (!is.null(dim(gene_lists))) {
		gene_lists_new <- list()
		for (i in 1:nrow(gene_lists)) {
			gene_lists_new[[i]] <- gene_lists[i,1]
		}
		gene_lists <- gene_lists_new
	}

	sim_mat <- matrix(0, nrow=length(gene_lists), ncol=length(gene_lists))
	for (i in 1:length(gene_lists)) {
		for (j in i:length(gene_lists)) {
			int <- length(intersect(unlist(gene_lists[i]), unlist(gene_lists[j])))
			uni <- length(union(unlist(gene_lists[i]), unlist(gene_lists[j])))
			sim_mat[i,j] <- int/uni
			sim_mat[j,i] <- int/uni
			colnames(sim_mat) <- unlist(res[,1])
			rownames(sim_mat) <- unlist(res[,1])
		}
	}
	require(igraph)
	G <- simplify(graph_from_adjacency_matrix(sim_mat > 0.1, mode="undirected"))
	plot(G, vertex.color=col_palette[colour+2], vertex.size=size*5, edge.width=2)
	res$cluster <- components(G)$membership
	return(list(rich=res_full, graph=G, vertex_col = col_palette[colour+2], vertex_size = size*5))
}


require(fgsea)
immune_path <- gmtPathways("C:/Users/tandrews/Documents/UHNSonya/ExternalData/MSigDb_immune_signatures_c7.all.v7.1.symbols.gmt")
Hallmark_path <- gmtPathways("C:/Users/tandrews/Documents/UHNSonya/ExternalData/MSigDbHalmarkPathways.gmt")
MSigAll <- gmtPathways("C:/Users/tandrews/Documents/UHNSonya/ExternalData/MSigDb_curated_c2.all.v7.1.symbols.gmt")
reactome <- gmtPathways("C:/Users/tandrews/Documents/UHNSonya/ExternalData/ReactomePathways.gmt")
BaderMSig <- gmtPathways("C:/Users/tandrews/Documents/UHNSonya/ExternalData/BaderLab25Aug2020/Human_MSigdb_August_01_2020_symbol.gmt.txt")
BaderWP <- gmtPathways("C:/Users/tandrews/Documents/UHNSonya/ExternalData/BaderLab25Aug2020/Human_WikiPathways_August_01_2020_symbol.gmt.txt")
BaderKegg <- gmtPathways("C:/Users/tandrews/Documents/UHNSonya/ExternalData/BaderLab25Aug2020/Human_KEGG_August_01_2020_symbol.gmt.txt")

