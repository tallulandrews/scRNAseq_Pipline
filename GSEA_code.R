#source: https://bioinformaticsbreakdown.com/how-to-gsea/

require(fgsea)
require(ggplot2)


read_gmt <- function(gmt_file) {
  require(fgsea)
  geneset <- gmtPathways(gmt_file);
  return(geneset)
}


my_collapsePathways <- function(fgseaRes,
                             pathways,
                             stats,
                             pval.threshold=0.05,
                             nperm=10/pval.threshold,
                             gseaParam=1) {
    require(fastmatch)
    universe <- names(stats)

    pathways <- pathways[fgseaRes$pathway]
    pathways <- lapply(pathways, intersect, universe)

    parentPathways <- setNames(rep(NA, length(pathways)), names(pathways))

    for (i in seq_along(pathways)) {
        p <- names(pathways)[i]
        if (!is.na(parentPathways[p])) {
            next
        }

        pathwaysToCheck <- setdiff(names(which(is.na(parentPathways))), p)
        pathwaysUp <- fgseaRes[fgseaRes$pathway %fin% pathwaysToCheck & fgseaRes$ES >= 0, "pathway"]
        pathwaysDown <- fgseaRes[fgseaRes$pathway %fin% pathwaysToCheck & fgseaRes$ES < 0, "pathway"]

        if (length(pathwaysToCheck) == 0) {
            break
        }

        minPval <- setNames(rep(1, length(pathwaysToCheck)), pathwaysToCheck)

        u1 <- setdiff(universe, pathways[[p]])

        fgseaResUp1 <- fgseaSimple(pathways = pathways[pathwaysUp], stats=stats[u1],
                                   nperm=nperm, maxSize=length(u1)-1, nproc=1,
                                   gseaParam=gseaParam, scoreType = "pos")
        fgseaResDown1 <- fgseaSimple(pathways = pathways[pathwaysDown], stats=stats[u1],
                                     nperm=nperm, maxSize=length(u1)-1, nproc=1,
                                     gseaParam=gseaParam, scoreType = "neg")
        fgseaRes1 <- rbindlist(list(fgseaResUp1, fgseaResDown1), use.names = TRUE)

        minPval[fgseaRes1$pathway] <- pmin(minPval[fgseaRes1$pathway], fgseaRes1$pval)

        u2 <- pathways[[p]]

        fgseaResUp2 <- fgseaSimple(pathways = pathways[pathwaysUp], stats=stats[u2],
                                   nperm=nperm, maxSize=length(u2)-1, nproc=1,
                                   gseaParam=gseaParam, scoreType = "pos")
        fgseaResDown2 <- fgseaSimple(pathways = pathways[pathwaysDown], stats=stats[u2],
                                     nperm=nperm, maxSize=length(u2)-1, nproc=1,
                                     gseaParam=gseaParam, scoreType = "neg")
        fgseaRes2 <- rbindlist(list(fgseaResUp2, fgseaResDown2), use.names = TRUE)

        minPval[fgseaRes2$pathway] <- pmin(minPval[fgseaRes2$pathway], fgseaRes2$pval)

        parentPathways[names(which(minPval > pval.threshold))] <- p
    }

    return(list(mainPathways=names(which(is.na(parentPathways))),
                parentPathways=parentPathways))
}


	

GSEA = function(gene_list, gene_set_list=read_gmt("/cluster/projects/macparland/TA/ExternalData/GeneSets/MSigdb_Hallmark_17Aug2020.gmt"), pval=0.05) {
  set.seed(54321)
  require(dplyr)
  require(gage)
  require(fgsea)
  require(data.table)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
 # myGO = fgsea::gmtPathways(GO_file)
  myGO = gene_set_list
  
  # Run GSEA
  fgRes <- fgsea::fgsea(pathways = myGO, 
                           stats = gene_list,
                           minSize=15,
                           maxSize=1000,
                           nperm=10000) %>% 
                  as.data.frame() %>% 
                  dplyr::filter(padj < !!pval)
  #print(dim(fgRes))
    
## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
  gaRes = gage::gage(gene_list, gsets=myGO, same.dir=TRUE, set.size =c(15,1000))
  
  ups = as.data.frame(gaRes$greater) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  downs = as.data.frame(gaRes$less) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  #print(dim(rbind(ups,downs)))
  keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
  keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]
  
  ### Collapse redundant pathways
  #Up = fgsea::collapsePathways(keepups, pathways = myGO, stats = gene_list,  nperm = 500, pval.threshold = 0.05)
  #Down = fgsea::collapsePathways(keepdowns, myGO, gene_list,  nperm = 500, pval.threshold = 0.05) 
  
  Up = my_collapsePathways(keepups, pathways = myGO, stats = gene_list,  nperm = 500, pval.threshold = 0.05)
  Down = my_collapsePathways(keepdowns, myGO, gene_list,  nperm = 500, pval.threshold = 0.05)

  fgRes = fgRes[ !is.na(match(fgRes$pathway, 
           c( Up$mainPathways, Down$mainPathways))), ] %>% 
    arrange(desc(NES))
  fgRes$pathway = stringr::str_replace(fgRes$pathway, "GO_" , "")
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
  geom_point( size=5, aes( fill = Enrichment),
              shape=21, stroke=2) +
    scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
                      "Up-regulated" = "firebrick") ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="GSEA - Biological Processes") + 
    theme_minimal()
  
  output = list("Results" = fgRes, "Plot" = g)
  return(output)
}
