# scRNAseq_Pipline
A full set of scripts for professional quality analysis of multiple-condition single-cell RNAseq data. 

# Data Compatibility
We assume 10X Chromium data for single-cell / single-nucleus RNAseq. The pipeline is compatible with either 3' or 5' data, but assumes all data has generated with the same technology. In addition, each biological sample should be a separate run (i.e. no pooling of samples).

## Pooled samples
Pooled samples require additional upstream analysis to demultiplex. In addition, doublet identification can take advantage of the nature of pooled samples and we recommend using SouporCell to complement DoubletFinder. Once demultiplexed these data are suitable for using in this pipeline.


## Mixed platforms
If your data contains a mixture of 3' and 5' or sc and sn RNAseq then we recommend adjusting the normalization and data integration steps to scale data for each platform separately before integration. In addition, differential expression should be performed for each technology independently, e.g. healthy 5' vs disease 5' and healthy 3' vs disease 3'.

# Pipeline Steps & Current Script Status

1. QC - EmptyDrops + manual quality metric thresholds
2. Normalization/Scaling - scTransform, LogCPM, scran, M3Drop Pearson Residuals
3. Doublet detection - 
4. Integration - CCA, RPCA, harmony, 
5. Clustering - Seurat + metrics for determining appropriate resolution
6. Annotation - Manual + Reference-based
7. Differential Expression
8. Pathway Analysis
9. Cell-cell interactions
10. Gene Regulatory Network Analysis

# Packages Used / Installation
DropletUtils
Seurat
ggplot2
DoubletFinder
edgeR
RColorBrewer
CellChat
M3Drop
scran
scater
SingleCellExperiment
harmony
scmap
pProfiler
fgsea

# Citation
If you use the scripts in this repository for your analysis we request you cite one of the papers below:

Andrews, T.S., Atif, J., Liu, J.C., Perciani, C.T., Ma, X.-Z., Thoeni, C., Slyper, M., Eraslan, G., Segerstolpe, A., Manuel, J., Chung, S., Winter, E., Cirlan, I., Khuu, N., Fischer, S., Rozenblatt-Rosen, O., Regev, A., McGilvray, I.D., Bader, G.D. and MacParland, S.A. (2022), Single-Cell, Single-Nucleus, and Spatial RNA Sequencing of the Human Liver Identifies Cholangiocyte and Mesenchymal Heterogeneity. Hepatol Commun, 6: 821-840. https://doi.org/10.1002/hep4.1854

Andrews, T.S., Kiselev, V.Y., McCarthy, D. et al. Tutorial: guidelines for the computational analysis of single-cell RNA sequencing data. Nat Protoc 16, 1–9 (2021). https://doi.org/10.1038/s41596-020-00409-w

Clarke, Z.A., Andrews, T.S., Atif, J. et al. Tutorial: guidelines for annotating single-cell transcriptomic maps using automated and manual methods. Nat Protoc 16, 2749–2764 (2021). https://doi.org/10.1038/s41596-021-00534-0

## Other reference materials

Tallulah S. Andrews, Martin Hemberg, Identifying cell populations with scRNASeq, Molecular Aspects of Medicine, Volume 59, 2018, Pages 114-122, ISSN 0098-2997,
https://doi.org/10.1016/j.mam.2017.07.002.

Lun, A., Riesenfeld, S., Andrews, T. et al. EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. Genome Biol 20, 63 (2019). https://doi.org/10.1186/s13059-019-1662-y

Andrews TS, Hemberg M. False signals induced by single-cell imputation. F1000Res. 2018 Nov 2;7:1740. doi: 10.12688/f1000research.16613.2. PMID: 30906525; PMCID: PMC6415334.

Kiselev, V.Y., Andrews, T.S. & Hemberg, M. Challenges in unsupervised clustering of single-cell RNA-seq data. Nat Rev Genet 20, 273–282 (2019). https://doi.org/10.1038/s41576-018-0088-9
