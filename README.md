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
