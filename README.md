# Overview - WGCNA


# Data format

# WGCN Build

## Step 1: Import data

## Step 2: Soft thresholiding power select

## Step 3: Build a co-expression network

A co-expression nework is a m * m matrix record gene-to-gene relationships. R package WGCNA offers the function 'blockwiseModules()' to do automatic network construction and module detection. This function performs automatic network construction and module detection on large expression datasets in a block-wise manner.

### blockwiseModules()
#### important index
- **blocks:** optional specification of blocks in which hierarchical clustering and module detection should be performed. If given, must be a numeric vector with one ectry per column (gene) of exprData giving the number of the block ot which the corresponding gene belongs.
- **maxBlockSize:** maximum block size for module detection.
- **nPreclusteringCenters:** number of centers for pre-clustering. Large numbers typically results in better but slower pre-clustering.
- **loadTOM:** logical: should Topological Overlap Matrices be loaded from previously saved files (TRUE) or calculated (FALSE)? It may be useful to load previously saved TOM matrices if these have been calculated previously, since TOM calculation is often the most computationally expensive part of network construction and module identification. See saveTOM and saveTOMFileBase below for when and how TOM files are saved, and what the file names are. If loadTOM is TRUE but the files cannot be found, or do not contain the correct TOM data, TOM will be recalculated.
- **corType:** character string specifying the correlation to be used. Allowed values are (unique abbreviations of) 'pearson' and 'bicor', corresponding to Pearson and bidweight midcorrelation, respectively. Missing values are handled using the pairwise.complete.obs options.
- **maxPOutliers:** only used for cor Type == "bicor". Specifies the maximum percentile of data that can be considered outliers on either side of the median separately. For each side of the median, if higher percentile than maxPOutliers is considered an outlier by the weight function based on 9 * mad(x), the width of the weight function is increased such that the percentile of outliers on taht side of the median equals maxPOutliers. Using maxPOutliers = 1 will effectively disable all weight function broadning; using maxPOutliers=0 wil give rsults that are quite similar(but not equal to) Pearson correlation.
- **power:** soft-thresholding power for network construction.
- **networkType:** network type. Allowed values are (unique abbreviations of) **"unsigned"**, **"signed"** and **"signed hybrid"**.
- **TOMType:** one of 'none', 'unsigned', 'signed'. If 'none', adjacency will be used for clustering . If 'unsigned', the standard TOM will be used. If 'signed', TOM will keep track of the sign of correlations between nerghbors.
- **saveTOMs:** logical: should the consensus topological overlap matrices for each block be saved and returned?
- **saveTOMFileBase:** character string containing the file name base for files containing the consensus topological overlaps. The full file names have "block.1.RData", "block.2.RData" etc. appenden. These files are R data files and can be loaded using the load function.
- **deepSplit:** integer value between 0 and 4. Provides a simplified control over how sensitive module detection should be to module splitting, with 0 least and 4 most sensitive.
- **detectCutHeight:** dendrogram cut height for module detection.
- **minModuleSize:** minimum module size for module detection.
- **numericLabels:** logical: should returned modules be labeled by colors (FALSE), or by numbers (TRUE)?
#### Details
Before module detection start, genes and samples are optionally checked for the presence of NAs. Genes and/or samples that have too many NAs are flagged as bad and removed from the analysis; bad genes will be automotically labeled as unassigned, while the returned eigengenes will have NA entries for all bad samples.

If blocks is not given and the number of genes exceeds maxBlockSize, genes are pre-clustered into blocks using the function projectiveKMeans; otherwise all genes are treated in a single block.

For each block of genes, the network is constructed and (if requested) topological overlap is calculated. If requested, the topological overlaps are returned as part of the return value list. Genes and then clustered using average linkage hierarchical clustering and modules are indentified in the resulting dendrogram by the Dynamic Hybrid tree cut. Found modules are trimmed of genes whose correlation with module eigengene(KME) is less than minKMEtoStay. Modules in which fewer than minCoreKMESize genes have KME higher than minCoreKME are disbanded, i.e., their constituent genes are pronounced unassigned.

After all blocks have been processed, the function check whether there are genes whose KME in the module they assigned is lower than KME to another module. If p-values of the higher correlations are smaller than those of the native module by the factor reassignThresholdPS, the gene is reassigned to the closer module.

In the last step, modules whose eigengenes are highly correlated are merged. This is achieved by clustering module eigengenes using the dissimilarity given by one minus their correlation, cutting the dendrogram at the height mergeCutHeight and merging all modules on each branch. The process is iterated until no modules are merged.

#### Output
The function blockwiseModules returns a list with the following components:

- **colors:** a vector of color or numeric module labels for all genes
- **unmergedColors:** a vector of color or numeric module labels for all genes before module merging.
- **MEs:** a data frame containing module eigengenes of the found modules (given by colors).
- **goodSamples:** numeric vector giving indices of good samples, that is samples that do not have too many missing entries
- **goodGenes:** numeric vector giving indices of good genes, that is genes that do not have too many missing entries.
- **dendrograms:** a list whose components conatain hierarchical clustering dendrograms of genes in each block.
- **TOMFiles:** if saveTOMs == TRUE, a vector of character strings, one string per block giving the file names of files (relative to current directory) in which blockwise topological overlaps were saved.
- **blockGenes:** a list whose components give the indices of genes in each block.
- **blocks:** if input blocks were given, its copy; otherwise a vector of length equal number of genes giving the block label for each gene. Note that block label. Note that block labels are not necessarilly sorted in the order in which the blocks were processed (since we do not require this for the input blocks).
- **blockOrder:** a vector giving the order in which blocks were processed and in which blockGenes above is returned. For example, blockOrder contains the label of the first-processed block.

- **MEsOK:** logical indicating whether the module eigengenes were calculated without eroors.

#### Note
If the input dataset has a large number of genes, consider carefully the maxBlockSize as it significantly affects the memory footprint (and whether the function will fail with a memory allocation error). From a theoretical point of view it is advantageous to use blocks as large as possible; on the other hand, using smaller blocks is substantially faster and often the only way to work with large numbers of genes. As a rough guide, it is unlikely a standard desktoop computer with 4GB memory or less will be able to work with blocks larger than 8000 genes.

## Step 4: Module visualize

## Step 5: Module merge

## Step 6: Module Enrichment
