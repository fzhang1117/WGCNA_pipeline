# Overview - WGCNA


# Data format

# WGCN Build

## Step 1: Import data

## Step 2: Soft thresholiding power select

## Step 3: Build a co-expression network

A co-expression nework is a m * m matrix record gene-to-gene relationships. R package WGCNA offers the function 'blockwiseModules()' to do automatic network construction and module detection. This function performs automatic network construction and module detection on large expression datasets in a block-wise manner.

### blcokwishModules()

maxBlockSize: maximum block size for module detection.

nPreclusteringCenters: number of centers for pre-clustering. Large numbers typically results in better but slower pre-clustering

loadTOM: logical: should Topological Overlap Matrices be loaded from previously saved files (TRUE) or calculated (FALSE)? It may be useful to load previously saved TOM matrices if these have been calculated previously, since TOM calculation is often the most computationally expensive part of network construction and module identification. See saveTOM and saveTOMFileBase below for when and how TOM files are saved, and what the file names are. If loadTOM is TRUE but the files cannot be found, or do not contain the correct TOM data, TOM will be recalculated.

corType: character string specifying the correlation to be used. Allowed values are (unique abbreviations of) 'pearson' and 
bicor', corresponding to Pearson and bidweight midcorrelation, respectively. Missing values are handled using the pairwise.complete.obs options.

maxPOutliers: only used for cor Type == 'bicor'


## Step 4: Module visualize

## Step 5: Module merge

## Step 6: Module Enrichment
