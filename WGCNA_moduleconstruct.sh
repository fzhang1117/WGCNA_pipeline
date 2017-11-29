#! /bin/bash

TOMmatrix=$1
cor_method=$2
minModuleSize=$3
deepSplit=$4
MEDissThres=$5

#Rscript WGCNA_dissTOM.R

## WGCNA clustering
## example:
##
## Rscript WGCNA_TOMclustering.R ./data/dissTOM.pearson.Rdata pearson 30 2 0.10

Rscript WGCNA_TOMclustering.R $TOMmatrix $cor_method $minModuleSize $deepSplit $MEDissThres
mv ./result/*dissTOM*pdf ./result/visualization/
mv ./result/merged_modes*txt ./result/modules/
mv ./result/MEs.dissTOM*txt ./result/eigengene/

Rscript WGCNA_SumModule.R
