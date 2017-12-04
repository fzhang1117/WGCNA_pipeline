## WGCNA_demo.R ##
## this script is used to examine functions ##
## Zhang Fei ##
## 2017-12-02 ##

library(WGCNA)

load("./data/ave_expression.Rdata")

color.JS <- read.table("./result/modules/merged_modes.dissTOM.JS.pearson.deepsplit=2.MEDissThres=0.10.txt", header = F, sep = "\t")
mergeColors = color.JS$V3

MEList.new <- moduleEigengenes(t(ave.JS), colors = mergeColors, nPC = 1, excludeGrey = T)
