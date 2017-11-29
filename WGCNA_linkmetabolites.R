## WGCNA_linkmetabolites.R ##
## build the link of GCN modules to metabolites ##
## Zhang Fei ##
## 2017-11-29 ##
library(WGCNA)
library(reshape2)

load("./data/metabolic.Rdata")


## eigen import ##
eigen.JS <- read.table("./result/eigengene/MEs.dissTOM.JS.pearson.deepsplit=2MEDissThres=0.1.txt", header = T, sep = "\t", row.names = 1)
eigen.JC <- read.table("./result/eigengene/MEs.dissTOM.JC.pearson.deepsplit=2MEDissThres=0.1.txt", header = T, sep = "\t", row.names = 1)
eigen.HS <- read.table("./result/eigengene/MEs.dissTOM.HS.pearson.deepsplit=2MEDissThres=0.1.txt", header = T, sep = "\t", row.names = 1)
eigen.HC <- read.table("./result/eigengene/MEs.dissTOM.HC.pearson.deepsplit=2MEDissThres=0.1.txt", header = T, sep = "\t", row.names = 1)


## PM extract ##
GC.JS <- GC.ave[c(46: 60), ]
GC.JC <- GC.ave[c(31: 45), ]
GC.HS <- GC.ave[c(16: 30), ]
GC.HC <- GC.ave[c(1: 15), ]

rownames(GC.JS) <- GC.JS[c(1: 15), 1]
rownames(GC.JC) <- GC.JC[c(1: 15), 1]
rownames(GC.HS) <- GC.HS[c(1: 15), 1]
rownames(GC.HC) <- GC.HC[c(1: 15), 1]

GC.JS <- GC.JS[ , -1]
GC.JC <- GC.JC[ , -1]
GC.HS <- GC.HS[ , -1]
GC.HC <- GC.HC[ , -1]

## Link Calculate ##
moduleTraitCor.JS <- cor(eigen.JS, GC.JS, use = "p")
moduleTraitCor.JC <- cor(eigen.JC, GC.JC, use = "p")
moduleTraitCor.HS <- cor(eigen.HS, GC.HS, use = "p")
moduleTraitCor.HC <- cor(eigen.HC, GC.HC, use = "p")

## data melt ##
melt.JS <- melt(moduleTraitCor.JS)
melt.JC <- melt(moduleTraitCor.JC)
melt.HS <- melt(moduleTraitCor.HS)
melt.HC <- melt(moduleTraitCor.HC)

write.table(melt.JS, "./result/cor_metabolites/PM_JS.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(melt.JC, "./result/cor_metabolites/PM_JC.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(melt.HS, "./result/cor_metabolites/PM_HS.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(melt.HC, "./result/cor_metabolites/PM_HC.txt", quote = F, sep = "\t", row.names = F, col.names = F)