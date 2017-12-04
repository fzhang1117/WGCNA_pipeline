## WGCNA_linkmetabolites.R ##
## build the link of GCN modules to metabolites ##
## Zhang Fei ##
## 2017-11-29 ##
library(WGCNA)
library(reshape2)
library(pheatmap)

load("./data/metabolic.Rdata")

## function import ##
my.reshape <- function(data){
    data <- dcast(data, Var1 ~ Var2)
    rownames(data) <- data[ , 1]
    data <- data[ , -1]
    return(data)
}


## eigen import AE import ##
# eigen.JS <- read.table("./result/eigengene/MEs.dissTOM.JS.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t", row.names = 1)
# eigen.JC <- read.table("./result/eigengene/MEs.dissTOM.JC.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t", row.names = 1)
# eigen.HS <- read.table("./result/eigengene/MEs.dissTOM.HS.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t", row.names = 1)
# eigen.HC <- read.table("./result/eigengene/MEs.dissTOM.HC.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t", row.names = 1)

AE.JS <- read.table("./result/modules/AEs.dissTOM.JS.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t", row.names = 1)
AE.JC <- read.table("./result/modules/AEs.dissTOM.JC.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t", row.names = 1)
AE.HS <- read.table("./result/modules/AEs.dissTOM.HS.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t", row.names = 1)
AE.HC <- read.table("./result/modules/AEs.dissTOM.HC.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t", row.names = 1)
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
# moduleTraitCor.JS <- cor(eigen.JS, GC.JS, use = "p")
# moduleTraitCor.JC <- cor(eigen.JC, GC.JC, use = "p")
# moduleTraitCor.HS <- cor(eigen.HS, GC.HS, use = "p")
# moduleTraitCor.HC <- cor(eigen.HC, GC.HC, use = "p")
moduleTraitCor.JS <- cor(AE.JS, GC.JS, use = "p")
moduleTraitCor.JC <- cor(AE.JC, GC.JC, use = "p")
moduleTraitCor.HS <- cor(AE.HS, GC.HS, use = "p")
moduleTraitCor.HC <- cor(AE.HC, GC.HC, use = "p")



## data melt ##
melt.JS <- melt(moduleTraitCor.JS)
melt.JC <- melt(moduleTraitCor.JC)
melt.HS <- melt(moduleTraitCor.HS)
melt.HC <- melt(moduleTraitCor.HC)

write.table(melt.JS, "./result/cor_metabolites/PM_JS.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(melt.JC, "./result/cor_metabolites/PM_JC.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(melt.HS, "./result/cor_metabolites/PM_HS.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(melt.HC, "./result/cor_metabolites/PM_HC.txt", quote = F, sep = "\t", row.names = F, col.names = F)

## link select ##
JS.select <- na.omit(melt.JS[melt.JS[ , 3] >= 0.9 | melt.JS[ , 3] <= -0.9, ])
JS.select[ , 1] <- paste('JS', JS.select[ , 1], sep = ".")
JC.select <- na.omit(melt.JC[melt.JC[ , 3] >= 0.9 | melt.JC[ , 3] <= -0.9, ])
JC.select[ , 1] <- paste('JC', JC.select[ , 1], sep = ".")

HS.select <- na.omit(melt.HS[melt.HS[ , 3] >= 0.9 | melt.HS[ , 3] <= -0.9, ])
HS.select[ , 1] <- paste('HS', HS.select[ , 1], sep = ".")
HC.select <- na.omit(melt.HC[melt.HC[ , 3] >= 0.9 | melt.HC[ , 3] <= -0.9, ])
HC.select[ , 1] <- paste('HC', HC.select[ , 1], sep = ".")

JI.select <- rbind(JS.select, JC.select)
HZ.select <- rbind(HS.select, HC.select)

link.select <- rbind(JI.select, HZ.select)
link.select[ , 2] <- as.character(link.select[ , 2])

## fill these not significant 
cpd.select <- sort(link.select[ , 2][!duplicated(link.select[ , 2])])
JS.fill <- melt.JS[melt.JS[ , 2] %in% cpd.select, ]
JS.fill[ , 1] <- paste('JS', JS.fill[ , 1], sep = ".")
JC.fill <- melt.JC[melt.JC[ , 2] %in% cpd.select, ]
JC.fill[ , 1] <- paste('JC', JC.fill[ , 1], sep = ".")
HS.fill <- melt.HS[melt.HS[ , 2] %in% cpd.select, ]
HS.fill[ , 1] <- paste('HS', HS.fill[ , 1], sep = ".")
HC.fill <- melt.HC[melt.HC[ , 2] %in% cpd.select, ]
HC.fill[ , 1] <- paste('HC', HC.fill[ , 1], sep = ".")

JI.fill <- rbind(JS.fill, JC.fill)
HZ.fill <- rbind(HS.fill, HC.fill)

link.fill <- rbind(JI.fill, HZ.fill)

link.wide <- dcast(link.fill, Var1 ~ Var2)
rownames(link.wide) <- link.wide[ , 1]
link.wide <- link.wide[ , -1]
link.wide[is.na(link.wide)] <- 0

min.link <- apply(link.wide, 1, min)
max.link <- apply(link.wide, 1, max)

link.wide <- link.wide[min.link <= -0.9 | max.link >= 0.9, ]

## pheatmap ##
annotation_row <- data.frame(class = rownames(link.wide))
rownames(annotation_row) <- rownames(link.wide)

### row color extract ##
my.substr <- function(str){
    str <- substr(str, 6, length(str))
    return(str)
}

module_colors <- as.character((apply(t(rownames(link.wide)), 1, my.substr)))
names(module_colors) <- rownames(link.wide)
ann_colors <- list(class = module_colors)
#pheatmap(link.wide, cluster_rows = T, color = c('#0c8918', rep('#f0f0f4', 7), '#ff4c00'), cellwidth = 10, cellheight = 10, legend = F, gaps_row = c(10, 21, 31), 
#         annotation_row = annotation_row, annotation_legend = F, annotation_colors = ann_colors, border_color = '#f0f0f4')

pheatmap(link.wide, cluster_rows = T, cellwidth = 10, cellheight = 10, legend = T, gaps_row = c(10, 21, 31),
         annotation_row = annotation_row, annotation_legend = F, annotation_colors = ann_colors, border_color = '#f0f0f4',
         filename = "./result/cor_metabolites/PM_cor.pdf", height = 10, width = 16.18)

