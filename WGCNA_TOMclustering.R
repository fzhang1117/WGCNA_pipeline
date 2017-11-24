## WGCNA_TOMclustering.R ##
## Clustering dissTOM matrix to basic module ##
## Zhang Fei ##
## 2017-11-24 ##

# We use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes. WGCNA using the same name function hclust that provides a much faster hierarchical clusering routine than the standard hclust function.

# Rscript WGCNA_TOMclustering.R <TOM matrix> <cor_method> <minModuleSize> <deepSplit(from 0 to 4)>
argv <- commandArgs(T)

library(WGCNA)

## load the dissTOM matrix
load(argv[1])


fl.diss <- ls(pattern = "dissTOM")

for(i in 1: length(fl.diss)){
	geneTree = hclust(as.dist(get(fl.diss[i])), method = "average")
	print(fl.diss[i])
	print(dim(get(fl.diss[i])))
	minModuleSize <- as.integer(argv[3])
	dynamicMods <- cutreeDynamic(dendro = geneTree, distM = get(fl.diss[i]), deepSplit = as.integer(argv[4]), minClusterSize = minModuleSize)
	dynamicColors <- labels2colors(dynamicMods)
	write.table(dynamicMods, paste("./result/mods.", fl.diss[i], ".", argv[2], ".", "deepsplit=", argv[4], ".txt", sep = ""), quote = F, sep = "\t", col.names = F)
	write.table(dynamicColors, paste("./result/colors.", fl.diss[i], ".", argv[2], ".", "deepsplit=", argv[4], ".txt", sep = ""), quote = F, sep = "\t", col.names = F)
	pdf(paste("./result/", fl.diss[i], ".", argv[2], ".", "deepsplit=", argv[4], ".pdf", sep = ""), height = 9, width = 12)
		#plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
		plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = F, hang = 0.03, addGuide = T, guideHang = 0.05, main = "Gene dendrogram and module colors")
	dev.off()
}

