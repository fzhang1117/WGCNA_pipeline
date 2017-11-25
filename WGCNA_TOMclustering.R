## WGCNA_TOMclustering.R ##
## Clustering dissTOM matrix to basic module ##
## Zhang Fei ##
## 2017-11-24 ##

# We use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes. WGCNA using the same name function hclust that provides a much faster hierarchical clusering routine than the standard hclust function.

# Rscript WGCNA_TOMclustering.R <TOM matrix> <cor_method> <minModuleSize> <deepSplit(from 0 to 4)> <MEDissThres>
argv <- commandArgs(T)

library(WGCNA)

## load the dissTOM matrix
load(argv[1])
load("./data/ave_expression.Rdata")

fl.diss <- ls(pattern = "dissTOM")
fl.expression <- ls(pattern = "ave")
#print(fl.diss)
#print(fl.expression)
index.HZ <- list.index$HZ
index.JI <- list.index$JI


for(i in 1: length(fl.diss)){
	geneTree <- hclust(as.dist(get(fl.diss[i])), method = "average")
	print(fl.diss[i])
	print(fl.expression[i])
	minModuleSize <- as.integer(argv[3])
	dynamicMods <- cutreeDynamic(dendro = geneTree, distM = get(fl.diss[i]), deepSplit = as.integer(argv[4]), minClusterSize = minModuleSize)
	dynamicColors <- labels2colors(dynamicMods)
	
	#write.table(dynamicMods, paste("./result/mods.", fl.diss[i], ".", argv[2], ".", "deepsplit=", argv[4],".txt", sep = ""), quote = F, sep = "\t", col.names = F)
	#write.table(dynamicColors, paste("./result/colors.", fl.diss[i], ".", argv[2], ".", "deepsplit=", argv[4], ".txt", sep = ""), quote = F, sep = "\t", col.names = F)
	
	# merge similar modules
	# The dynamic tree cut may identify modules whose expression profiles are very similar. It may be prudent to merge such modules since their genes are highly co-expressed. To quantify co-expression similarity of entire modules, we calculate their eigengenes and cluster them on their correlation:
	# Calculate eigengenes
	MEList <- moduleEigengenes(t(get(fl.expression[i])), colors = dynamicColors, nPC = 1, excludeGrey = F)
	MEs <- MEList$eigengenes
	write.table(MEs, paste("./result/MEs.", fl.diss[i], ".", argv[2], ".", "deepsplit=", argv[4], ".txt", sep = ""), quote = F, sep = "\t", col.names = T)	
	# Calculate dissimilarity of module eigengenes
	MEDiss <- 1 - cor(MEs)
	
	# Cluster module eigengenes
	METree = hclust(as.dist(MEDiss), method = "average")
	print(dim(get(fl.diss[i])))
	MEDissThres = as.numeric(argv[5])

	# Call an automatic merging function
	merge = mergeCloseModules(t(get(fl.expression[i])), dynamicColors, cutHeight = argv[5], verbose = 3)
	mergedColors = merge$colors

	# new merged modules
	mergeMEs = merge$newMEs
	if(length(mergedColors) == 3707){
		mergedColors.out <- data.frame(Gene.id = index.JI, module = mergedColors)
	}else if(length(mergedColors) == 5768){
		mergedColors.out <- data.frame(Gene.id = index.HZ, module = mergedColors)
	}
	#mergeColors.output = data.frame(GeneName = list.index)
	write.table(mergedColors.out, paste("./result/merged_modes.", fl.diss[i], ".", argv[2], ".", "deepsplit=", argv[4], ".", "MEDissThres=", argv[5], ".txt", sep = ""), quote = F, sep = "\t", col.names = F)
	pdf(paste("./result/", fl.diss[i], ".", argv[2], ".", "deepsplit=", argv[4], ".MEDissThres=", argv[5],".pdf", sep = ""), height = 9, width = 12)
		#plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
		plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = F, hang = 0.03, addGuide = T, guideHang = 0.05, main = "Gene dendrogram and module colors")
		plot(METree, main = "Clustering of module eigengenes (PC1))", xlab = "", sub = "")
		abline(h = MEDissThres, col = 'red')
		plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLables = F, hang = 0.03, addGuide = T, guideHang = 0.05)
	dev.off()
}

