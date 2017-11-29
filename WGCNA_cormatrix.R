## WGCNA ##
## Zhang Fei ##
## 2017-11-22 ##

library(WGCNA)

# step 1: import data
# different data format, for hclust
my.data <- t(read.table("./demo_data/Demo.txt", header = T, row.names = 1, sep = "\t"))

# setp 2: soft thresholding power 'beta'
powers = c(c(1: 10), seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(my.data, powerVector = powers, verbose = 5) ## The best sft is sft$powerEstimate 

pdf("./result/sft.pdf")
    mat <- matrix(c(1, 2), 1, 2)
    layout(mat)
    plot(sft$fitIndices[ , 1], -sign(sft$fitIndices[ , 3])*sft$fitIndices[ , 2], 
        xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", 
        main = "Scale independence")
    text(sft$fitIndices[ , 1], -sign(sft$fitIndices[ , 3])*sft$fitIndices[ , 2], labels = powers, cex = 1, col = 'red')
    abline(h = 0.90, col = 'red')
    plot(sft$fitIndices[ , 1], sft$fitIndices[ , 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
    text(sft$fitIndices[ , 1], sft$fitIndices[ , 5], labels = powers, cex = 1, col = 'red')
    dev.off()

# step 3: build co-expression Network
net = blockwiseModules(
    my.data, corType = "bicor", deepSplit = 4,
    power = sft$powerEstimate,
    maxBlockSize = 1000,
    TOMType = "unsigned", minModuleSize = 30,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = T, pamRespectsDendro = F,
    saveTOMs = T, saveTOMFileBase = 'Test-TPM-TOM',
    verbose = 3
)



# step 4: visualize the moudle
pdf("./result/bicor_module_deep4.pdf")
    mergedColors = labels2colors(net$colors)
    plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = F, hang = 0.03, addGuide = T, guideHang = 0.05)
dev.off()

# step 5: merge moudles


# step 6: Enrichment