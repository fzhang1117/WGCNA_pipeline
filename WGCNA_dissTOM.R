## WGCNA network construct ##
## From expression Matrix to TOM ##
## Zhang Fei ##
## 2017-11-24 ##

## From expression matrix to cor matrix to TOM matrix

library(WGCNA)
load("./data/ave_expression.Rdata")

my.dissTOM <- function(expression, filebase = "myWGCNA", method.cor = "pearson", type.net = 'unsigned', verbose = 5){
    # This function soft the softthreshold and calculate dissTOM from expression data
    # cor.method could be either 'pearson', 'bicor' and 'spearman'
    # output list
        # cor matrix
        # TOM matrix
        # network
    
    # This is the softhreshold select codes
    expression <- t(expression)
    powers <- c(c(1: 10), seq(from = 12, to = 20, by = 2))
    sft <- pickSoftThreshold(expression, powerVector = powers, verbose = verbose)
    pdf(paste("./result/",  filebase, ".pdf", sep = ""))
        mat <- matrix(c(1, 2), 1, 2)
        layout(mat)
        plot(sft$fitIndices[ , 1], -sign(sft$fitIndices[ , 3])*sft$fitIndices[ , 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = "Scale independence")
        text(sft$fitIndices[ , 1], -sign(sft$fitIndices[ , 3])*sft$fitIndices[ , 2], labels = powers, cex = 1, col = 'red')
        abline(h = 0.90, col = 'red')
        plot(sft$fitIndices[ , 1], sft$fitIndices[ , 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
        text(sft$fitIndices[ , 1], sft$fitIndices[ , 5], labels = powers, cex = 1, col = 'red')
    dev.off()
    
    # Co-expression simlarity and adjacency
    softPower <- sft$powerEstimate
    if(method.cor == 'pearson'){
        adjacency <- adjacency(datExpr = expression, type = type.net, power = softPower, corFnc = "cor", corOptions = "use = 'p'")
    }else if(method.cor == 'spearman'){
        adjacency <- adjacency(datExpr = expression, type = type.net, power = softPower, conFnc = "cor", corOptions = "method = 'spearman")
    }else if(method.cor == 'bicor'){
        adjacency <- adjacency(datExpr = expression, type = type.net, power = softPower, conFnc = "bicor")
    }
    
    
    #Topological Overlap Matrix (TOM)
    
    TOM = TOMsimilarity(adjacency, TOMType = type.net, verbose = verbose)
    dissTOM = 1 - TOM
    
    return(dissTOM)    
}

