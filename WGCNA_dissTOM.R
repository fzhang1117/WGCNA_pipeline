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
    powers <- c(c(1: 10), seq(from = 12, to = 100, by = 2))
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
        adjacency <- adjacency(datExpr = expression, type = type.net, power = softPower, corFnc = "cor", corOptions = "use = 'p', method = 'spearman'")
    }else if(method.cor == 'bicor'){
        adjacency <- adjacency(datExpr = expression, type = type.net, power = softPower, corFnc = "bicor")
    }
    
    
    #Topological Overlap Matrix (TOM)
    
    TOM = TOMsimilarity(adjacency, TOMType = type.net, verbose = verbose)
    dissTOM = 1 - TOM
    
    return(dissTOM)    
}

#dissTOM.JC <- my.dissTOM(ave.JC, filebase = "JC.pcc", method.cor = "pearson", type.net = 'unsigned', verbose = 5)
#dissTOM.JS <- my.dissTOM(ave.JS, filebase = "JS.pcc", method.cor = "pearson", type.net = 'unsigned', verbose = 5)
#dissTOM.HC <- my.dissTOM(ave.HC, filebase = "HC.pcc", method.cor = "pearson", type.net = 'unsigned', verbose = 5)
#dissTOM.HS <- my.dissTOM(ave.HS, filebase = "HS.pcc", method.cor = "pearson", type.net = 'unsigned', verbose = 5)

#save('dissTOM.JC', 'dissTOM.JS', 'dissTOM.HC', 'dissTOM.HS', list = c('dissTOM.JC', 'dissTOM.JS', 'dissTOM.HC', 'dissTOM.HS'), file = "./data/dissTOM.pearson.Rdata")

dissTOM2.JC <- my.dissTOM(ave.JC, filebase = "JC.sp", method.cor = "spearman", type.net = 'unsigned', verbose = 5)
dissTOM2.JS <- my.dissTOM(ave.JS, filebase = "JS.sp", method.cor = "spearman", type.net = 'unsigned', verbose = 5)
dissTOM2.HC <- my.dissTOM(ave.HC, filebase = "HC.sp", method.cor = "spearman", type.net = 'unsigned', verbose = 5)
dissTOM2.HS <- my.dissTOM(ave.HS, filebase = "HS.sp", method.cor = "spearman", type.net = 'unsigned', verbose = 5)

save('dissTOM2.JC', 'dissTOM2.JS', 'dissTOM2.HC', 'dissTOM2.HS', list = c('dissTOM.JC', 'dissTOM.JS', 'dissTOM.HC', 'dissTOM.HS'), file = "./data/dissTOM.spearman.Rdata")

dissTOM3.JC <- my.dissTOM(ave.JC, filebase = "JC.bicor", method.cor = "bicor", type.net = 'unsigned', verbose = 5)
dissTOM3.JS <- my.dissTOM(ave.JS, filebase = "JS.bicor", method.cor = "bicor", type.net = 'unsigned', verbose = 5)
dissTOM3.HC <- my.dissTOM(ave.HC, filebase = "HC.bicor", method.cor = "bicor", type.net = 'unsigned', verbose = 5)
dissTOM3.HS <- my.dissTOM(ave.HS, filebase = "HS.bicor", method.cor = "bicor", type.net = 'unsigned', verbose = 5)

save('dissTOM3.JC', 'dissTOM3.JS', 'dissTOM3.HC', 'dissTOM3.HS', list = c('dissTOM.JC', 'dissTOM.JS', 'dissTOM.HC', 'dissTOM.HS'), file = "./data/dissTOM.bicor.Rdata")