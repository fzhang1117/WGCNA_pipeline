## WGCNA_HubSearching.R ##
## Searching hub genes of each module ##
## Zhang Fei ##
## 2017-12-13 ##

library(WGCNA)
load("./data/ave_expression.Rdata")
load("./data/modules.signed.pearson.Rdata")
load("./data/metabolic.Rdata")


#meta.test <- GC.ave$`N,N'-Diacetylchitobiose (1MEOX) (6TMS) MP`[46: 60]

my.hubsearch <- function(datExpr, module, sftpower = 1, phenotype, file = "genesearch.pdf", type = "HS"){
    if(type == 'HC'){
        phenotype <- phenotype[1: 15]
    }else if(type == 'HS'){
        phenotype <- phenotype[16: 30]
    }else if(type == 'JC'){
        phenotype <- phenotype[31: 45]
    }else if(type == 'JS'){
        phenotype <- phenotype[46: 60]
    }
    datExpr <- t(datExpr)
    adjacency <- abs(cor(datExpr))^sftpower
    colorh1 <- module[ , 2]
    alldegree <- intramodularConnectivity(adjacency, colorh1)
    GeneSignificance <- (cor(datExpr, phenotype, use = "p"))
    colorlevels <- unique(colorh1)
    pdf(file = file, width = 7, height = 7)
    #sizeGrWindow(9, 6)
    #par(mfrow = c(2, as.integer(0.5 + length(colorlevels)/2)))
    #par(mar = c(4, 5, 3, 1))
    for(i in c(1: length(colorlevels))){
        whichmodule <- colorlevels[[i]]
        restrict <- (colorh1 == whichmodule)
        verboseScatterplot(alldegree$kWithin[restrict],
        GeneSignificance[restrict], col = colorh1[restrict],
        main = whichmodule, 
        xlab = 'Connectivity', ylab = "Gene Significance", abline = F, pch = 16
        )
    }
    dev.off()
    
    return(alldegree)
}

my.resultout <- function(metamatrix, datExpr, module, type){
    for(i in 1: (dim(metamatrix)[2] - 1)){
        Degree <- my.hubsearch(datExpr = datExpr, module = module, sftpower = 1, phenotype = metamatrix[ , (i + 1)], file = paste("./result/GeneSearch/plot/", type, ".", colnames(metamatrix[i + 1]), ".pdf", sep = ""), type = type)
        write.table(Degree, file = paste("./result/GeneSearch/degree/", type, ".", colnames(metamatrix[i + 1]), ".txt", sep = ""), quote = F)   
    }
    gc()
}

my.resultout(GC.ave, ave.HS, module.HS, type = "HS")
my.resultout(GC.ave, ave.HC, module.HC, type = "HC")
my.resultout(GC.ave, ave.JS, module.JS, type = "JS")
my.resultout(GC.ave, ave.JC, module.JC, type = "JC")