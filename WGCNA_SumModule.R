## WGCNA_SumModule.R ##
## Summary basic informations of WGCNA modules ##
## Zhang Fei ##
## 2017-11-26 ##

net.list <- dir("./result/modules/", pattern = "merged_modes.dissTOM")

for(i in 1: length(net.list)){
    module <- read.table(paste("./result/modules/", net.list[i], sep = ""), header = F, sep = "\t")
    table <- as.data.frame(sort(table(module[ , 3]), decreasing = T))
    colnames(table) <- c("mergeColor", "Number")
    write.table(table, paste("./result/summary/module_size/", "summary_", net.list[i], sep = ""), quote = F, sep = "\t")
    gene.list <- as.data.frame(table(module[ , c(2, 3)]))
    gene.list <- gene.list[which(gene.list[ , 3] == 1), ]
    module.name <- levels(gene.list[ , 2])
    bingo.char <- c()
    #bingo.char <- append(bingo.char, 'batch')
    for(j in 1: length(module.name)){
        my.gene <- as.character(gene.list[which(gene.list[ , 2] == module.name[j]), 1])
        bingo.char <- append(bingo.char, module.name[j])
        bingo.char <- append(bingo.char, my.gene)
        bingo.char <- append(bingo.char, 'batch')
    }
    file.pools <- strsplit(net.list[i], split = "dissTOM", fixed = T)
    write.table(bingo.char, paste("./result/summary/module_genes/genes.bingo", file.pools[[1]][2], sep = "") , quote = F, row.names = F, col.names = F)
}

## GO_enrichment ##


## KEGG_enrichment ##