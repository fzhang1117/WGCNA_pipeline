## WGCNA_module_pattern.R ##
## Zhang Fei ##
## 2017-12-04 ##
## Plot the module expression pattern ##

load("./data/ave_expression.Rdata")

eigengene <- read.table("./result/eigengene/MEs.dissTOM.HS.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, row.names = 1, sep = "\t")
AVE <- read.table("./result/modules/AEs.dissTOM.HS.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, row.names = 1, sep = "\t")
module <- read.table("./result/modules/merged_modes.dissTOM.HS.pearson.deepsplit=2.MEDissThres=0.10.txt", header = F, row.names = 1, sep = "\t")
modulename <- levels(module[ , 2])
#modulename <- modulename[modulename != 'grey']
expression <- t(apply(ave.HS, 1, scale))

my.moduleExtract <- function(module, modulename, expression){
    output = list()
        for(i in 1: length(modulename)){
            module.gene <- module[module[ , 2] == modulename[i], ]
            module.expression <- expression[rownames(expression) %in% module.gene[ , 1], ]
            colnames(module.expression) <- c(1: 15)
            output[[i]] <- module.expression
    }
    names(output) <- modulename
    return(output)
}

list.module <- my.moduleExtract(module = module, modulename = modulename, expression = expression)

pdf("test2.pdf", width = 7, height = 7)
    for(i in 1: length(modulename)){
            plot(c(1, 15), c(min(expression), max(expression)), type = "n", xlab = "time point", ylab = "expression", main = modulename[i])
            module.AE <- paste("AE", modulename[i], sep = "")
            module.ME <- paste("ME", modulename[i], sep = "")
            for(j in 1: dim(list.module[[i]])[1]){
                lines(c(1: 15), list.module[[i]][j, ], col = "grey")
            }
            #lines(c(1: 15), eigengene[ , i], col = 'blue', type = "o", pch = 16)
            lines(c(1: 15), AVE[ , i], col = 'red', type = "o", pch = 18)
        }
        
        

dev.off()
# points(c(1: 15), eigengene$module.ME, pch = 16, col = 'red')
# lines(c(1: 15), eigengene$module.ME, col = 'red')
# 
# 
# lines(c(1: 15), list.module$black[1, ], col = "grey")
# lines(c(1: 15), eigengene$MEblack, col = 'red')
