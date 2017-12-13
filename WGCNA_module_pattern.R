## WGCNA_module_pattern.R ##
## Zhang Fei ##
## 2017-12-04 ##
## Plot the module expression pattern ##

load("./data/ave_expression.Rdata")

#eigengene <- read.table("./result/eigengene/MEs.dissTOM.HS.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, row.names = 1, sep = "\t")

## data prepare ##

### ave module ###
AVE.JS <- read.table("./result/modules/AEs.dissTOM.JS.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, row.names = 1, sep = "\t")
AVE.JC <- read.table("./result/modules/AEs.dissTOM.JC.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, row.names = 1, sep = "\t")
AVE.HS <- read.table("./result/modules/AEs.dissTOM.HS.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, row.names = 1, sep = "\t")
AVE.HC <- read.table("./result/modules/AEs.dissTOM.HC.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, row.names = 1, sep = "\t")

### module import ###
module.JS <- read.table("./result/modules/merged_modes.dissTOM.JS.pearson.deepsplit=2.MEDissThres=0.10.txt", header = F, row.names = 1, sep = "\t")
modulename.JS <- levels(module.JS[ , 2])
module.JC <- read.table("./result/modules/merged_modes.dissTOM.JC.pearson.deepsplit=2.MEDissThres=0.10.txt", header = F, row.names = 1, sep = "\t")
modulename.JC <- levels(module.JC[ , 2])
module.HS <- read.table("./result/modules/merged_modes.dissTOM.HS.pearson.deepsplit=2.MEDissThres=0.10.txt", header = F, row.names = 1, sep = "\t")
modulename.HS <- levels(module.HS[ , 2])
module.HC <- read.table("./result/modules/merged_modes.dissTOM.HC.pearson.deepsplit=2.MEDissThres=0.10.txt", header = F, row.names = 1, sep = "\t")
modulename.HC <- levels(module.HC[ , 2])

### expression scale ###
expression.JS <- t(apply(ave.JS, 1, scale))
expression.JC <- t(apply(ave.JC, 1, scale))
expression.HS <- t(apply(ave.HS, 1, scale))
expression.HC <- t(apply(ave.HC, 1, scale))


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

my.patterndrawing <- function(file, list.module, expression, modulename, AVE){
    pdf(file, width = 7, height = 7)
    for(i in 1: length(modulename)){
        plot(c(1, 15), c(min(expression), max(expression)), type = "n", xlab = "time point", ylab = "expression", main = modulename[i])
        module.AE <- paste("AE", modulename[i], sep = "")
        module.ME <- paste("ME", modulename[i], sep = "")
        for(j in 1: dim(list.module[[i]])[1]){
            lines(c(1: 15), list.module[[i]][j, ], col = "grey")
        }
        #lines(c(1: 15), eigengene[ , i], col = 'blue', type = "o", pch = 16)
        lines(c(1: 15), AVE[ , i], col = 'midnightblue', type = "o", pch = 16, lwd = 2)
    }
    
    dev.off()
    return(0)
}

lst.JS <- my.moduleExtract(module = module.JS, modulename = modulename.JS, expression = expression.JS)
lst.JC <- my.moduleExtract(module = module.JC, modulename = modulename.JC, expression = expression.JC)
lst.HS <- my.moduleExtract(module = module.HS, modulename = modulename.HS, expression = expression.HS)
lst.HC <- my.moduleExtract(module = module.HC, modulename = modulename.HC, expression = expression.HC)

my.patterndrawing("./result/module_pattern/pattern.JS.pdf", list.module = lst.JS, expression = expression.JS, modulename = modulename.JS, AVE = AVE.JS)
my.patterndrawing("./result/module_pattern/pattern.HS.pdf", list.module = lst.HS, expression = expression.HS, modulename = modulename.HS, AVE = AVE.HS)
my.patterndrawing("./result/module_pattern/pattern.JC.pdf", list.module = lst.JC, expression = expression.JC, modulename = modulename.JC, AVE = AVE.JC)
my.patterndrawing("./result/module_pattern/pattern.HC.pdf", list.module = lst.HC, expression = expression.HC, modulename = modulename.HC, AVE = AVE.HC)