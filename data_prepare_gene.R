## WGCNA data_prepare ##
## import and pre-processing data for next step analysis ##
## Zhang Fei ##
## 2017-11-23 ##

mydata <- read.table("H://inMP/time-series/expression_data/tpm.combat.txt", header = T, sep = " ", row.names = 1, check.names = F)
index.JI <- read.table("H://inMP/time-series/expression_data/index_JI.txt", header = F, sep = "\t", row.names = 1, check.names = F)
index.HZ <- read.table("H://inMP/time-series/expression_data/index_HZ.txt", header = F, sep = "\t", row.names = 1, check.names = F)

expression.JI <- mydata[rownames(mydata) %in% index.JI[ , 1], c(60: 1, 120: 61)][ , c(31: 60, 91: 120)]
expression.HZ <- mydata[rownames(mydata) %in% index.HZ[ , 1], c(60: 1, 120: 61)][ , c(1: 30, 61: 90)]

expression.JS <- expression.JI[ , c(16: 30, 46: 60)]
expression.JC <- expression.JI[ , c(1: 15, 31: 45)]
expression.HS <- expression.HZ[ , c(16: 30, 46: 60)]
expression.HC <- expression.HZ[ , c(1: 15, 31: 45)]

ave.JS <- (expression.JS[ , c(1: 15)] + expression.JS[ , c(16: 30)])/2
ave.JC <- (expression.JC[ , c(1: 15)] + expression.JC[ , c(16: 30)])/2
ave.HS <- (expression.HS[ , c(1: 15)] + expression.HS[ , c(16: 30)])/2
ave.HC <- (expression.HC[ , c(1: 15)] + expression.HC[ , c(16: 30)])/2

timepoint <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15')
colnames(ave.JS) <- paste('JS', timepoint, sep = "")
colnames(ave.JC) <- paste('JC', timepoint, sep = "")
colnames(ave.HS) <- paste('HS', timepoint, sep = "")
colnames(ave.HC) <- paste('HC', timepoint, sep = "")

save('ave.JS', 'ave.JC', 'ave.HS', 'ave.HC' ,list = c('ave.JS', 'ave.JC', 'ave.HS', 'ave.HC'), file = "./data/ave_expression.Rdata")
