## WGCNA_metabolicdatareshape.R ##
## Zhang Fei ##
## 2017-11-29 ##

library(reshape2)

## GC part ##
# 
# GC <- read.table("H://inMP/time-series/metabolic_data/GC_data_after_normalization.txt", header = T, check.names = F, row.names = 1, quote = "")
# location <- c(30,150,125,67,115,71,126,136,19,144,70,82,120,52,101,159,66,2,137,99,73,55,48,74,111,122,100,93,102,119,75,34,164,77,139,140,3,167,148,107,10,89,128,97,153,149,129,21,79,35,28,123,72,18,1,15,133,108,134,158,154,174,112,114,65,142,45,6,83,24,170,117,118,49,69,5,124,80,36,56,84,143,43,14,95,44,50,16,31,37,168,130,12,26,135,98,155,113,121,76,109,20,62,61,173,163,32,41,47,94,165,4,17,29,39,58,96,13,169,86,157,9,64,40,54,156,152,25,53,38,110,132,60,131,160,7,104,145,78,85,166,59,91,161,51,106,105,23,90,42,141,8,63,146,147,88,27,87,171)
# GC.new <- GC[ , location]
# 
# data.HC <- data.new[ ,c(1: 40)]
# data.HS <- data.new[ ,c(41: 80)]
# data.JC <- data.new[ ,c(81: 120)]
# data.JS <- data.new[ ,c(121: 159)]
# 
# timepoint.1 <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15)
# timepoint.2 <- c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15)
# 


### PM_part ###
GC <- read.table("./data/GC.new.txt", header = T, row.names = 1, quote = "", sep = "\t")
GC.melt <- melt(GC, id.vars = "Type", na.rm = T)
GC.ave <- dcast(GC.melt, Type ~ variable, mean)

### SM_part ###

SM_neg <- read.table("./data/SM.neg.txt", header = T, row.names = 1, quote = "", sep = "\t")
SM_neg.melt <- melt(SM_neg, id.vars = "Type", na.rm = T)
SM_neg.ave <- dcast(SM_neg.melt, Type ~ variable, mean)

SM_pos <- read.table("./data/SM.pos.txt", header = T, row.names = 1, quote = "", sep = "\t")
SM_pos.melt <- melt(SM_pos, id.vars = "Type", na.rm = T)
SM_pos.ave <- dcast(SM_pos.melt, Type ~ variable, mean)

### LP_Part ###

LP_neg.ave <- read.table("./data/LP.neg.txt", header = T, quote = "", sep = "\t")
LP_pos.ave <- read.table("./data/LP.pos.txt", header = T, quote = "", sep = "\t")

save(GC.ave, SM_neg.ave, SM_pos.ave, LP_neg.ave, LP_pos.ave, file = "./data/metabolic.Rdata")
