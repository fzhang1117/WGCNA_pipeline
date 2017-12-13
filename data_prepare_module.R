## data_prepare_module.R ##
## Zhang Fei ##
## 2017-12-13 ##
module.HS <- read.table("./result/modules/merged_modes.dissTOM.HS.pearson.deepsplit=2.MEDissThres=0.10.txt", header = F, row.names = 1, sep = "\t")
colnames(module.HS) <- c("gene", "module")
ME.HS <- read.table("./result/modules/MEs.dissTOM.HS.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t")
AE.HS <- read.table("./result/modules/AEs.dissTOM.HS.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t")

module.HC <- read.table("./result/modules/merged_modes.dissTOM.HC.pearson.deepsplit=2.MEDissThres=0.10.txt", header = F, row.names = 1, sep = "\t")
colnames(module.HC) <- c('gene', "module")
ME.HC <- read.table("./result/modules/MEs.dissTOM.HC.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t")
AE.HC <- read.table("./result/modules/AEs.dissTOM.HC.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t")

module.JS <- read.table("./result/modules/merged_modes.dissTOM.JS.pearson.deepsplit=2.MEDissThres=0.10.txt", header = F, row.names = 1, sep = "\t")
colnames(module.JS) <- c("gene", "module")
ME.JS <- read.table("./result/modules/MEs.dissTOM.JS.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t")
AE.JS <- read.table("./result/modules/AEs.dissTOM.JS.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t")

module.JC <- read.table("./result/modules/merged_modes.dissTOM.JC.pearson.deepsplit=2.MEDissThres=0.10.txt", header = F, row.names = 1, sep = "\t")
colnames(module.JC) <- c("gene", "module")
ME.JC <- read.table("./result/modules/MEs.dissTOM.JC.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t")
AE.JC <- read.table("./result/modules/AEs.dissTOM.JC.pearson.deepsplit=2MEDissThres=0.10.txt", header = T, sep = "\t")

save(module.HS, module.HC, module.JS, module.JC, ME.HS, ME.HC, ME.JS, ME.JC, AE.HS, AE.HC, AE.JS, AE.JC, file = "./data/modules.signed.pearson.Rdata")