## ANOVA_metabolites.R ##
## Zhang Fei ##
## 2017-12-13 ##

load("./data/metabolic.Rdata")

my.ANOVA <- function(metabolic.matrix){
    meta.new <- t(metabolic[ , -1])
    colnames(meta.new) <- metabolic[ , 1]
    meta.HZ <- meta.new[ , c(1: 30)]
    meta.JI <- meta.
    
}