args = commandArgs(trailingOnly = T)
name = args[1]

library(reshape2)
library(ggplot2)

#go to the gene directory
setwd(paste0('/fs/ess/PAS0854/Active_projects/', name))
cib <- read.csv(list.files()[grep('CIBERSORT', list.files())])
key <- read.csv("LOF.GOF.key.csv")
mut <- read.csv('mut.stat.key.csv')

key2 = cbind(key$V1, key$func2)
colnames(key2) = c('V1', 'V2')

coldata = rbind(key2, cbind(mut$V1[mut$V2 == 'control'], 
                            mut$V2[mut$V2 == 'control']))
coldata = as.data.frame(coldata)

coldata = coldata[coldata[, 1] %in% cib$Mixture, ]
coldata = coldata[match(cib$Mixture, coldata$V1), ]

cib = cbind(cib$T.cells.CD8, cib$T.cells.regulatory..Tregs., cib$Macrophages.M1, cib$Macrophages.M2, coldata$V2)
colnames(cib) = c("CD8+_Tcells", "T_regs", "M1_Macrophage", "M2_Macrophage", "Function")
cib = as.data.frame(cib)
cibm = melt(cib, id.vars = c('Function'))
cibm$value = as.numeric(as.character(cibm$value))

my_comparisons <- list(  c("control", "unknown"), c("control", "LOF"),c("control", "GOF")  )

ggplot(cibm, aes(x = Function, y= value)) + 
  geom_violin(aes(fill = Function)) + 
  facet_wrap(~ variable) +
  stat_compare_means(comparisons = my_comparisons, label.y = c(.4, .45, .5)) +
  stat_compare_means(label.y = .6)
