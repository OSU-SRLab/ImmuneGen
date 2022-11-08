# take in the command line arguments
args = commandArgs(trailingOnly = T)
name = as.character(args[1])

setwd(paste0(
'/fs/project/PAS0854/Active_projects/',
name,
'/bams'
))

library(TCGAutils)

#make a list of the bams in your folder
ls = list.files(pattern = 'bam')
ls = ls[!(grepl('bai', ls))]
#remove the bam extension
ls1 = gsub("\\..*", "", ls)
#conver to barcode IDs
ls2 = UUIDtoBarcode(unique(ls1), from_type = 'file_id')
#make a list of the barcodes from the manifest key
ls3 = unlist(lapply(ls1, function(x){
  id = which(ls2[, 1] == x);
  return(ls2[id, 2])
}))
#put the bam extension back on
ls3 = paste0(ls3, ".bam")
#rename
lapply(1:length(ls), function(x)
  file.rename(ls[x], ls3[x]))
