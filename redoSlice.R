#read in the arguments from command line
args = commandArgs(trailingOnly = T)
name = as.character(args[1])

#read in the manifest

library(GenomicDataCommons)

q = files() %>%
  GenomicDataCommons::select(available_fields('files')) %>%
  filter(~ data_type == 'Aligned Reads' &
           experimental_strategy == 'WXS' | experimental_strategy == 'WGS' &
           data_format == 'BAM')
file_ids = q %>%  ids()

#make a list of all the files to download
fl = unique(file_ids)
#go to the bams location
setwd(paste0(
  '/fs/ess/PAS0854/Active_projects/',
  name,
  '/bams'
))

#make a list of all the bams located in that folder
flout = list.files(pattern = '.bam')
#check if there are in fact bams in that folder
if (length(flout) > 0) {
  #record the file IDs of the bams
  flout = gsub('\\..*', "", flout)
  #records the indexes of the bams
  ids = unlist(lapply(flout, function(x)
    return(which(fl == x))))
  #get new starting point from the max index
  if (length(ids) > 0 ){
    unlink(paste0(flout[which(flout == fl[max(ids)])], '.bam'))   
    newid = max(ids)
  } else {
    newid = 1
  }
  #remove the last bam in case it didnt finish downloading
  #make new list of bams to download
  fl = fl[newid:length(fl)]
}

#go to the home directory of the project
setwd(paste0('/fs/ess/PAS0854/Active_projects/', name))

#check and see if there is already toSlice file
test = list.files(pattern="toSlice.txt")
if (length(test) > 0){
  #if there is, delete it
  unlink('toSlice.txt') 
}

#make a new toSlice file with the files in fl
writeLines(fl, 'toSlice.txt')
