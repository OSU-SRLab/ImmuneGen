#load the necesary libraries
library(GenomicDataCommons)
library(doParallel)

#load in the arguments from command line
args = commandArgs(trailingOnly=T)
name = as.character(args[1])
reg = as.character(args[2])

#set the token
Sys.setenv(GDC_TOKEN = readLines('/fs/ess/PAS0854/Reference_Data/GDC_Token/gdc-user-token.2022.10.19.txt'))

#get the list of files to download
fl = readLines(paste0('/fs/ess/PAS0854/Active_projects/', name, '/toSlice.txt'))

#download and slice

foreach(x=fl) %dopar% { 
         GenomicDataCommons::slicing(x,
                 regions = reg,
                 destination = file.path(paste0('/fs/ess/PAS0854/Active_projects/', name, '/bams/', x, '.bam'))
         )
         print(x)
       }
