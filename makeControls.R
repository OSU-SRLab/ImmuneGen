args = commandArgs(trailingOnly=TRUE)

setwd(paste0('/fs/project/PAS0854/Active_projects/', args[1]))

master = read.csv('Master.csv')

newTCGA = master$ID

rnaIDs = read.csv('/fs/ess/PAS0854/Raven/Immune_Checkpoints/TCGA_Cibersort_Results_Categorized.csv')
rnaIDs = rnaIDs$Input.Sample

offTCGA = strsplit(rnaIDs, "-")
offTCGA = unlist(lapply(offTCGA, function(x) return(paste0(x[1:3], collapse = "-"))))

offTCGA = offTCGA[!(offTCGA %in% newTCGA)]

offTCGA = offTCGA[sample(length(offTCGA), size=500)]

newTCGA = offTCGA
rm(offTCGA)

###################### do TCGA first #####################

library(TCGAbiolinks)
library(TCGAutils)
library(GenomicDataCommons)

caseIDs = barcodeToUUID(newTCGA)

rnaIDs1 = c()
for (x in newTCGA){
  if (T %in% grepl(x, rnaIDs)){
    id = paste0(rnaIDs[which(grepl(x, rnaIDs) == T)], collapse= ';')
  } else {
    id = 'none'
  }
  rnaIDs1 = c(rnaIDs1, id)
}

rnaIDs = rnaIDs1
rm(rnaIDs1)

wxgs = readLines(paste0('/fs/project/PAS0854/Active_projects/', args[1], '/CD40.PairedBams.txt'))
wxgs = unlist(lapply(wxgs, function(x) return(strsplit(x, ";")[[1]][2])))
wxgs = unique(wxgs)

wxgs1 = c()
for (x in newTCGA){
  if (T %in% grepl(x, wxgs)){
    id = paste0(wxgs[which(grepl(x, wxgs) == T)], collapse=";")
  } else {
    id = 'none'
  }
  wxgs1 = c(wxgs1, id)
}

wxgs = wxgs1
rm(wxgs1)

ls = list.files('/fs/ess/PAS0854/Active_projects/TCGA_novel_icktps/TCGA.SomaticMut.Mafs/')
ls = filenameToBarcode(ls)
ls1 = c()
for (x in newTCGA){
  if (T %in% grepl(x, ls$aliquots.submitter_id)){
    id = paste0(ls$file_name[which(grepl(x, ls$aliquots.submitter_id) == T)], 
                collapse=";")
  } else {
    id = 'none'
  }
  ls1 = c(ls1, id)
}

ls = ls1
rm(ls1)

TCGAdat= cbind(newTCGA, rnaIDs, wxgs, ls, 'control')

rm(rnaIDs)
rm(wxgs)
rm(ls)
rm(newTCGA)

########################## now get clinical data ########################

clin = gdc_clinical(caseIDs$case_id)
clin1 = clin$diagnoses
clin1 = clin1[match(caseIDs$case_id, clin1$case_id), ]
clin1 = clin1[, c("primary_diagnosis", 
                  "ajcc_pathologic_stage")]
clin2 = clin$main
clin2 = clin2[match(caseIDs$case_id, clin2$case_id), ]
clin2 = clin2[, c('primary_site')]
TCGAdat = cbind(TCGAdat, clin2, clin1)

rm(clin)
rm(clin2)
rm(clin1)

######################### Combine ##########################

colnames(TCGAdat) = c('ID', 'RNAseq', 'WXS.WGS', 'SNV.file', 'Type', 
                      'TissueOfOrigin', "Histology", "Stage")

write.csv(TCGAdat, 'Controls.csv')

