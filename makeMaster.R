args = commandArgs(trailingOnly=TRUE)

setwd(paste0('/fs/project/PAS0854/Active_projects/', args[1]))

library(GenomicDataCommons)

os = read.table('Orien.SNV.maf', sep='\t', header=T, fill = T)
ts = read.table('TCGA.SNV.maf', sep="\t", header=T, fill = T)

dam = (grepl('delet', ts$SIFT) & grepl('damag', ts$PolyPhen) & !(grepl('benign', ts$CLIN_SIG))) | 
  grepl("Frame_Shift", ts$Variant_Classification) 
ts = ts[dam, ]

dam = (grepl('D', os$SIFT_pred) & grepl('D', os$Polyphen2_HDIV_pred)) | 
  grepl("Frame_Shift", os$ExonicFunc.refGene) 
os = os[dam, ]

rm(dam)

fus = read.table('new.fusion.txt', sep="\t", header=T, fill = T)
newnames = gsub("\\..*", "", fus$Sample)
newnames = unlist(lapply(newnames, function(x) UUIDtoBarcode(x, from_type = "file_id")[, 2]))
fus$Sample = newnames

copy = read.table('copy.num.txt', sep="\t", fill = T)
copy = copy[copy[, 7] >= 6, ]
newnames = gsub("\\..*", "", copy[, 1])
newnames = unlist(lapply(newnames, function(x) UUIDtoBarcode(x, from_type = "file_id")[1, 2]))
copy[,1] = newnames

imp = read.csv('Imprecise.G.D.M.res.csv', row.names=NULL)
sv = read.csv('G.D.M.res.csv', row.names=NULL)

library('TCGAutils')

names = c(ts$Tumor_Sample_Barcode, os[, 1], imp$Sample, sv$Sample, copy[, 1], fus[, 1])
names = unique(names)
new = c()
for (x in names){
  e = strsplit(x, "-")
  if (length(e[[1]]) == 5){
    e = paste(e[[1]][3:length(e[[1]])], collapse="-")
  } else if (length(e[[1]]) == 2){
    e = e[[1]][1]
  } else {
    e = paste(e[[1]][1:3], collapse="-")
  }
  new = c(new, e)
}

new = unique(new)

new1 = c()
for (x in new){
  if (grepl('SL', x)){
    x = gsub("\\_.*", "", x)
  }
  new1 = c(new1, x)
}

new1 = unique(new1)
paste('Number of positive cases:', length(new1))

newTCGA = new1[grepl('TCGA', new1)]

newOrien = new1[!(grepl('TCGA', new1))]

rm(new)
rm(new1)
rm(names)
rm(e)
rm(x)

###################### do TCGA first #####################

library(TCGAbiolinks)

caseIDs = barcodeToUUID(newTCGA)

rnaIDs = read.csv('/fs/ess/PAS0854/Raven/Immune_Checkpoints/TCGA_Cibersort_Results_Categorized.csv')
rnaIDs = rnaIDs$Input.Sample

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

TCGAdat= cbind(newTCGA, rnaIDs, wxgs, ls)

rm(rnaIDs)
rm(wxgs)
rm(ls)
rm(newTCGA)

########################## now get mut data ########################

clas = c()
for (x in TCGAdat[, 1]){
  if (T %in% grepl(x, ts$Tumor_Sample_Barcode)){
    clas = c(clas, 'SNV')
  } else {
    clas = c(clas, 'no SNV')
  }
}

clas2 = c()
for (x in TCGAdat[, 1]){
  if (T %in% grepl(x, sv$Sample)){
    clas2 = c(clas2, 'Precise_SV')
  } else if (T %in% grepl(x, imp$Sample)){
    clas2 = c(clas2, 'Imprecise_SV')
  } else {
    clas2 = c(clas2, 'no SV')
  }
}

clas3 = c()
for (x in TCGAdat[, 1]){
  if (T %in% grepl(x, copy[, 1])){
    clas3 = c(clas3, 'CNV')
  } else {
    clas3 = c(clas3, 'no CNV')
  }
}

clas4 = c()
for (x in TCGAdat[, 1]){
  if (T %in% grepl(x, fus[, 1])){
    clas4 = c(clas4, 'Fusion')
  } else {
    clas4 = c(clas4, 'no Fusion')
  }
}

clas1 = cbind(clas, clas2, clas3, clas4)
type = c()
for (x in 1:nrow(clas1)){
  row = clas1[x, ]
  row[grepl('no', row)] = NA
  row = row[!(is.na(row))]
  type = c(type, paste0(row, collapse = ";"))
}
TCGAdat = cbind(TCGAdat, type)

rm(clas)
rm(clas2)
rm(clas3)
rm(clas4)
rm(clas1)
rm(type)

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

######################### now do Orien ################################

newOrien = newOrien[grepl("SL", newOrien)]
OrienClin = read.csv('/fs/ess/PAS0854/Raven/fusionsClin.csv')
OrienClin = OrienClin[OrienClin$WES %in% newOrien, ]
Oriendat = cbind(OrienClin$ORIENAvatarKey, OrienClin$RNASeq, OrienClin$WES, NA)

clas = c()
for (x in Oriendat[, 3]){
  if (T %in% grepl(x, os$Sample)){
    clas = c(clas, 'SNV')
  } else {
    clas = c(clas, 'no SNV')
  }
}

clas2 = c()
for (x in Oriendat[, 3]){
  if (T %in% grepl(x, sv$Sample)){
    clas2 = c(clas2, 'Precise_SV')
  } else if (T %in% grepl(x, imp$Sample)){
    clas2 = c(clas2, 'Imprecise_SV')
  } else {
    clas2 = c(clas2, 'no SV')
  }
}

clas3 = c()
for (x in Oriendat[, 3]){
  if (T %in% grepl(x, copy[, 1])){
    clas3 = c(clas3, 'CNV')
  } else {
    clas3 = c(clas3, 'no CNV')
  }
}

clas4 = c()
for (x in Oriendat[, 3]){
  if (T %in% grepl(x, fus[, 1])){
    clas4 = c(clas4, 'Fusion')
  } else {
    clas4 = c(clas4, 'no Fusion')
  }
}

clas1 = cbind(clas, clas2, clas3, clas4)
type = c()
for (x in 1:nrow(clas1)){
  row = clas1[x, ]
  row[grepl('no', row)] = NA
  row = row[!(is.na(row))]
  type = c(type, paste0(row, collapse = ";"))
}
Oriendat = cbind(Oriendat, type)

rm(clas)
rm(clas2)
rm(clas3)
rm(clas4)
rm(type)
rm(x)
rm(clas1)

Oriendat = cbind(Oriendat, OrienClin[, c("SpecimenSiteOfOrigin", 'Histology.Behavior')], NA)

rm(OrienClin)

######################### Combine ##########################

colnames(TCGAdat) = c('ID', 'RNAseq', 'WXS.WGS', 'SNV.file', 'Type', 
                      'TissueOfOrigin', "Histology", "Stage")
colnames(Oriendat) = c('ID', 'RNAseq', 'WXS.WGS', 'SNV.file', 'Type', 
                      'TissueOfOrigin', "Histology", "Stage")

master = rbind(TCGAdat, Oriendat)

write.csv(master, 'Master.csv')
