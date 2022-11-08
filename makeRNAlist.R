library(TCGAutils)
library(GenomicDataCommons)
library(dplyr)

mas = read.csv('Master.csv')
con = read.csv('Controls.csv')

q = files() %>%
  GenomicDataCommons::select(available_fields('files')) %>%
  GenomicDataCommons::filter(~ data_type == 'Aligned Reads' &
           experimental_strategy == 'RNA-Seq' &
           data_format == 'BAM' &
           analysis.workflow_type == 'STAR 2-Pass Genome')
file_ids = q %>%  ids()

set = UUIDtoUUID(file_ids, to_type = 'case_id')

our = c(mas$ID[grepl('TCGA', mas$ID)], con$ID[grepl('TCGA', con$ID)])

our = barcodeToUUID(our)

set = set[set$cases.case_id %in% our$case_id, ]

writeLines(set$file_id, 'RNAseqBams.TCGA.txt')

#clin = read.csv('/fs/ess/PAS0854/Raven/fusionsClin.csv')
#clin = clin[clin$ORIENAvatarKey %in% mas$ID, ]
#clin = clin[!(clin$RNASeq == ""), ]

#our = c(paste0('OSU_IMP_Roychowdhury:/', clin$Project, '/RNAseq/alignment_crams/', clin$RNASeq, '.genome.cram'))

#writeLines(our, 'RNAseqBams.Orien.txt')
