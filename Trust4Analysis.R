args = commandArgs(trailingOnly = T)

name = args[1]

setwd(paste0('/fs/ess/PAS0854/Active_projects/', name, "/trust4"))

library(seqinr)
library(ape)
library(msa)
library(data.table)
library(reshape2)
library(plyr)
library(dplyr)
library(stringr)
library(tidyr)

##function of BCR clustering 
SeqDist <- function(x, y) {
  if (nchar(x) != nchar(y))
    return(NA)
  x.l = unlist(strsplit(x, ''))
  y.l = unlist(strsplit(y, ''))
  nn1 = length(which(x.l != y.l))
  nn2 = length(which(x.l == '-' & y.l != '-'))
  nn3 = length(which(x.l != '-' & y.l == '-'))
  return(nn1 - nn2 - nn3)
}
SeqDist.AA <- function(x, y) {
  if (nchar(x) != nchar(y))
    return(NA)
  x.l = unlist(strsplit(x, ''))
  y.l = unlist(strsplit(y, ''))
  tmp.vv = which(x.l != '-' & y.l != '-')
  tmp.st.x = min(which(x.l != '-'))
  tmp.ed.x = max(which(x.l != '-'))
  tmp.st.y = min(which(y.l != '-'))
  tmp.ed.y = max(which(y.l != '-'))
  gap.vv = which(x.l == '-' &
                   y.l != '-' | x.l != '-' & y.l == '-')## needs work
  gap.vv = gap.vv[gap.vv >= max(c(tmp.st.x, tmp.st.y)) &
                    gap.vv <= min(c(tmp.ed.x, tmp.ed.y))]
  tmp = sum(diag(BLOSUM62[x.l[tmp.vv], y.l[tmp.vv]]))
  tmp0 = sum(diag(BLOSUM62[x.l[tmp.vv], x.l[tmp.vv]]))
  score = (tmp0 - tmp + length(gap.vv)) / length(tmp.vv)
  return(score)
}
MergeSeqs <- function(seqs, dd) {
  
  tmp.seqs = gsub('-', '', seqs)
  tmp.vv = order(nchar(tmp.seqs))
  seqs = seqs[tmp.vv]
  dd = dd[tmp.vv, ]
  tmp.seqs = tmp.seqs[tmp.vv]
  newseqs = tmp.seqs
  seq.labs = rep(0, length(newseqs))
  nn = length(seq.labs)
  while (T) {
    for (ii in 1:nn) {
      if (length(grep(newseqs[ii], newseqs[-ii])) > 0)
        seq.labs[ii] = 1
    }
    if (sum(seq.labs) == 0)
      break
    vv0 = which(seq.labs == 0)
    newseqs = newseqs[vv0]
    seq.labs = seq.labs[vv0]
    nn = length(newseqs)
    seqs = seqs[vv0]
    dd = dd[vv0, ]
  }
  return(list(SS = seqs, DD = dd))
}
CreateMotifList <- function(mm) {
  tmp.mat = matrix(unlist(strsplit(rep(mm, 8, ''), '')), nrow = 8, byrow =
                     T)
  diag(tmp.mat) = '.'
  mm.list = apply(tmp.mat, 1, paste, collapse = '')
  return(mm.list)
}
MergeMotifs <- function(motif.list) {
  ## Merge motifs by allowing one mismatch
  unique.motifs = c()
  for (mm in motif.list) {
    mm.list = CreateMotifList(mm)
    sign = 0
    for (tmp.mm in mm.list) {
      if (length(grep(tmp.mm, unique.motifs)) > 0) {
        sign = 1
        break
      }
    }
    if (sign == 0)
      unique.motifs = c(unique.motifs, mm)
  }
  return(unique.motifs)
}

BuildBCRlineage <- function(sampleID, Bdata = BCRdata, start=3, end=10) {
  ## Given sample ID, start and end position of complete CDR3,  return all the lineages in the sample
  # sampleID <- "SRR3184301"
  # Bdata = cdr3.bcr.heavy
  # start=3
  # end=10
  tmp.dd.ss = subset(Bdata, sample == sampleID)
  tmp.dd.ss = tmp.dd.ss[!duplicated(tmp.dd.ss[,"CDR3nt"]),]
  if (is.null(dim(tmp.dd.ss)))
    return(NA)
  tmp.comp.vv <- which(tmp.dd.ss[, "is_complete"] == "Y")
  comp.CDR3.ss = data.frame(CDR3aa = tmp.dd.ss[tmp.comp.vv, "CDR3aa"])
  if (length(comp.CDR3.ss) == 0)
    return(NA)
  tmp.tt = table(substr(comp.CDR3.ss$CDR3aa, start, end))
  tmp.tt = sort(tmp.tt, decreasing = T)
  tmp.tt <- tmp.tt[which(nchar(names(tmp.tt))==(end-start+1))]
  tmp.motifs = MergeMotifs(names(tmp.tt))
  count = 0
  BCRlineage = c() ## a list of BCR lineage trees
  kept.motifs = c()
  for (mm in tmp.motifs) {
    mm.list = CreateMotifList(mm)
    tmp.vv.ss = c()
    for (tmp.mm in mm.list) {
      tmp.vv.ss = c(tmp.vv.ss, grep(tmp.mm, tmp.dd.ss$CDR3aa))
    }
    tmp.vv.ss = unique(tmp.vv.ss)
    if (length(tmp.vv.ss) < 2)
      next
    SEQs = unique(as.character(tmp.dd.ss[tmp.vv.ss, "CDR3nt"]))
    #SEQs = SEQs$CDR3nt
    tmp.dd0 = tmp.dd.ss[tmp.vv.ss, ]
    setDF(tmp.dd0)   ###format as dataframe
    rownames(tmp.dd0) = tmp.dd0$CDR3nt   ####same cdr3dna, same cdr3aa, different Ig gene and totaldna
    tmp.dd0 = tmp.dd0[SEQs, ]
    if (length(SEQs) < 3)
      next
    MSAalign = msa(DNAStringSet(SEQs), 'ClustalW')
    seqs = as.character(attributes(MSAalign)$unmasked)
    seqs0 = gsub('-', '', seqs)
    tmp.dd0 = tmp.dd0[match(seqs0, SEQs),]
    tmp = MergeSeqs(seqs, tmp.dd0)
    seqs = tmp$SS
    tmp.dd0 = tmp$DD
    if (is.null(dim(tmp.dd0)))
      next
    nn = nrow(tmp.dd0)
    if (nn <= 3)
      next
    sDist = matrix(0, nn, nn)
    for (ii in 1:nn) {
      for (jj in ii:nn) {
        if (jj == ii)
          next
        tmp.dist = SeqDist(seqs[ii], seqs[jj])
        sDist[ii, jj] = sDist[ii, jj] + tmp.dist
      }
    }
    kept.motifs = c(kept.motifs, mm)
    rownames(tmp.dd0) = NULL
    lineage.obj = list(distMat = sDist,
                       Sequences = seqs,
                       data = tmp.dd0)
    BCRlineage = c(BCRlineage, list(lineage.obj))
    count = count + 1
  }
  names(BCRlineage) = kept.motifs
  return(BCRlineage)
}
##function of computing clonality
getClonality <- function(sampleID, Bdata = BCRdata, start=3, end=10) {
  ## Given sample ID, start and end position of complete CDR3,  return all the lineages in the sample
  # sampleID <- "SRR3184301"
  # Bdata = cdr3.bcr.heavy
  # start=3
  # end=10
  cluster.ID <- c()
  cluster.list <- list()
  tmp.dd.ss = subset(Bdata, sample == sampleID)
  tmp.dd.ss = tmp.dd.ss[!duplicated(tmp.dd.ss[,"CDR3nt"]),]
  if (is.null(dim(tmp.dd.ss)))
    return(NA)
  tmp.comp.vv <- which(tmp.dd.ss[, "is_complete"] == "Y")
  comp.CDR3.ss = data.frame(CDR3aa = tmp.dd.ss[tmp.comp.vv, "CDR3aa"])
  if (length(comp.CDR3.ss) == 0)
    return(NA)
  tmp.tt = table(substr(comp.CDR3.ss$CDR3aa, start, end))
  tmp.tt = sort(tmp.tt, decreasing = T)
  tmp.tt <- tmp.tt[which(nchar(names(tmp.tt))==(end-start+1))]
  tmp.motifs = MergeMotifs(names(tmp.tt))
  count = 0
  BCRlineage = c() ## a list of BCR lineage trees
  kept.motifs = c()
  for (mm in tmp.motifs) {
    mm.list = CreateMotifList(mm)
    
    tmp.vv.ss = c()
    for (tmp.mm in mm.list) {
      tmp.vv.ss = c(tmp.vv.ss, grep(tmp.mm, tmp.dd.ss$CDR3aa))
    }
    tmp.vv.ss = unique(tmp.vv.ss)
    if (length(tmp.vv.ss) < 2)
      next
    SEQs = unique(as.character(tmp.dd.ss[tmp.vv.ss, "CDR3nt"]))
    tmp.dd0 = tmp.dd.ss[tmp.vv.ss, ]
    setDF(tmp.dd0)   ###format as dataframe
    rownames(tmp.dd0) = tmp.dd0$CDR3nt   ####same cdr3dna, same cdr3aa, different Ig gene and totaldna
    tmp.dd0 = tmp.dd0[SEQs, ]
    if (length(SEQs) < 3)
      next
    MSAalign = msa(DNAStringSet(SEQs), 'ClustalW')
    seqs = as.character(attributes(MSAalign)$unmasked)
    seqs0 = gsub('-', '', seqs)
    tmp.dd0 = tmp.dd0[match(seqs0, SEQs),]
    tmp = MergeSeqs(seqs, tmp.dd0)
    seqs = tmp$SS
    tmp.dd0 = tmp$DD
    if (is.null(dim(tmp.dd0)))
      next
    nn = nrow(tmp.dd0)
    if (nn <= 3)
      next
    cluster.ID <- c(cluster.ID,seqs0)
    cluster.freq <- cbind.data.frame(sample = sampleID, CDR3aa = mm, frequency = sum(tmp.dd0$frequency))
    cluster.list[[mm]] <- cluster.freq 
  }
  ###get non-cluster clones
  non.clonal <- subset(Bdata, sample == sampleID & !(CDR3nt %in% cluster.ID),
                       select = c(sample,CDR3aa,frequency))
  ###total cluster clones
  clonal <- do.call("rbind",cluster.list)
  ###combine
  clone.frequency <- rbind(clonal,non.clonal) %>% mutate(normalized_est_clonal_exp = frequency/sum(frequency) )
  ###calculate clonality
  entropy <- -sum(clone.frequency$normalized_est_clonal_exp*log2(clone.frequency$normalized_est_clonal_exp))
  clonality <- 1 - entropy/log2(dim(clone.frequency)[1])
  res <- c(sample = sampleID,clonality = clonality)
  print(sampleID)
  return(res)
}
##function of SHM ratio
getSHMratio <-function(tmp){
  all.dist <- 0
  mm <- 0
  for (kk in tmp){
    no.gap <- sapply(kk$Sequences,function(x) nchar(gsub("-","",x)))
    base.l <- max(no.gap)                
    dist=kk$distMat
    dist[dist > 1] <- 0
    dist.s <- sum(dist)
    all.dist <- all.dist + dist.s
    mm <- mm + base.l
    print(base.l);print(dist.s)
  }
  s.ratio <- all.dist/mm
  name <- tmp[[1]]$data$sample[1]
  return (paste0(name, ";", s.ratio))
}
##function of computing TCR clonality
getClonalityTCR <- function(sampleID, Tdata = TCRdata) {
  ## Given sample ID, start and end position of complete CDR3,  return all the lineages in the sample
  #   sampleID <- "FZ-100Aligned.sortedByCoord.out.bam"
  #   Bdata = BCRdata
  ###get clones and frequency
  clone.frequency <- subset(Tdata, sample == sampleID ,select = c(sample,CDR3aa,frequency)) %>% 
    mutate(normalized_est_clonal_exp = frequency/sum(frequency) )
  ###calculate clonality
  entropy <- -sum(clone.frequency$normalized_est_clonal_exp*log2(clone.frequency$normalized_est_clonal_exp))
  clonality <- 1 - entropy/log2(dim(clone.frequency)[1])
  res <- c(sample = sampleID,clonality = clonality)
  return(res)
}
### BCR cluster & isotype class switch in each sample
get.bcr.cluster.classswitch <- function(bcr_clusters){
  bcr.cluster.isotypes <- NULL
  for(i in 1:length(bcr_clusters)){
    tmp <- bcr_clusters[[i]]
    if(is.null(tmp)==T){
      next
    }
    if(is.na(tmp)==T){
      next
    }
    tmp.cw <- matrix(0, nrow=length(tmp), ncol=12)
    colnames(tmp.cw) <- c('filename', 'motif', 'IGHA1','IGHA2','IGHE','IGHD','IGHM','IGHG1','IGHG2','IGHG3','IGHG4', 'Unidentified')
    tmp.cw[,'filename'] <- names(bcr_clusters)[i]
    tmp.cw[,'motif'] <- names(tmp)
    for(j in 1:length(tmp)){
      # j <- 1
      tmp_cluster <- tmp[[j]]
      tmp.is <- as.character(tmp_cluster$data$C)
      id <- which(tmp.is %in% c('IGHA1','IGHA2','IGHE','IGHD','IGHM','IGHG1','IGHG2','IGHG3','IGHG4'))
      if(length(id) > 0){
        tmp.is[-id] = 'Unidentified'
      }else{
        tmp.is = rep('Unidentified', length(tmp.is))
      }
      isotype.count <- table(tmp.is)
      tmp.cw[j, names(isotype.count)]=isotype.count
    }
    bcr.cluster.isotypes <- rbind(bcr.cluster.isotypes, tmp.cw)
  }
  return(bcr.cluster.isotypes)
}
######### function of computing Jacard Similarity
getCDR3Jacard <- function(meta, cdr3.complete){
  CDR3Jacard <- lapply(rownames(meta), function(x){
    n.bcr <- subset(cdr3.complete,sample == x)
    TBCR <- lapply(rownames(meta), function(y){
      if(x != y){
        t.bcr <- subset(cdr3.complete, sample == y)
        if(dim(t.bcr)[1] == 0){
          nt <- 0
        }
        else{
          nt <- signif(length(intersect(n.bcr$CDR3aa,t.bcr$CDR3aa))/length(na.omit(unique(t.bcr$CDR3aa))),5)
        }
        share <- cbind.data.frame(s1 = x, s2 = y, jacard = nt, 
                                  tIg = paste0(intersect(n.bcr$C,t.bcr$C),collapse = ","),
                                  nIg = paste0(intersect(n.bcr$C,t.bcr$C),collapse = ","),
                                  share = paste0(intersect(n.bcr$CDR3aa,t.bcr$CDR3aa),collapse = ","))
        return (share)
      }
    })
    similarity <- do.call("rbind",TBCR)
    return (similarity)
  })
  share.mat <- do.call("rbind",CDR3Jacard)
  return(share.mat)
}

fl = list.files(pattern = "*_report.tsv")
fl2 = gsub("TRUST_", "", gsub('\\_r.*', "", fl))

mas = read.csv(paste0('/fs/ess/PAS0854/Active_projects/', name, "/Master.csv"))
con = read.csv(paste0('/fs/ess/PAS0854/Active_projects/', name, "/Controls.csv"))

fl = lapply(fl, function(x){
    tryCatch({
      read.table(x, sep="\t")
    },
    error=function(cond){
      message('No lines in: ', x)
      return(NA)
    })
  })

fl2 = fl2[-which(is.na(fl) == T)]
fl = fl[-which(is.na(fl) == T)]

for (x in 1:length(fl)){
  samp = fl2[x]
  if (grepl('SL', samp) == T){
    if (samp %in% mas$RNAseq){
      stat = mas$Type[which(mas$RNAseq == samp)]
    } else {
      stat = 'out'
    }
  } else {
    case = UUIDtoUUID(samp, to_type = 'case_id')[, 2]
    case = UUIDtoBarcode(case, from_type = 'case_id')[, 2]
    if (case %in% mas$ID){
      stat = mas$Type[which(mas$ID == case)]
    } else {
      stat = 'control'
    }
  }
  fl[[x]] = cbind(fl[[x]], rep(samp, nrow(fl[[x]])), rep(stat, nrow(fl[[x]])))
  colnames(fl[[x]]) = c('count', 'frequency', 'CDR3nt', 'CDR3aa', 'V', 'D', 'J', 'C', 'cid', 'cid_full_length', 'sample', 'clinic')
}
fl = do.call(rbind, fl)
fl = fl[fl$clinic != 'out', ]

write.csv(fl, paste0('/fs/ess/PAS0854/Active_projects/', name, '/cdr3.out.processed.csv'))

cdr3 = fl

rm(fl)

cdr3$is_complete <- sapply(cdr3$CDR3aa, function(x) 
  ifelse(x != "partial" && x != "out_of_frame" && !grepl("^_",x) && !grepl("^\\?", x),"Y","N"))

cdr3.bcr <- subset(cdr3, grepl("^IG",V) | grepl("^IG",J) | grepl("^IG",C))
cdr3.tcr <- subset(cdr3, grepl("^TR",V) | grepl("^TR",J) | grepl("^TR",C))

cdr3.bcr <- cdr3.bcr %>% mutate(lib.size = sum(count)) 
cdr3.tcr <- cdr3.tcr %>% mutate(lib.size = sum(count))

cdr3.bcr.heavy <- subset(cdr3.bcr, grepl("^IGH",V) | grepl("^IGH",J) | grepl("^IGH",C))
cdr3.bcr.light <- subset(cdr3.bcr, grepl("^IG[K|L]",V) | grepl("^IG[K|L]",J) | grepl("^IG[K|L]",C))

sample_bcr_cluster <- lapply(unique(cdr3.bcr.heavy$sample), function(sample){
  BuildBCRlineage(sampleID = sample, Bdata = cdr3.bcr.heavy, start=3, end=10)})

cpk <- aggregate(CDR3aa ~ sample+clinic+lib.size, cdr3.tcr, function(x) length(unique(x))) %>%
  mutate(CPK = signif(CDR3aa/(lib.size/1000),4))

single_sample_bcr_clonality <- lapply(unique(cdr3.bcr.heavy$sample), function(sample){
  getClonality(sample, cdr3.bcr.heavy, start=3, end=10)})
single_sample_bcr_clonality <- do.call(rbind, single_sample_bcr_clonality)
single_sample_bcr_clonality <- single_sample_bcr_clonality[single_sample_bcr_clonality[, 2] != "NaN", ]

single_sample_tcr_clonality <- lapply(unique(cdr3.tcr$sample), function(sample){
  getClonalityTCR(sample,cdr3.tcr)})
single_sample_tcr_clonality <- do.call(rbind, single_sample_tcr_clonality)
single_sample_tcr_clonality <- single_sample_tcr_clonality[single_sample_tcr_clonality[, 2] != "NaN", ]

SHM.ratio <- do.call(rbind, lapply(sample_bcr_cluster, function(x) getSHMratio(x)))
SHM.ratio <- matrix(unlist(strsplit(SHM.ratio, ";")), ncol = 2, byrow = T)
SHM.ratio <- SHM.ratio[SHM.ratio[, 1] != "", ]

st.Ig <- cdr3.bcr.heavy %>% 
  group_by(clinic, sample) %>%
  mutate(est_clonal_exp_norm = frequency/sum(frequency)) %>%   #as.numeric(sample.clones[filename,2])
  dplyr::filter(C != "*" & C !=".") %>%
  group_by(sample, C, clinic) %>% 
  dplyr::summarise(Num.Ig = sum(est_clonal_exp_norm))

st.Ig = as.data.frame(st.Ig)
single_sample_bcr_clonality = as.data.frame(single_sample_bcr_clonality)
single_sample_tcr_clonality = as.data.frame(single_sample_tcr_clonality)
SHM.ratio = as.data.frame(SHM.ratio)
cpk = as.data.frame(cpk)

samps = c(SHM.ratio[, 1], single_sample_bcr_clonality$sample, single_sample_tcr_clonality$sample)
samps = unique(samps)
sampst = samps[!grepl('SL', samps)]
sampso = samps[grepl('SL', samps)]

sampst = UUIDtoBarcode(sampst, from_type = 'file_id')

sampst$associated_entities.entity_submitter_id = unlist(lapply(strsplit(sampst$associated_entities.entity_submitter_id, "-"), 
                                                        function(x) return(paste0(x[1:3], collapse = "-"))))
sampso = cbind(sampso, unlist(lapply(1:length(sampso), function(x) mas$Type[mas$RNAseq == sampso[x]])))

sampst = cbind(as.character(sampst$file_id), unlist(lapply(1:length(sampst), function(x){
  if (sampst$associated_entities.entity_submitter_id[x] %in% mas$ID){
    return(mas$Type[which(mas$ID == sampst$associated_entities.entity_submitter_id[x])])
  } else {
    return(con$Type[which(con$ID == sampst$associated_entities.entity_submitter_id[x])])
  }
})))

colnames(sampst) = c('Sample', 'Type')
colnames(sampso) = c('Sample', 'Type')

key = rbind('Sample', 'Type')

colnames(SHM.ratio) = c('sample', 'SHM.ratio')
new = c()
for (x in SHM.ratio$sample){
  new = c(new, 
          key[which(key[, 1] == x), 2])
}
SHM.ratio$clinic = new

single_sample_bcr_clonality = single_sample_bcr_clonality[single_sample_bcr_clonality$sample %in% 
                                                            key[, 1], ]
single_sample_tcr_clonality = single_sample_tcr_clonality[single_sample_tcr_clonality$sample %in% 
                                                            key[, 1], ]
new = c()
for (x in single_sample_bcr_clonality$sample){
  new = c(new, 
          key[which(key[, 1] == x), 2])
}
single_sample_bcr_clonality$clinic = new

new = c()
for (x in single_sample_tcr_clonality$sample){
  new = c(new, 
          key[which(key[, 1] == x), 2])
}
single_sample_tcr_clonality$clinic = new

setwd(paste0('/fs/ess/PAS0854/Active_projects/', name))

SHM.ratio$SHM.ratio = as.numeric(as.character(SHM.ratio$SHM.ratio))

pval = kruskal.test(SHM.ratio ~ clinic, data = SHM.ratio)

SHM <- ggplot(SHM.ratio, aes(x = clinic, y = SHM.ratio)) + 
  geom_boxplot(aes(fill = clinic)) +
  annotate(geom="text", x=3, y=1.05, label=paste0("Kruskal Wallis: ", round(pval$p.value, 2)),
             color="red") +
  theme(legend.position = 'none')

single_sample_tcr_clonality$clonality = as.numeric(as.character(single_sample_tcr_clonality$clonality))

pval = kruskal.test(clonality ~ clinic, data = single_sample_tcr_clonality)

tcr <- ggplot(single_sample_tcr_clonality, aes(x = clinic, y = clonality)) + 
  geom_boxplot(aes(fill = clinic)) +
  annotate(geom="text", x=3, y=1.05, label=paste0("Kruskal Wallis: ", round(pval$p.value, 2)),
           color="red") + 
  theme(legend.position = 'none')

single_sample_bcr_clonality$clonality = as.numeric(as.character(single_sample_bcr_clonality$clonality))

pval = kruskal.test(clonality ~ clinic, data = single_sample_bcr_clonality)

bcr <- ggplot(single_sample_bcr_clonality, aes(x = clinic, y = clonality)) + 
  geom_boxplot(aes(fill = clinic)) +
  annotate(geom="text", x=3, y=1.05, label=paste0("Kruskal Wallis: ", round(pval$p.value, 2)),
           color="red") + 
  theme(legend.position = 'none')

st.Ig$Num.Ig = as.numeric(as.character(st.Ig$Num.Ig))

pvals = lapply(unique(st.Ig$C), function(x){
  new = subset(st.Ig, C == x)
  pval = kruskal.test(Num.Ig ~ clinic, data = new)
  return(paste0(x, ";", pval$p.value))
})

pvals = unlist(pvals)
pvals = unlist(strsplit(pvals, ";"))
pvals = matrix(pvals, ncol = 2, byrow = T)
colnames(pvals) = c('C', 'p.value')
pvals = as.data.frame(pvals)
pvals$p.value = paste0('Kruskal Wallis: ', round(as.numeric(as.character(pvals$p.value)), 2))

Ig <- ggplot(st.Ig, aes(x = clinic, y = Num.Ig)) + 
  geom_boxplot(aes(fill = clinic)) +
  facet_wrap(~ C, scales = 'free') + 
  theme(legend.position = 'none')
Ig <- Ig + geom_text(
  data    = pvals,
  mapping = aes(x = -Inf, y = -Inf, label = p.value),
  hjust   = -1,
  vjust   = -15
)

pval = kruskal.test(clonality ~ clinic, data = single_sample_tcr_clonality)

cpkg <- ggplot(cpk, aes(x = clinic, y = CPK)) + 
  geom_boxplot(aes(fill = clinic)) +
  theme(legend.position = 'none') +
  annotate(geom="text", x=3, y=20, label=paste0("Kruskal Wallis: ", round(pval$p.value, 2)),
           color="red")
#save plots somehow

write.csv(st.Ig, 'IgPercent.csv')
write.csv(single_sample_bcr_clonality, 'BCRclon.csv')
write.csv(single_sample_tcr_clonality, 'TCRclon.csv')
write.csv(SHM.ratio, 'SHM.csv')
write.csv(cpk, 'CPK.csv')

