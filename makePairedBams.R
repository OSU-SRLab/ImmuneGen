args = commandArgs(trailingOnly = T)
nameOG = as.character(args[1])

setwd(paste0(
  '/fs/ess/PAS0854/Active_projects/',
  nameOG,
  '/bams'
))

fl = list.files(pattern=".bam")
fl = fl[grepl('TCGA', fl)]
fl = fl[!grepl('bai', fl)]
fls = strsplit(fl, "-")
fls = do.call(rbind, fls)
fls = fls[, 1:3]
fls = apply(fls, 1, function(x){ return(paste0(x, collapse="-"))})
fls = unique(fls)
st <- read.csv("/fs/ess/PAS0854/Raven/Immune_Checkpoints/TCGA_SampleType_Code.csv", header=T)
st = as.data.frame(st)
st$Code = paste0('0', as.character(st$Code))
st$Code = lapply(st$Code, function(x){ if(nchar(x) == 3){ return(substring(x, 2, 3))} else { return(x)}})
out = st$Code[grep('normal', st$Definition, ignore.case = T)]
l = c()
for (x in fls){
  name = x
  set = fl[grep(x, fl)]
  if (length(set) == 2){
    l = c(l, paste0(c(name, set), collapse=";"))
  } else if (length(set) > 1) {
    set1 = strsplit(set, "-")
    set1 = do.call(rbind, set1)
    set1 = as.data.frame(set1)
    set1$Type = grepl(paste0(unlist(out), collapse="|"), set1[, 4])
    n = set[which(set1$Type == 'TRUE')]
    t = set[which(set1$Type == 'FALSE')]
    set = expand.grid(t, n)
    set = as.matrix(set)
    name1 = as.matrix(set)
    name1 = cbind(unlist(lapply(name1[, 1], function(r) strsplit(r, "-")[[1]][4])),
                         unlist(lapply(name1[, 2], function(r) strsplit(r, "-")[[1]][4])),
                         name)
    name1 = unlist(lapply(1:nrow(name1), function(v) paste0(name1[v, ], collapse="-")))
    l = c(l, unlist(lapply(1:nrow(set), function(y) paste0(c(name1[y], set[y, ]), collapse =";"))))
  }
}

setwd(paste0(
  '/fs/ess/PAS0854/Active_projects/',
  nameOG))

fl = file('CD40.PairedBams.txt')
writeLines(l, fl)
close(fl)
