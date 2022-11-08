#read in all the files, the format is wrong at first, need to fix it 
#and write the files out with the right format for future use.

args = commandArgs(trailingOnly=TRUE)

setwd(paste0('/fs/project/PAS0854/Active_projects/', args[1]))

ma = readLines('MantaAllTog.txt')
ma = lapply(ma, function(x) strsplit(x, "\t"))
ma <- data.frame(matrix(unlist(ma), nrow=length(ma), byrow=TRUE))
colnames(ma) = ma[1,]
ma = ma[-1,]
write.table(ma, 'MantaAllTog.txt', sep="\t")

gr = readLines('GridssAllTog.txt')
gr = lapply(gr, function(x) strsplit(x, "\t"))
gr <- data.frame(matrix(unlist(gr), nrow=length(gr), byrow=TRUE))
colnames(gr) = gr[1,]
gr = gr[-1,]
write.table(gr, 'GridssAllTog.txt', sep="\t")

de = readLines('DellyAllTog.txt')
de = lapply(de, function(x) strsplit(x, "\t"))
de <- data.frame(matrix(unlist(de), nrow=length(de), byrow=TRUE))
colnames(de) = de[1,]
de = de[-1,]
write.table(de, 'DellyAllTog.txt', sep="\t")

#only take the passes from gridss and delly

dep = de[de$FILTER == 'PASS', ]
grp = gr[gr$FILTER == 'PASS',]

#make an ID with the sample and the chromosome

grp$sum = paste0(grp$Sample, grp$`#CHROM`)
dep$sum = paste0(dep$Sample, dep$`#CHROM`)
ma$sum = paste0(gsub("\\..*", "", ma$Sample),ma$`#CHROM`)

#make sure the positions are numeric values

grp$POS = as.numeric(as.character(grp$POS))
dep$POS = as.numeric(as.character(dep$POS))
ma$POS = as.numeric(as.character(ma$POS))

#make a dataframe with the same samples

dem = dep[dep$sum %in% grp$sum & dep$sum %in% ma$sum,]

#go through the same sample data frames and then make a matrix 
#with all the positions in that dataframe by sample ID

r= matrix(NA, nrow = 1, ncol = 13)
colnames(r) = colnames(dem)
for (x in 1:nrow(dem)){
  gset = grp[grp$sum == dem$sum[x], ]
  mset = ma[ma$sum == dem$sum[x], ]
#make data frames with the same sample id as dem row
  set = expand.grid(gset$POS, mset$POS)
#make a dataframe with all the positions
  num = as.numeric(as.character(set[1, ])) - as.numeric(as.character(set[2, ]))
  set = set[which(abs(num) < 200), ]
  num = unlist(set)
#store the positions that are within 200 bp of eachother
  ref = as.numeric(dem$POS[x])
  dif = num - ref
  dif = abs(dif)
  dif = dif < 200
#check the position of the row to see if is within 200 bp 
#of at least one of the positions found before.
  if (TRUE %in% dif){
    r = rbind(r, dem[x, ])
  }
#put it in a data frame
}
r = as.data.frame(r)
r = r[-1, ]
if (nrow(r) >0){
    r$sum = 'gridss&delly&manta'   
}

#same as above but only between GRP and DEP

a= matrix(NA, nrow = 1, ncol = 13)
colnames(a) = colnames(dep)
for (x in 1:nrow(dep)){
  set = grp[grp$sum == dep$sum[x], ]
  num = as.numeric(as.character(set$POS))
  ref = as.numeric(dep$POS[x])
  dif = num - ref
  dif = abs(dif)
  dif = dif < 200
  if (TRUE %in% dif){
    a = rbind(a, dep[x, ])
  }
}
a = as.data.frame(a)
if (nrow(a) >0){
    a$sum = 'delly&gridss'   
}

#same as above but only between DEP and MA

t= matrix(NA, nrow = 1, ncol = 13)
colnames(t) = colnames(dep)
for (x in 1:nrow(dep)){
  set = ma[ma$sum == dep$sum[x], ]
  num = as.numeric(as.character(set$POS))
  ref = as.numeric(dep$POS[x])
  dif = num - ref
  dif = abs(dif)
  dif = dif < 200
  if (TRUE %in% dif){
    t = rbind(t, dep[x, ])
  }
}
t = as.data.frame(t)
if (nrow(r) >0){
    r$sum = 'delly&manta'   
}

#put them all together and take out duplicates

alltog = rbind(r, a, t)

alltog = alltog[!duplicated(rownames(alltog)), ]

#filter by Median mapping quality of paired-ends > 20
fil = unlist(lapply(1:nrow(alltog), function(x) return(as.numeric(as.character(strsplit(strsplit(alltog$INFO[x], ";")[[1]][6], "=")[[1]][2])) >20 )))
alltog = alltog[fil, ]

#filter by Raw high-quality read counts for the SV > 5
fil = unlist(lapply(1:nrow(alltog), function(x) return(as.numeric(as.character(strsplit(alltog[x, 11], ":")[[1]][12])) > 5)))
alltog = alltog[fil, ]

write.table(alltog, 'Imprecise.G.D.M.res.csv', sep=",")

#filter for only precise SVs
fil = unlist(lapply(1:nrow(alltog), function(x) return(strsplit(alltog$INFO[x], ";")[[1]][1] == "PRECISE")))
alltog = alltog[fil, ]

write.table(alltog, 'G.D.M.res.csv', sep=",")