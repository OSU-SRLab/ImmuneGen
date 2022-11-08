library('TxDb.Hsapiens.UCSC.hg38.knownGene')
library('DEGobj.utils')
library('quantiseqr')
library("AnnotationDbi")
library("org.Hs.eg.db")

rna = read.csv('RNAseq.csv', row.names = 1)

tx = TxDb.Hsapiens.UCSC.hg38.knownGene

genelen = transcripts(tx, columns=c("tx_id", "tx_name", "gene_id"))

genelen$symbol = mapIds(org.Hs.eg.db,
                        keys=as.character(genelen$gene_id), 
                        column="SYMBOL",
                        keytype="ENTREZID",
                        multiVals="first")
genelen2 = as.data.frame(genelen)
genelen2 = genelen2[!(is.na(names(genelen2$symbol))), ]
genelen2 = genelen2[genelen2$symbol %in% rownames(rna), ]

ma = matrix(NA, ncol=2, nrow = 1)

for (x in unique(genelen2$symbol)){
  set = genelen2[genelen2$symbol == x, ]
  ma = rbind(ma, cbind(x, max(set$width)))
}

ma = ma[-1, ]
ma = as.data.frame(ma)
colnames(ma) = c('Gene_Symbol', 'Length')
ma$Length = as.numeric(as.character(ma$Length))
rna = rna[rownames(rna) %in% ma$Gene_Symbol, ]
ma = ma[match(rownames(rna), ma$Gene_Symbol), ]

rnatmp = convertCounts(as.matrix(rna), 'TPM', ma$Length)

rnatmp = as.data.frame(rnatmp)

ti_racle <- quantiseqr::run_quantiseq(
  expression_data = rnatmp,
  signature_matrix = "TIL10",
  is_arraydata = FALSE,
  is_tumordata = TRUE,
  scale_mRNA = TRUE
)

key <- read.csv("LOF.GOF.key.csv")
mut <- read.csv('mut.stat.key.csv')

key2 = cbind(key$V1, key$func2)
colnames(key2) = c('V1', 'V2')

coldata = rbind(key2, cbind(mut$V1[mut$V2 == 'control'], 
                            mut$V2[mut$V2 == 'control']))
coldata = as.data.frame(coldata)

coldata = coldata[gsub('\\-', ".", coldata[, 1]) %in% ti_racle$Sample, ]
coldata = coldata[match(ti_racle$Sample, gsub('\\-', ".", coldata[, 1])), ]
quant = cbind(ti_racle[2:ncol(ti_racle)], coldata$V2)
colnames(quant)[ncol(quant)] = 'Function'
quantm = melt(quant, id.vars = c('Function'))
quantm$value = as.numeric(as.character(quantm$value))

agquant = aggregate(. ~ variable + Function, quantm, mean)

maov = manova(as.matrix(quant[, 2:(ncol(quant)-1)]) ~ quant$Function, data = quant)
pval = summary(maov)
pval = pval$stats[1, 6]

if (pval < .05){
  ma = matrix(NA, 1, 3)
  for (x in unique(quantm$variable)){
    set = quant[, colnames(quant) %in% c(x, 'Function')]
    colnames(set) = c('Var', "Fun")
    set = as.data.frame(set)
    set$Var = as.numeric(as.character(set$Var))
    k = kruskal.test(set[, 1], set[, 2])
    k = k$p.value
    if (k < .05){
      post.hoc = c()
      p = pairwise.wilcox.test(set[, 1], set[, 2])
      p = p$p.value
      id = which(colnames(p) == 'control')
      p = p[, id]
      out = names(p)[which(p < .05)]
      ag = aggregate(. ~ Fun, set, mean)
      for (y in out){
        if (ag[which(ag[, 1] == y), 2] > ag[which(ag[, 1] == 'control'), 2]){
          post.hoc = c(post.hoc, paste0(y, ':greater'))
        } else {
          post.hoc = c(post.hoc, paste0(y, ':lesser'))
        }
      }
      post.hoc = paste0(post.hoc, collapse = ";")
    } else {
      post.hoc = ""
    }
    ma = rbind(ma, cbind(x, k, post.hoc))
  }
}
ma = ma[-1,]
ma = ma[ma[, 3] != "", ]

l = list()
for (x in 1:length(ma[, 1])){
  my_comparisons <- list(  c("control", "unknown"), c("control", "LOF"),c("control", "GOF")  )
  new = quantm[quantm$variable == ma[x, 1], ]
  mx = max(new[, 3])
  l[[x]]= ggplot(new, aes(x = Function, y= value)) + 
    geom_violin(aes(fill = Function)) + 
    facet_wrap(~ variable) +
    stat_compare_means(comparisons = my_comparisons, label.y = c(mx, mx + .05, mx + .1)) +
    stat_compare_means(label.y = mx + .2) + theme_bw()
}
mx = length(l) + 1

l[[mx]] = ggplot(agquant, aes(fill = variable, x = Function, y = value)) +
  geom_bar(position="stack", stat="identity") + 
  annotate(geom="text", x=1, y=1.05, label=paste0('Manova: ', round(pval, 3))) + theme_bw()

pdf('QuantiseqRes.pdf')

for (x in l){
  print(x)
}

dev.off()
