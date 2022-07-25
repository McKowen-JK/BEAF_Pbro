GeneList = read.csv("/home/chart/scripts/gene_lists/finalsetMJG_BedFormat_LacZ_50bpPause.csv", header=T)
OutputDirectory = '/home/chart/PROseq2/DESeq/DESeq_gene_lists/'

#Use the get counts function from when I plotted replicates.  I'll run DESeq on these same data (counts of total reads in the gene body, per gene).  USE THE NON-NORMALIZED BIGWIGS
#2016oct01 I don't know why this should count the promter region as opposed to the entire gene (with offsets set at 0)
library(bigWig)
library(DESeq2)
basepath = '/home/chart/PROseq2/bigwig/'    #I made non-normalized bigWigs for running DESeq

tss.off = 0
geneEnd.off = 0

# PromReads.tss.off = -150
# PromReads.geneEnd.off = 150

getCounts <- function(wig.p, wig.m, genes, tss.off, geneEnd.off) {
  N = dim(genes)[1]
  plusStrand = genes[,6] == '+'    #vector of all genes returns true if the gene is on the plus strand, false if on the minus strand.  6 is column with strand information in Jay's gene list  [,6] means all of column 6
  
  result = vector(mode="integer", length=N)
  
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 2])
    end = as.integer(genes[i, 3])
    
    if (plusStrand[i]) {  #if the gene is on the plus strand...
      wig = wig.p    #use the plus strand bigWig file
      
      qStart = start + tss.off    #start from the genelist plus the offset defined at the beginning
      qEnd = end - geneEnd.off   #for Jay's list, we don't want the end of genes (his goes all the way to the end of the gene annotation), so define the end as 200 bp before the end 
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(wig, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))  #if data was returned (is not null), add the sum of all the reads to the results list defined at the beginning (data is a table, third row is the read counts at each base)
    } else {
      wig = wig.m
      if (is.null(wig.m))
        wig = wig.p
      
      qStart = start + geneEnd.off
      qEnd = end - tss.off
      
      if (qStart > qEnd) {
        result[i] = NA
        next
      }
      
      data = query.bigWig(wig, chrom, qStart, qEnd)
      
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    }
  }
  
  return(cbind(result))
}#returns the number of counts in gene body

wig.table = rbind(
  c("BF1_plus_noMito.bw", "BF1_minus_noMito.bw", "BEAF_1"),
  c("BF2_plus_noMito.bw", "BF2_minus_noMito.bw", "BEAF_2"),
  c("PB1_plus_noMito.bw", "PB1_minus_noMito.bw", "PBro_1"),
  c("PB2_plus_noMito.bw", "PB2_minus_noMito.bw", "PBro_2"),
  c("LZ1_plus_noMito.bw", "LZ1_minus_noMito.bw", "LacZ_1"),
  c("LZ2_plus_noMito.bw", "LZ2_minus_noMito.bw", "LacZ_2"))

load.wigset <- function(row) {
  file = wig.table[row, 1]
  wig.p = NULL
  if (file != "")
    wig.p = load.bigWig(paste(basepath, file, sep=''))
  file = wig.table[row, 2]
  wig.m = NULL
  if (file != "")
    wig.m = load.bigWig(paste(basepath, file, sep=''))
  
  return(list(wig.p, wig.m, wig.table[row, 3]))
}

unload.wigset <- function(set) {
  if (!is.null(set[[1]]))
    unload.bigWig(set[[1]])
  if (!is.null(set[[2]]))
    unload.bigWig(set[[2]])
}

N = dim(wig.table)[1]
gb.res = vector(mode="list", length=N)
for (i in 1:N) {
  cat("* loading", i, "\n")
  wigs = load.wigset(i)#this should give all the wig files
  cat("* computing ...\n")
  gb.res[[i]] = getCounts(wigs[[1]], wigs[[2]], GeneList, tss.off, geneEnd.off)
  colnames(gb.res[[i]]) <- c(paste(wigs[[3]], "counts", sep="_"))
  cat("* unloading.\n")
  unload.wigset(wigs)
}

output <- cbind(GeneList, gb.res[[1]], gb.res[[2]], gb.res[[3]], gb.res[[4]], gb.res[[5]], gb.res[[6]])

BF_LZ_DESeq = output[,c(7:8,11:12)]
rownames(BF_LZ_DESeq) <- output[,4]

coldata = data.frame(row.names=colnames(BF_LZ_DESeq),condition=as.factor(c(rep("BEAF",2),rep("LacZ",2))))    #I don't understand this... just copy from Fabs

dds.BF_LZ_DESeq <- DESeqDataSetFromMatrix(
  countData = BF_LZ_DESeq,
  colData = coldata,
  design = ~ condition)   #build the DESeq matrix.  Uses the metadata from coldata and the counts from BF_LZ_DESeq.  Design lets it know that the condition column in the coldata frame indicate which are the treated samples (points to the column in coldata that has the options for treatment... I only have experimental or control RNAi, but there could be multiple variables like time, concentration, and temperature, etc)

dds.BF_LZ_DESeq$condition <- relevel(dds.BF_LZ_DESeq$condition, "LacZ")  #makes sure that the control samples are first (so that fold change is right when DESeq calculates it).  Probably unnecessary as it is already in this order.

DESeq.BF_LZ_DESeq=DESeq(dds.BF_LZ_DESeq)
res.BF_LZ_DESeq = results(DESeq.BF_LZ_DESeq)
write.csv(as.data.frame(res.BF_LZ_DESeq), file=paste(c(OutputDirectory,'aBF_LZ_50bpPause_DESeq2.csv'), collapse = ''))
write.csv(BF_LZ_DESeq, '/home/chart/PROseq2/DESeq/DESeq_gene_lists/bf_lz_pr_counts')

plotMA(DESeq.BF_LZ_DESeq)
title(main='Effect of BEAF RNAi', cex=1.5)


PB_LZ_DESeq = output[,c(9:10,11:12)]
rownames(PB_LZ_DESeq) <- output[,4]

coldata_PB_LZ = data.frame(row.names=colnames(PB_LZ_DESeq),condition=as.factor(c(rep("PBro",2),rep("LacZ",2))))    #I don't understand this... just copy from Fabs

dds.PB_LZ_DESeq <- DESeqDataSetFromMatrix(
  countData = PB_LZ_DESeq,
  colData = coldata_PB_LZ,
  design = ~ condition)   #build the DESeq matrix.  Uses the metadata from coldata and the counts from PB_LZ_DESeq.  Design lets it know that the condition column in the coldata frame indicate which are the treated samples (points to the column in coldata that has the options for treatment... I only have experimental or control RNAi, but there could be multiple variables like time, concentration, and temperature, etc)

dds.PB_LZ_DESeq$condition <- relevel(dds.PB_LZ_DESeq$condition, "LacZ")  #makes sure that the control samples are first (so that fold change is right when DESeq calculates it).  Probably unnecessary as it is already in this order.

DESeq.PB_LZ_DESeq=DESeq(dds.PB_LZ_DESeq)
res.PB_LZ_DESeq = results(DESeq.PB_LZ_DESeq)

write.csv(as.data.frame(res.PB_LZ_DESeq), file=paste(c(OutputDirectory,'aPB_LZ_50bpPause_DESeq2.csv'), collapse = ''))
write.csv(PB_LZ_DESeq, '/home/chart/PROseq2/DESeq/DESeq_gene_lists/pb_lz_pr_counts')


plotMA(DESeq.PB_LZ_DESeq)
title(main='Effect of PBro RNAi', cex=1.5)

#################################################Now do DESeq2 to look for changes between BEAF RNAi and polybromo RNAi treated cells


BF_PB_DESeq = output[,c(7:8, 9:10)]
rownames(BF_PB_DESeq) <- output[,4]

coldata_BF_PB = data.frame(row.names=colnames(BF_PB_DESeq),condition=as.factor(c(rep("BF",2),rep("PB",2))))    #I don't understand this... just copy from Fabs

dds.BF_PB_DESeq <- DESeqDataSetFromMatrix(
  countData = BF_PB_DESeq,
  colData = coldata_BF_PB,
  design = ~ condition)   #build the DESeq matrix.  Uses the metadata from coldata and the counts from BF_PB_DESeq.  Design lets it know that the condition column in the coldata frame indicate which are the treated samples (points to the column in coldata that has the options for treatment... I only have BF RNAi or PB RNAi, but there could be multiple variables like time, concentration, and temperature, etc)

dds.BF_PB_DESeq$condition <- relevel(dds.BF_PB_DESeq$condition, "PBro RNAi")  #makes sure that the PB RNAi samples are first (so that fold change is right when DESeq calculates it).  Probably unnecessary as it is already in this order.

DESeq.BF_PB_DESeq=DESeq(dds.BF_PB_DESeq)
res.BF_PB_DESeq = results(DESeq.BF_PB_DESeq)

write.csv(as.data.frame(res.BF_PB_DESeq), file=paste(c(OutputDirectory,'aBF_PB_50bpPause_DESeq2.csv'), collapse = ''))


plotMA(DESeq.BF_PB_DESeq)
title(main='BEAF RNAi vs PBro RNAi', cex=1.5)

