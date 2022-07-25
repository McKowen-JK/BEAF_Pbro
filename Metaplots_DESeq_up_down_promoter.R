library(bigWig)

setwd('/home/chart/metaplots/DESeq_genes')

Genes_BF_pr = read.csv('/home/chart/PROseq2/DESeq/DESeq_gene_lists/BF_LZ_50bpPause_DESeq2.csv', header=TRUE)
Genes_BF_gb = read.csv('/home/chart/PROseq2/DESeq/DESeq_gene_lists/BF_LZ_GeneBody_DESeq2.csv', header=TRUE)
Genes_PB_pr = read.csv('/home/chart/PROseq2/DESeq/DESeq_gene_lists/PB_LZ_50bpPause_DESeq2.csv', header=TRUE)
Genes_PB_gb = read.csv('/home/chart/PROseq2/DESeq/DESeq_gene_lists/PB_LZ_GeneBody_DESeq2.csv', header=TRUE)
GeneList_MJG = read.table('/home/chart/scripts/gene_lists/finalsetMJG_BedFormat.dat', header=TRUE)

Genes_BF_pr_up = subset(Genes_BF_pr, padj < 0.05 & log2FoldChange > 0)
Genes_BF_gb_up = subset(Genes_BF_gb, padj < 0.05 & log2FoldChange > 0)
Genes_PB_pr_up = subset(Genes_PB_pr, padj < 0.05 & log2FoldChange > 0)
Genes_PB_gb_up = subset(Genes_PB_gb, padj < 0.05 & log2FoldChange > 0)

Genes_BF_pr_down = subset(Genes_BF_pr, padj < 0.05 & log2FoldChange < 0)
Genes_BF_gb_down = subset(Genes_BF_gb, padj < 0.05 & log2FoldChange < 0)
Genes_PB_pr_down = subset(Genes_PB_pr, padj < 0.05 & log2FoldChange < 0)
Genes_PB_gb_down = subset(Genes_PB_gb, padj < 0.05 & log2FoldChange < 0)

Genes_BF_pr_up = merge(GeneList_MJG, Genes_BF_pr_up, by.x = 'name', by.y = 'X')[,c(2:5,1,6,12,7:11)]
Genes_BF_gb_up = merge(GeneList_MJG, Genes_BF_gb_up, by.x = 'name', by.y = 'X')[,c(2:5,1,6,12,7:11)]
Genes_PB_pr_up = merge(GeneList_MJG, Genes_PB_pr_up, by.x = 'name', by.y = 'X')[,c(2:5,1,6,12,7:11)]
Genes_PB_gb_up = merge(GeneList_MJG, Genes_PB_gb_up, by.x = 'name', by.y = 'X')[,c(2:5,1,6,12,7:11)]

Genes_BF_pr_down = merge(GeneList_MJG, Genes_BF_pr_down, by.x = 'name', by.y = 'X')[,c(2:5,1,6,12,7:11)]
Genes_BF_gb_down = merge(GeneList_MJG, Genes_BF_gb_down, by.x = 'name', by.y = 'X')[,c(2:5,1,6,12,7:11)]
Genes_PB_pr_down = merge(GeneList_MJG, Genes_PB_pr_down, by.x = 'name', by.y = 'X')[,c(2:5,1,6,12,7:11)]
Genes_PB_gb_down = merge(GeneList_MJG, Genes_PB_gb_down, by.x = 'name', by.y = 'X')[,c(2:5,1,6,12,7:11)]

#if GeneList columns need to be rearranged or deleted, here is an example
#GenesWithNBE_Downstream = GenesWithNBE_Downstream[, c(3,2,4:14)]

LZ_12_plus = load.bigWig('/home/chart/PROseq2/bigwig/LacZ_12_plus_noMito.bw')
LZ_12_minus = load.bigWig('/home/chart/PROseq2/bigwig/LacZ_12_minus_noMito.bw')

BF_12_plus = load.bigWig('/home/chart/PROseq2/bigwig/BEAF_12_plus_noMito.bw')
BF_12_minus = load.bigWig('/home/chart/PROseq2/bigwig/BEAF_12_minus_noMito.bw')

PB_12_plus = load.bigWig('/home/chart/PROseq2/bigwig/PBro_12_plus_noMito.bw')
PB_12_minus = load.bigWig('/home/chart/PROseq2/bigwig/PBro_12_minus_noMito.bw')

#copy Leighton's powerpoint's code... figure out what it does after it works

meta.overlay <- function(result, step, col1=rgb(0, 0.5, 0), col2=rgb(0.1, 0, 0.2), strand = "")
{
  N=length(result[[4]])
  x=((1:N)-N/2)*step
  
  if(strand == "+")
    #draw shade area
  {polygon(c(x, rev(x)), c(result[[2]], rev(result[[3]])), col=col2, border=NA)
    
    #redraw main plot line at eop
    lines(x, result[[4]], col=col1, lwd=3)}
  else if(strand=="-")
    #draw shade area
  {polygon(c(x, rev(x)), c(-result[[2]], rev(-result[[3]])), col=col2, border=NA)
    
    #redraw main plot line on top
    
    lines(x, -result[[4]], col=col1, lwd=3)}
  else
  {print("ERROR: need strand")}
}

#####add meta.subsample because default sampleFrac is 0.1, but I have <20 genes in some lists 
##### but the subsample needs to be at least two genes. 
##### So sampleFrac needs to be increased (try 0.5 or 0.2)

meta.subsample <- function (bed, bigWig.plus, bigWig.minus, halfWindow, step, name = "Untitled", at.TSS = TRUE, do.sum = TRUE) {
#meta.subsample <- function (bed, bigWig.plus, bigWig.minus, halfWindow, step, name, at.TSS = TRUE) {
  N = dim(bed)[1]
  nPermut = 1000
  sampleFrac = 0.2
  windowSize = (2 * halfWindow) %/% step
  result = matrix(nrow = nPermut, ncol = windowSize)
  M = as.integer(round(N * sampleFrac, 0))

    values = collect.many(bed, bigWig.plus, bigWig.minus, halfWindow, step, at.TSS = at.TSS, do.sum = do.sum)
  for (i in 1:nPermut) {
    idx <- sample(N, size = M, replace = T)
    result[i, ] = colSums(values[idx, ])
  }
#changed from ci9 = 0.875 to ci9 = 0.5 and ci1 = 0.125 to ci1 = 0.5 to eliminate ci plot
  ci9 = sapply(1:windowSize, function(idx) quantile(result[, idx], 0.5))
  ci1 = sapply(1:windowSize, function(idx) quantile(result[, idx], 0.5))
  ci5 = sapply(1:windowSize, function(idx) quantile(result[, idx], 0.5))
  return(list(result, ci9, ci1, ci5, name, values))
}

##########################################################################Start plotting
BinSize = 3
halfWindow = 200

BF_pr_up_anti_meta <- meta.subsample(Genes_BF_pr_up, BF_12_minus, BF_12_plus,halfWindow,BinSize, at.TSS=T)
BF_pr_up_sense_meta <- meta.subsample(Genes_BF_pr_up, BF_12_plus, BF_12_minus,halfWindow,BinSize, at.TSS=T)

BF_gb_up_anti_meta <- meta.subsample(Genes_BF_gb_up, BF_12_minus, BF_12_plus,halfWindow,BinSize, at.TSS=T)
BF_gb_up_sense_meta <- meta.subsample(Genes_BF_gb_up, BF_12_plus, BF_12_minus,halfWindow,BinSize, at.TSS=T)

PB_pr_up_anti_meta <- meta.subsample(Genes_PB_pr_up, PB_12_minus, PB_12_plus,halfWindow,BinSize, at.TSS=T)
PB_pr_up_sense_meta <- meta.subsample(Genes_PB_pr_up, PB_12_plus, PB_12_minus,halfWindow,BinSize, at.TSS=T)

PB_gb_up_anti_meta <- meta.subsample(Genes_PB_gb_up, PB_12_minus, PB_12_plus,halfWindow,BinSize, at.TSS=T)
PB_gb_up_sense_meta <- meta.subsample(Genes_PB_gb_up, PB_12_plus, PB_12_minus,halfWindow,BinSize, at.TSS=T)

LZBF_pr_up_anti_meta <- meta.subsample(Genes_BF_pr_up, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZBF_pr_up_sense_meta <- meta.subsample(Genes_BF_pr_up, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

LZBF_gb_up_anti_meta <- meta.subsample(Genes_BF_gb_up, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZBF_gb_up_sense_meta <- meta.subsample(Genes_BF_gb_up, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

LZPB_pr_up_anti_meta <- meta.subsample(Genes_PB_pr_up, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZPB_pr_up_sense_meta <- meta.subsample(Genes_PB_pr_up, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

LZPB_gb_up_anti_meta <- meta.subsample(Genes_PB_gb_up, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZPB_gb_up_sense_meta <- meta.subsample(Genes_PB_gb_up, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

BF_pr_down_anti_meta <- meta.subsample(Genes_BF_pr_down, BF_12_minus, BF_12_plus,halfWindow,BinSize, at.TSS=T)
BF_pr_down_sense_meta <- meta.subsample(Genes_BF_pr_down, BF_12_plus, BF_12_minus,halfWindow,BinSize, at.TSS=T)

BF_gb_down_anti_meta <- meta.subsample(Genes_BF_gb_down, BF_12_minus, BF_12_plus,halfWindow,BinSize, at.TSS=T)
BF_gb_down_sense_meta <- meta.subsample(Genes_BF_gb_down, BF_12_plus, BF_12_minus,halfWindow,BinSize, at.TSS=T)

PB_pr_down_anti_meta <- meta.subsample(Genes_PB_pr_down, PB_12_minus, PB_12_plus,halfWindow,BinSize, at.TSS=T)
PB_pr_down_sense_meta <- meta.subsample(Genes_PB_pr_down, PB_12_plus, PB_12_minus,halfWindow,BinSize, at.TSS=T)

PB_gb_down_anti_meta <- meta.subsample(Genes_PB_gb_down, PB_12_minus, PB_12_plus,halfWindow,BinSize, at.TSS=T)
PB_gb_down_sense_meta <- meta.subsample(Genes_PB_gb_down, PB_12_plus, PB_12_minus,halfWindow,BinSize, at.TSS=T)

LZBF_pr_down_anti_meta <- meta.subsample(Genes_BF_pr_down, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZBF_pr_down_sense_meta <- meta.subsample(Genes_BF_pr_down, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

LZBF_gb_down_anti_meta <- meta.subsample(Genes_BF_gb_down, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZBF_gb_down_sense_meta <- meta.subsample(Genes_BF_gb_down, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

LZPB_pr_down_anti_meta <- meta.subsample(Genes_PB_pr_down, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZPB_pr_down_sense_meta <- meta.subsample(Genes_PB_pr_down, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

LZPB_gb_down_anti_meta <- meta.subsample(Genes_PB_gb_down, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZPB_gb_down_sense_meta <- meta.subsample(Genes_PB_gb_down, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

PBBF_pr_up_anti_meta <- meta.subsample(Genes_BF_pr_up, PB_12_minus, PB_12_plus,halfWindow,BinSize, at.TSS=T)
PBBF_pr_up_sense_meta <- meta.subsample(Genes_BF_pr_up, PB_12_plus, PB_12_minus,halfWindow,BinSize, at.TSS=T)

PBBF_gb_up_anti_meta <- meta.subsample(Genes_BF_gb_up, PB_12_minus, PB_12_plus,halfWindow,BinSize, at.TSS=T)
PBBF_gb_up_sense_meta <- meta.subsample(Genes_BF_gb_up, PB_12_plus, PB_12_minus,halfWindow,BinSize, at.TSS=T)

BFPB_pr_up_anti_meta <- meta.subsample(Genes_PB_pr_up, BF_12_minus, BF_12_plus,halfWindow,BinSize, at.TSS=T)
BFPB_pr_up_sense_meta <- meta.subsample(Genes_PB_pr_up, BF_12_plus, BF_12_minus,halfWindow,BinSize, at.TSS=T)

BFPB_gb_up_anti_meta <- meta.subsample(Genes_PB_gb_up, BF_12_minus, BF_12_plus,halfWindow,BinSize, at.TSS=T)
BFPB_gb_up_sense_meta <- meta.subsample(Genes_PB_gb_up, BF_12_plus, BF_12_minus,halfWindow,BinSize, at.TSS=T)

LZPBBF_pr_up_anti_meta <- meta.subsample(Genes_BF_pr_up, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZPBBF_pr_up_sense_meta <- meta.subsample(Genes_BF_pr_up, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

LZPBBF_gb_up_anti_meta <- meta.subsample(Genes_BF_gb_up, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZPBBF_gb_up_sense_meta <- meta.subsample(Genes_BF_gb_up, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

LZBFPB_pr_up_anti_meta <- meta.subsample(Genes_PB_pr_up, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZBFPB_pr_up_sense_meta <- meta.subsample(Genes_PB_pr_up, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

LZBFPB_gb_up_anti_meta <- meta.subsample(Genes_PB_gb_up, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZBFPB_gb_up_sense_meta <- meta.subsample(Genes_PB_gb_up, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

PBBF_pr_down_anti_meta <- meta.subsample(Genes_BF_pr_down, PB_12_minus, PB_12_plus,halfWindow,BinSize, at.TSS=T)
PBBF_pr_down_sense_meta <- meta.subsample(Genes_BF_pr_down, PB_12_plus, PB_12_minus,halfWindow,BinSize, at.TSS=T)

PBBF_gb_down_anti_meta <- meta.subsample(Genes_BF_gb_down, PB_12_minus, PB_12_plus,halfWindow,BinSize, at.TSS=T)
PBBF_gb_down_sense_meta <- meta.subsample(Genes_BF_gb_down, PB_12_plus, PB_12_minus,halfWindow,BinSize, at.TSS=T)

BFPB_pr_down_anti_meta <- meta.subsample(Genes_PB_pr_down, BF_12_minus, BF_12_plus,halfWindow,BinSize, at.TSS=T)
BFPB_pr_down_sense_meta <- meta.subsample(Genes_PB_pr_down, BF_12_plus, BF_12_minus,halfWindow,BinSize, at.TSS=T)

BFPB_gb_down_anti_meta <- meta.subsample(Genes_PB_gb_down, BF_12_minus, BF_12_plus,halfWindow,BinSize, at.TSS=T)
BFPB_gb_down_sense_meta <- meta.subsample(Genes_PB_gb_down, BF_12_plus, BF_12_minus,halfWindow,BinSize, at.TSS=T)

LZPBBF_pr_down_anti_meta <- meta.subsample(Genes_BF_pr_down, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZPBBF_pr_down_sense_meta <- meta.subsample(Genes_BF_pr_down, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

LZPBBF_gb_down_anti_meta <- meta.subsample(Genes_BF_gb_down, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZPBBF_gb_down_sense_meta <- meta.subsample(Genes_BF_gb_down, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

LZBFPB_pr_down_anti_meta <- meta.subsample(Genes_PB_pr_down, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZBFPB_pr_down_sense_meta <- meta.subsample(Genes_PB_pr_down, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

LZBFPB_gb_down_anti_meta <- meta.subsample(Genes_PB_gb_down, LZ_12_minus, LZ_12_plus,halfWindow,BinSize, at.TSS=T)
LZBFPB_gb_down_sense_meta <- meta.subsample(Genes_PB_gb_down, LZ_12_plus, LZ_12_minus,halfWindow,BinSize, at.TSS=T)

##Plot the PRO-seq counts around promoters: DESeq2 BF or PB RNAi up or down regulated genes relative to LacZ RNAi; mjg gene list
##Promoter PRO-seq counts: DESeq2 BF RNAi pr up regulated genes, BF RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZBF_pr_up_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 5000))
meta.overlay(BF_pr_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1)), y.intersp = 1.5)
text(-100, 1000,labels=bquote('n'==.(dim(Genes_BF_pr_up)[1])), cex=1)
title(main='BEAF RNAi: Promoter up', cex.main=0.75)

##Promoter PRO-seq counts: DESeq2 BF RNAi pr up regulated genes, PB RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZPBBF_pr_up_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 5000))
meta.overlay(PBBF_pr_up_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topleft", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 1000,labels=bquote('n'==.(dim(Genes_BF_pr_up)[1])), cex=1)
title(main='BEAF RNAi: Promoter up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_BFpr_up_PBrnai_p0_05.png')
#dev.off()

##Promoter PRO-seq counts: DESeq2 PB RNAi pr up regulated genes, PB RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZPB_pr_up_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 5000))
meta.overlay(PB_pr_up_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
#changed from col2=rgb(0, 1, 0, 0.25) to col2=rgb(0, 1, 0, 0) for all plots
legend("topleft", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 840,labels=bquote('n'==.(dim(Genes_PB_pr_up)[1])), cex=1)
title(main='Pbro RNAi: Promoter up', cex.main=0.75)

##Promoter PRO-seq counts: DESeq2 PB RNAi pr up regulated genes, BF RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZBFPB_pr_up_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 5000))
meta.overlay(BFPB_pr_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1)), y.intersp = 1.5)
text(-100, 840,labels=bquote('n'==.(dim(Genes_PB_pr_up)[1])), cex=1)
title(main='Pbro RNAi: Promoter up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_PBpr_up_BFrnai_p0_05.png')
#dev.off()

##Promoter PRO-seq counts: DESeq2 BF RNAi gb up regulated genes, BF RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZBF_gb_up_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 5000))
meta.overlay(BF_gb_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1)), y.intersp = 1.5)
text(-100, 840,labels=bquote('n'==.(dim(Genes_BF_gb_up)[1])), cex=1)
title(main='BEAF RNAi: Gene body up', cex.main=0.75)

##Promoter PRO-seq counts: DESeq2 BF RNAi gb up regulated genes, PB RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZPBBF_gb_up_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 5000))
meta.overlay(PBBF_gb_up_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topleft", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1, 0)), y.intersp = 1.5)
text(-100, 840,labels=bquote('n'==.(dim(Genes_BF_gb_up)[1])), cex=1)
title(main='BEAF RNAi: Gene body up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_BFgb_up_PBrnai_p0_05.png')
#dev.off()

##Promoter PRO-seq counts: DESeq2 PB RNAi gb up regulated genes, PB RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZPB_gb_up_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 10000))
meta.overlay(PB_gb_up_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topleft", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 1670,labels=bquote('n'==.(dim(Genes_PB_gb_up)[1])), cex=1)
title(main='Pbro RNAi: Gene body up', cex.main=0.75)

##Promoter PRO-seq counts: DESeq2 PB RNAi gb up regulated genes, BF RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZBFPB_gb_up_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 10000))
meta.overlay(BFPB_gb_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1)), y.intersp = 1.5)
text(-100, 1670,labels=bquote('n'==.(dim(Genes_PB_gb_up)[1])), cex=1)
title(main='Pbro RNAi: Gene body up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_PBgb_up_BFrnai_p0_05.png')
#dev.off()

##Promoter PRO-seq counts: DESeq2 BF RNAi pr dn regulated genes, BF RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZBF_pr_down_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 2500))
meta.overlay(BF_pr_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1)), y.intersp = 1.5)
text(-100, 420,labels=bquote('n'==.(dim(Genes_BF_pr_down)[1])), cex=1)
title(main='BEAF RNAi: Promoter down', cex.main=0.75)

##Promoter PRO-seq counts: DESeq2 BF RNAi pr dn regulated genes, PB RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZPBBF_pr_down_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 2500))
meta.overlay(PBBF_pr_down_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topleft", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1, 0)), y.intersp = 1.5)
text(-100, 420,labels=bquote('n'==.(dim(Genes_BF_pr_down)[1])), cex=1)
title(main='BEAF RNAi: Promoter down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_BFpr_dn_PBrnai_p0_05.png')
#dev.off()

##Promoter PRO-seq counts: DESeq2 PB RNAi pr dn regulated genes, PB RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZPB_pr_down_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 30000))
meta.overlay(PB_pr_down_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topleft", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 5000,labels=bquote('n'==.(dim(Genes_PB_pr_down)[1])), cex=1)
title(main='Pbro RNAi: Promoter down', cex.main=0.75)

##Promoter PRO-seq counts: DESeq2 PB RNAi pr dn regulated genes, BF RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZBFPB_pr_down_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 30000))
meta.overlay(BFPB_pr_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1)), y.intersp = 1.5)
text(-100, 5000,labels=bquote('n'==.(dim(Genes_PB_pr_down)[1])), cex=1)
title(main='Pbro RNAi: Promoter down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_PBpr_dn_BFrnai_p0_05.png')
#dev.off()

##Promoter PRO-seq counts: DESeq2 BF RNAi gb dn regulated genes, BF RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZBF_gb_down_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 2500))
meta.overlay(BF_gb_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1)), y.intersp = 1.5)
text(-100, 420,labels=bquote('n'==.(dim(Genes_BF_gb_down)[1])), cex=1)
title(main='BEAF RNAi: Gene body down', cex.main=0.75)

##Promoter PRO-seq counts: DESeq2 BF RNAi gb dn regulated genes, PB RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZPBBF_gb_down_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 2500))
meta.overlay(PBBF_gb_down_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topleft", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1, 0)), y.intersp = 1.5)
text(-100, 420,labels=bquote('n'==.(dim(Genes_BF_gb_down)[1])), cex=1)
title(main='BEAF RNAi: Gene body down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_BFgb_dn_PBrnai_p0_05.png')
#dev.off()

##Promoter PRO-seq counts: DESeq2 PB RNAi  dgbn regulated genes, PB RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZPB_gb_down_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 25000))
meta.overlay(PB_gb_down_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topleft", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 4200,labels=bquote('n'==.(dim(Genes_PB_gb_down)[1])), cex=1)
title(main='Pbro RNAi: Gene body down', cex.main=0.75)

##Promoter PRO-seq counts: DESeq2 PB RNAi  dgbn regulated genes, BF RNAi relative to LacZ RNAi; mjg gene list
meta.plot(LZBFPB_gb_down_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 25000))
meta.overlay(BFPB_gb_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1)), y.intersp = 1.5)
text(-100, 4200,labels=bquote('n'==.(dim(Genes_PB_gb_down)[1])), cex=1)
title(main='Pbro RNAi: Gene body down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_PBgb_dn_BFrnai_p0_05.png')
#dev.off()

##Make combined metaplots of both RNAi treatments with LacZ control
##Plot the PRO-seq counts around promoters: DESeq2 BF or PB RNAi up or down regulated genes relative to LacZ RNAi; mjg gene list
meta.plot(LZPBBF_pr_up_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 6000))
meta.overlay(PBBF_pr_up_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
meta.overlay(BF_pr_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 1000,labels=bquote('n'==.(dim(Genes_BF_pr_up)[1])), cex=1)
title(main='BEAF RNAi: Promoter up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_BFpr_up_BFPBrnai_p0_05.png')
#dev.off()

meta.plot(LZBFPB_pr_up_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 5000))
meta.overlay(PB_pr_up_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
meta.overlay(BFPB_pr_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 840,labels=bquote('n'==.(dim(Genes_PB_pr_up)[1])), cex=1)
title(main='Pbro RNAi: Promoter up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_PBpr_up_BFPBrnai_p0_05.png')
#dev.off()

meta.plot(LZPBBF_gb_up_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 6000))
meta.overlay(PBBF_gb_up_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
meta.overlay(BF_gb_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1, 0)), y.intersp = 1.5)
text(-100, 840,labels=bquote('n'==.(dim(Genes_BF_gb_up)[1])), cex=1)
title(main='BEAF RNAi: Gene body up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_BFgb_up_BFPBrnai_p0_05.png')
#dev.off()

meta.plot(LZBFPB_gb_up_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 10000))
meta.overlay(PB_gb_up_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
meta.overlay(BFPB_gb_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 1670,labels=bquote('n'==.(dim(Genes_PB_gb_up)[1])), cex=1)
title(main='Pbro RNAi: Gene body up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_PBgb_up_BFPBrnai_p0_05.png')
#dev.off()

##Plot the PRO-seq counts around promoters: DESeq2 BF or PB RNAi up or down regulated genes relative to LacZ RNAi; BF_1818 genes
meta.plot(LZPBBF_pr_down_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 2000))
meta.overlay(PBBF_pr_down_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
meta.overlay(BF_pr_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1, 0)), y.intersp = 1.5)
text(-100, 420,labels=bquote('n'==.(dim(Genes_BF_pr_down)[1])), cex=1)
title(main='BEAF RNAi: Promoter down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_BFpr_dn_BFPBrnai_p0_05.png')
#dev.off()

meta.plot(LZBFPB_pr_down_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 30000))
meta.overlay(PB_pr_down_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
meta.overlay(BFPB_pr_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 5000,labels=bquote('n'==.(dim(Genes_PB_pr_down)[1])), cex=1)
title(main='Pbro RNAi: Promoter down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_PBpr_dn_BFPBrnai_p0_05.png')
#dev.off()

meta.plot(LZPBBF_gb_down_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 2500))
meta.overlay(PBBF_gb_down_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
meta.overlay(BF_gb_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1, 0)), y.intersp = 1.5)
text(-100, 420,labels=bquote('n'==.(dim(Genes_BF_gb_down)[1])), cex=1)
title(main='BEAF RNAi: Gene body down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_BFgb_dn_BFPBrnai_p0_05.png')
#dev.off()

meta.plot(LZBFPB_gb_down_sense_meta, BinSize, xlab="Distance from TSS (bp)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 25000))
meta.overlay(PB_gb_down_sense_meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
meta.overlay(BFPB_gb_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0))
legend("topleft", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 4200,labels=bquote('n'==.(dim(Genes_PB_gb_down)[1])), cex=1)
title(main='Pbro RNAi: Gene body down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/DESeq_PBgb_dn_BFPBrnai_p0_05.png')
#dev.off()


#####The end
