#Modify /media/jacobmtome/HDD/JMT_Lis_Lab/2015_07_PRO_seq/Quantify_Shift/KernelDensity/KernelDensity_AllGenes/WT_vs_Apt/SeminarFigure/151018_CompareWTvsApt_hMapForSeminar.R
  #to use the new data, afterm more sequencing


#Base off of /media/jacobmtome/HDD/JMT_Lis_Lab/2015_07_PRO_seq/Quantify_Shift/KernelDensity/KernelDensity_AllGenes/Test_Parameters/150930_UseKernelDensity_AllGenes_TestParameters.R

#After picking a minimum read count and bandwidth, do the comparison of WT cells to aptamer... both copper treated
  #I saw what I was hoping to see with just the copper treated aptamer containing cells... I am thinking I will see the same thing but amplified by looking this way
  

library(bigWig)
library(gplots)

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


LZ_12_plus = load.bigWig('/home/chart/PROseq2/bigwig/LacZ_12_plus_noMito.bw')
LZ_12_minus = load.bigWig('/home/chart/PROseq2/bigwig/LacZ_12_minus_noMito.bw')

BF_12_plus = load.bigWig('/home/chart/PROseq2/bigwig/BEAF_12_plus_noMito.bw')
BF_12_minus = load.bigWig('/home/chart/PROseq2/bigwig/BEAF_12_minus_noMito.bw')

PB_12_plus = load.bigWig('/home/chart/PROseq2/bigwig/PBro_12_plus_noMito.bw')
PB_12_minus = load.bigWig('/home/chart/PROseq2/bigwig/PBro_12_minus_noMito.bw')


##################################


collect.scaled_gbody = function (bed, bigWig.plus, bigWig.minus, NumberBins, do.sum = TRUE) {
  #a collect.many like function to collect counts around the promoter
  #goes from -50 to plus 250 of the TSS
  N = dim(bed)[1]
  result = matrix(nrow = N, ncol = NumberBins)
  for (i in 1:N) {
    step = (bed[i, 3] - bed[i, 2] - 400)%/%NumberBins
    chrom = as.character(bed[i, 1])
    strand = as.character(bed[i, 6])
    if (strand == "+") {
      bigWig = bigWig.plus
      start = bed[i, 2]+200    #start 200 bp ds of TSS
#      start = bed[i, 2]+400    #start 400 bp ds of TSS
      end = start + step * NumberBins
      row = collect.counts(bigWig, chrom, start, end, step, do.sum)
      result[i, ] = abs(row)/step  #divide counts by step size to normalize... so reads per base
    }
    else
    { #genes on the minus strand will be different.... remainder part from bin calculation is taken off of 5' end, rather than 3' as for plus strand genes
      bigWig = bigWig.minus
      start = bed[i, 2]+200    #start 200 bp ds of TSS
      end = start + step * NumberBins
      row = collect.counts(bigWig, chrom, start, end, step, do.sum)
      result[i, ] = abs(rev(row))/step    #for the minus strand, reverse
    }
  }
  result
}



meta.subsample <- function (bed, bigWig.plus, bigWig.minus, NumberBins, at.TSS = TRUE, do.sum = TRUE) {
  N = dim(bed)[1]
  nPermut = 1000
  sampleFrac = 0.2
  result = matrix(nrow = nPermut, ncol = NumberBins)
  M = as.integer(round(N * sampleFrac, 0))
  
  values = collect.scaled_gbody(bed, bigWig.plus, bigWig.minus, NumberBins)
  for (i in 1:nPermut) {
    idx <- sample(N, size = M, replace = T)
    result[i, ] = colSums(values[idx, ])
  }
#changed ci9 from 0.875 to 0.5 and ci1 from 0.125 to 0.5 to hide ci
  ci9 = sapply(1:NumberBins, function(idx) quantile(result[, idx], 0.5))
  ci1 = sapply(1:NumberBins, function(idx) quantile(result[, idx], 0.5))
  ci5 = sapply(1:NumberBins, function(idx) quantile(result[, idx], 0.5))
  return(list(result, ci9, ci1, ci5, values))
}

##Try adding meta.overlay to plot experimental with control
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

##experimental: is meta.overlay working?
BinSize = 1
#BEAF pr up gene list, BEAF and LacZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_BF_pr_up, BF_12_plus, BF_12_minus, 300)
BF_pr_up_Meta = meta.subsample(Genes_BF_pr_up, BF_12_plus, BF_12_minus, 300)
LZBF_pr_up_Meta = meta.subsample(Genes_BF_pr_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_pr_up_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 20))
#meta.plot(LZBF_pr_up_Meta, 1, ylim=c(0, 20))
meta.overlay(BF_pr_up_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topright", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1)), y.intersp = 1.5)
text(-100, 16.7,labels=bquote('n'==.(dim(Genes_BF_pr_up)[1])), cex=1)
title(main='BEAF RNAi: Promoter up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_BFpr_up_BFrnai_p0_05.png')
#dev.off()

#BEAF pr up gene list, Pbro and LacZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_BF_pr_up, BF_12_plus, BF_12_minus, 300)
PBBF_pr_up_Meta = meta.subsample(Genes_BF_pr_up, PB_12_plus, PB_12_minus, 300)
LZBF_pr_up_Meta = meta.subsample(Genes_BF_pr_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_pr_up_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 20))
#meta.plot(LZBF_pr_up_Meta, 1, ylim=c(0, 20))
meta.overlay(PBBF_pr_up_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topright", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1, 0)), y.intersp = 1.5)
text(-100, 16.7,labels=bquote('n'==.(dim(Genes_BF_pr_up)[1])), cex=1)
title(main='BEAF RNAi: Promoter up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_BFpr_up_PBrnai_p0_05.png')
#dev.off()

#BEAF pr up gene list, BEAF and Pbro and LacZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_BF_pr_up, BF_12_plus, BF_12_minus, 300)
BF_pr_up_Meta = meta.subsample(Genes_BF_pr_up, BF_12_plus, BF_12_minus, 300)
PBBF_pr_up_Meta = meta.subsample(Genes_BF_pr_up, PB_12_plus, PB_12_minus, 300)
LZBF_pr_up_Meta = meta.subsample(Genes_BF_pr_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_pr_up_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 20))
#meta.plot(LZBF_pr_up_Meta, 1, ylim=c(0, 20))
meta.overlay(BF_pr_up_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
meta.overlay(PBBF_pr_up_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topright", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1, 0)), y.intersp = 1.5)
text(-100, 16.7,labels=bquote('n'==.(dim(Genes_BF_pr_up)[1])), cex=1)
title(main='BEAF RNAi: Promoter up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_BFpr_up_BFPBrnai_p0_05.png')
#dev.off()

#BEAF gb up gene list, BEAF and LacZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_BF_gb_up, BF_12_plus, BF_12_minus, 300)
BF_gb_up_Meta = meta.subsample(Genes_BF_gb_up, BF_12_plus, BF_12_minus, 300)
LZBF_gb_up_Meta = meta.subsample(Genes_BF_gb_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_gb_up_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 10))
meta.overlay(BF_gb_up_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 8.3,labels=bquote('n'==.(dim(Genes_BF_gb_up)[1])), cex=1)
title(main='BEAF RNAi: Gene body up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_BFgb_up_BFrnai_p0_05.png')
#dev.off()

#BEAF gb up gene list, Pbro and LacZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_BF_gb_up, BF_12_plus, BF_12_minus, 300)
PBBF_gb_up_Meta = meta.subsample(Genes_BF_gb_up, PB_12_plus, PB_12_minus, 300)
LZBF_gb_up_Meta = meta.subsample(Genes_BF_gb_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_gb_up_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 10))
meta.overlay(PBBF_gb_up_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0.25))
legend("topright", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1, 0)), y.intersp = 1.5)
text(-100, 8.3,labels=bquote('n'==.(dim(Genes_BF_gb_up)[1])), cex=1)
title(main='BEAF RNAi: Gene body up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_BFgb_up_PBrnai_p0_05.png')
#dev.off()

#BEAF gb up gene list, BEAF and Pbro and LacZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_BF_gb_up, BF_12_plus, BF_12_minus, 300)
BF_gb_up_Meta = meta.subsample(Genes_BF_gb_up, BF_12_plus, BF_12_minus, 300)
PBBF_gb_up_Meta = meta.subsample(Genes_BF_gb_up, PB_12_plus, PB_12_minus, 300)
LZBF_gb_up_Meta = meta.subsample(Genes_BF_gb_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_gb_up_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 10))
meta.overlay(BF_gb_up_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0.25))
meta.overlay(PBBF_gb_up_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0.25))
legend("topright", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1, 0)), y.intersp = 1.5)
text(-100, 8.3,labels=bquote('n'==.(dim(Genes_BF_gb_up)[1])), cex=1)
title(main='BEAF RNAi: Gene body up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_BFgb_up_BFPBrnai_p0_05.png')
#dev.off()

#BEAF pr down gene list, BEAF and LacZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_BF_pr_down, BF_12_plus, BF_12_minus, 300)
BF_pr_down_Meta = meta.subsample(Genes_BF_pr_down, BF_12_plus, BF_12_minus, 300)
LZBF_pr_down_Meta = meta.subsample(Genes_BF_pr_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_pr_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 15))
meta.overlay(BF_pr_down_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topright", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 12.5,labels=bquote('n'==.(dim(Genes_BF_pr_down)[1])), cex=1)
title(main='BEAF RNAi: Promoter down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_BFpr_dn_BFrnai_p0_05.png')
#dev.off()

#BEAF pr down gene list, Pbro and LacZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_BF_pr_down, BF_12_plus, BF_12_minus, 300)
PBBF_pr_down_Meta = meta.subsample(Genes_BF_pr_down, PB_12_plus, PB_12_minus, 300)
LZBF_pr_down_Meta = meta.subsample(Genes_BF_pr_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_pr_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 15))
meta.overlay(PBBF_pr_down_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topright", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1, 0)), y.intersp = 1.5)
text(-100, 12.5,labels=bquote('n'==.(dim(Genes_BF_pr_down)[1])), cex=1)
title(main='BEAF RNAi: Promoter down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_BFpr_dn_PBrnai_p0_05.png')
#dev.off()

#BEAF pr down gene list, BEAF and Pbro and LacZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_BF_pr_down, BF_12_plus, BF_12_minus, 300)
BF_pr_down_Meta = meta.subsample(Genes_BF_pr_down, BF_12_plus, BF_12_minus, 300)
PBBF_pr_down_Meta = meta.subsample(Genes_BF_pr_down, PB_12_plus, PB_12_minus, 300)
LZBF_pr_down_Meta = meta.subsample(Genes_BF_pr_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_pr_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 15))
meta.overlay(BF_pr_down_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
meta.overlay(PBBF_pr_down_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topright", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1, 0)), y.intersp = 1.5)
text(-100, 12.5,labels=bquote('n'==.(dim(Genes_BF_pr_down)[1])), cex=1)
title(main='BEAF RNAi: Promoter down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_BFpr_dn_BFPBrnai_p0_05.png')
#dev.off()

#BEAF gb down gene list, BEAF and LacZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_BF_gb_down, BF_12_plus, BF_12_minus, 300)
BF_gb_down_Meta = meta.subsample(Genes_BF_gb_down, BF_12_plus, BF_12_minus, 300)
LZBF_gb_down_Meta = meta.subsample(Genes_BF_gb_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_gb_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 45))
#meta.plot(LZBF_gb_down_Meta, 1, ylim=c(0, 30))
meta.overlay(BF_gb_down_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topright", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 37.5,labels=bquote('n'==.(dim(Genes_BF_gb_down)[1])), cex=1)
title(main='BEAF RNAi: Gene body down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_BFgb_dn_BFrnai_p0_05.png')
#dev.off()

#BEAF gb down gene list, Pbro and LacZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_BF_gb_down, BF_12_plus, BF_12_minus, 300)
PBBF_gb_down_Meta = meta.subsample(Genes_BF_gb_down, PB_12_plus, PB_12_minus, 300)
LZBF_gb_down_Meta = meta.subsample(Genes_BF_gb_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_gb_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 45))
meta.overlay(PBBF_gb_down_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topright", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1, 0)), y.intersp = 1.5)
text(-100, 37.5,labels=bquote('n'==.(dim(Genes_BF_gb_down)[1])), cex=1)
title(main='BEAF RNAi: Gene body down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_BFgb_dn_PBrnai_p0_05.png')
#dev.off()

#BEAF gb down gene list, BEAF and Pbro and LacZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_BF_gb_down, BF_12_plus, BF_12_minus, 300)
BF_gb_down_Meta = meta.subsample(Genes_BF_gb_down, BF_12_plus, BF_12_minus, 300)
PBBF_gb_down_Meta = meta.subsample(Genes_BF_gb_down, PB_12_plus, PB_12_minus, 300)
LZBF_gb_down_Meta = meta.subsample(Genes_BF_gb_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_gb_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 45))
meta.overlay(BF_gb_down_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
meta.overlay(PBBF_gb_down_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topright", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1, 0)), y.intersp = 1.5)
text(-100, 37.5,labels=bquote('n'==.(dim(Genes_BF_gb_down)[1])), cex=1)
title(main='BEAF RNAi: Gene body down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_BFgb_dn_BFPBrnai_p0_05.png')
#dev.off()

#PBro pr up gene list, Pbro and LacZ
gb_scaled_matrix = collect.scaled_gbody(Genes_PB_pr_up, PB_12_plus, PB_12_minus, 300)
PB_pr_up_Meta = meta.subsample(Genes_PB_pr_up, PB_12_plus, PB_12_minus, 300)
LZPB_pr_up_Meta = meta.subsample(Genes_PB_pr_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_pr_up_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 30))
meta.overlay(PB_pr_up_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topright", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 25,labels=bquote('n'==.(dim(Genes_PB_pr_up)[1])), cex=1)
title(main='Pbro RNAi: Promoter up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_PBpr_up_PBrnai_p0_05.png')
#dev.off()

#PBro pr up gene list, BEAF and LacZ
gb_scaled_matrix = collect.scaled_gbody(Genes_PB_pr_up, PB_12_plus, PB_12_minus, 300)
BFPB_pr_up_Meta = meta.subsample(Genes_PB_pr_up, BF_12_plus, BF_12_minus, 300)
LZPB_pr_up_Meta = meta.subsample(Genes_PB_pr_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_pr_up_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 30))
meta.overlay(BFPB_pr_up_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topright", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1)), y.intersp = 1.5)
text(-100, 25,labels=bquote('n'==.(dim(Genes_PB_pr_up)[1])), cex=1)
title(main='Pbro RNAi: Promoter up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_PBpr_up_BFrnai_p0_05.png')
#dev.off()

#PBro pr up gene list, BEAF and Pbro and LacZ
gb_scaled_matrix = collect.scaled_gbody(Genes_PB_pr_up, PB_12_plus, PB_12_minus, 300)
PB_pr_up_Meta = meta.subsample(Genes_PB_pr_up, PB_12_plus, PB_12_minus, 300)
BFPB_pr_up_Meta = meta.subsample(Genes_PB_pr_up, BF_12_plus, BF_12_minus, 300)
LZPB_pr_up_Meta = meta.subsample(Genes_PB_pr_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_pr_up_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 30))
meta.overlay(PB_pr_up_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
meta.overlay(BFPB_pr_up_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topright", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 25,labels=bquote('n'==.(dim(Genes_PB_pr_up)[1])), cex=1)
title(main='Pbro RNAi: Promoter up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_PBpr_up_BFPBrnai_p0_05.png')
#dev.off()

#PBro gb up gene list, Pbro and LacZ
gb_scaled_matrix = collect.scaled_gbody(Genes_PB_gb_up, PB_12_plus, PB_12_minus, 300)
PB_gb_up_Meta = meta.subsample(Genes_PB_gb_up, PB_12_plus, PB_12_minus, 300)
LZPB_gb_up_Meta = meta.subsample(Genes_PB_gb_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_gb_up_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 30))
meta.overlay(PB_gb_up_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0.25))
legend("topright", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 25,labels=bquote('n'==.(dim(Genes_PB_gb_up)[1])), cex=1)
title(main='Pbro RNAi: Gene body up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_PBgb_up_PBrnai_p0_05.png')
#dev.off()

#PBro gb up gene list, BEAF and LacZ
gb_scaled_matrix = collect.scaled_gbody(Genes_PB_gb_up, PB_12_plus, PB_12_minus, 300)
BFPB_gb_up_Meta = meta.subsample(Genes_PB_gb_up, BF_12_plus, BF_12_minus, 300)
LZPB_gb_up_Meta = meta.subsample(Genes_PB_gb_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_gb_up_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 30))
meta.overlay(BFPB_gb_up_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topright", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1)), y.intersp = 1.5)
text(-100, 25,labels=bquote('n'==.(dim(Genes_PB_gb_up)[1])), cex=1)
title(main='Pbro RNAi: Gene body up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_PBgb_up_BFrnai_p0_05.png')
#dev.off()

#PBro gb up gene list, BEAF and Pbro and LacZ
gb_scaled_matrix = collect.scaled_gbody(Genes_PB_gb_up, PB_12_plus, PB_12_minus, 300)
PB_gb_up_Meta = meta.subsample(Genes_PB_gb_up, PB_12_plus, PB_12_minus, 300)
BFPB_gb_up_Meta = meta.subsample(Genes_PB_gb_up, BF_12_plus, BF_12_minus, 300)
LZPB_gb_up_Meta = meta.subsample(Genes_PB_gb_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_gb_up_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 30))
meta.overlay(PB_gb_up_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
meta.overlay(BFPB_gb_up_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topright", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 25,labels=bquote('n'==.(dim(Genes_PB_gb_up)[1])), cex=1)
title(main='Pbro RNAi: Gene body up', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_PBgb_up_BFPBrnai_p0_05.png')
#dev.off()

#PBro pr down gene list, Pbro and LacZ
gb_scaled_matrix = collect.scaled_gbody(Genes_PB_pr_down, PB_12_plus, PB_12_minus, 300)
PB_pr_down_Meta = meta.subsample(Genes_PB_pr_down, PB_12_plus, PB_12_minus, 300)
LZPB_pr_down_Meta = meta.subsample(Genes_PB_pr_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_pr_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 25))
meta.overlay(PB_pr_down_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topright", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 20.8,labels=bquote('n'==.(dim(Genes_PB_pr_down)[1])), cex=1)
title(main='Pbro RNAi: Promoter down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_PBpr_dn_PBrnai_p0_05.png')
#dev.off()

#PBro pr down gene list, BEAF and LacZ
gb_scaled_matrix = collect.scaled_gbody(Genes_PB_pr_down, PB_12_plus, PB_12_minus, 300)
BFPB_pr_down_Meta = meta.subsample(Genes_PB_pr_down, BF_12_plus, BF_12_minus, 300)
LZPB_pr_down_Meta = meta.subsample(Genes_PB_pr_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_pr_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 25))
meta.overlay(BFPB_pr_down_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topright", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1)), y.intersp = 1.5)
text(-100, 20.8,labels=bquote('n'==.(dim(Genes_PB_pr_down)[1])), cex=1)
title(main='Pbro RNAi: Promoter down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_PBpr_dn_BFrnai_p0_05.png')
#dev.off()

#PBro pr down gene list, BEAF and Pbro and LacZ
gb_scaled_matrix = collect.scaled_gbody(Genes_PB_pr_down, PB_12_plus, PB_12_minus, 300)
BFPB_pr_down_Meta = meta.subsample(Genes_PB_pr_down, BF_12_plus, BF_12_minus, 300)
PB_pr_down_Meta = meta.subsample(Genes_PB_pr_down, PB_12_plus, PB_12_minus, 300)
LZPB_pr_down_Meta = meta.subsample(Genes_PB_pr_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_pr_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 25))
meta.overlay(PB_pr_down_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
meta.overlay(BFPB_pr_down_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topright", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1,0)), y.intersp = 1.5)
text(-100, 20.8,labels=bquote('n'==.(dim(Genes_PB_pr_down)[1])), cex=1)
title(main='Pbro RNAi: Promoter down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_PBpr_dn_BFPBrnai_p0_05.png')
#dev.off()

#PBro gb down gene list, Pbro and LacZ
gb_scaled_matrix = collect.scaled_gbody(Genes_PB_gb_down, PB_12_plus, PB_12_minus, 300)
PB_gb_down_Meta = meta.subsample(Genes_PB_gb_down, PB_12_plus, PB_12_minus, 300)
LZPB_gb_down_Meta = meta.subsample(Genes_PB_gb_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_gb_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 350))
#meta.plot(LZPB_gb_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 120))
#meta.plot(LZPB_gb_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 50))
meta.overlay(PB_gb_down_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
legend("topright", c("LacZ RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,1,0)), y.intersp = 1.5)
text(-50, 292,labels=bquote('n'==.(dim(Genes_PB_gb_down)[1])), cex=1)
#text(-50, 100,labels=bquote('n'==.(dim(Genes_PB_gb_down)[1])), cex=1)
#text(-50, 42,labels=bquote('n'==.(dim(Genes_PB_gb_down)[1])), cex=1)
title(main='Pbro RNAi: Gene body down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_PBgb_dn_PBrnai_p0_05.png')
#dev.off()

#PBro gb down gene list, BEAF and LacZ
gb_scaled_matrix = collect.scaled_gbody(Genes_PB_gb_down, PB_12_plus, PB_12_minus, 300)
BFPB_gb_down_Meta = meta.subsample(Genes_PB_gb_down, BF_12_plus, BF_12_minus, 300)
LZPB_gb_down_Meta = meta.subsample(Genes_PB_gb_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_gb_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 350))
#meta.plot(LZPB_gb_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 120))
#meta.plot(LZPB_gb_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 50))
meta.overlay(BFPB_gb_down_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topright", c("LacZ RNAi", "BEAF RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1)), y.intersp = 1.5)
text(-50, 292,labels=bquote('n'==.(dim(Genes_PB_gb_down)[1])), cex=1)
#text(-50, 100,labels=bquote('n'==.(dim(Genes_PB_gb_down)[1])), cex=1)
#text(-50, 42,labels=bquote('n'==.(dim(Genes_PB_gb_down)[1])), cex=1)
title(main='Pbro RNAi: Gene body down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_PBgb_dn_BFrnai_p0_05.png')
#dev.off()

#PBro gb down gene list, BEAF and Pbro and LacZ
gb_scaled_matrix = collect.scaled_gbody(Genes_PB_gb_down, PB_12_plus, PB_12_minus, 300)
PB_gb_down_Meta = meta.subsample(Genes_PB_gb_down, PB_12_plus, PB_12_minus, 300)
BFPB_gb_down_Meta = meta.subsample(Genes_PB_gb_down, BF_12_plus, BF_12_minus, 300)
LZPB_gb_down_Meta = meta.subsample(Genes_PB_gb_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_gb_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 350))
#meta.plot(LZPB_gb_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 120))
#meta.plot(LZPB_gb_down_Meta, 1, xlab="Distance from gene center (scaled)", cex.lab=1, cex.axis=1, ylab="PRO-seq Signal", main=NULL, ylim=c(0, 50))
meta.overlay(PB_gb_down_Meta, BinSize, strand='+', col1=rgb(0, 1, 0, 0.75), col2=rgb(0, 1, 0, 0))
meta.overlay(BFPB_gb_down_Meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.75), col2=rgb(0, 0, 1, 0))
legend("topright", c("LacZ RNAi", "BEAF RNAi", "Pbro RNAi"), cex=0.75, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0,1), rgb(0,1,0)), y.intersp = 1.5)
text(-50, 292,labels=bquote('n'==.(dim(Genes_PB_gb_down)[1])), cex=1)
#text(-50, 100,labels=bquote('n'==.(dim(Genes_PB_gb_down)[1])), cex=1)
#text(-50, 42,labels=bquote('n'==.(dim(Genes_PB_gb_down)[1])), cex=1)
title(main='Pbro RNAi: Gene body down', cex.main=0.75)

#dev.copy(png, '/home/chart/metaplots/DESeq_genes/2022may_switchedBfPb/switchedBfPb_gb/gbDESeq_PBgb_dn_BFPBrnai_p0_05.png')
#dev.off()


#not modified below this point (april 2016) because I do not need single plots
#This works for single plots:BF up genes, BF and LZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_BF_pr_up, BF_12_plus, BF_12_minus, 300)
BF_pr_up_Meta = meta.subsample(Genes_BF_pr_up, BF_12_plus, BF_12_minus, 300)
meta.plot(BF_pr_up_Meta, 1)
#meta.overlay(BF_pr_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("BEAF RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 30,labels=bquote('n'==.(dim(Genes_BF_pr_up)[1])), cex=1.25)
title(main='Scaled gene bodies of promoter up \nregulated genes after BEAF RNAi', cex=1)

gb_scaled_matrix = collect.scaled_gbody(Genes_BF_gb_up, BF_12_plus, BF_12_minus, 300)
BF_gb_up_Meta = meta.subsample(Genes_BF_gb_up, BF_12_plus, BF_12_minus, 300)
meta.plot(BF_gb_up_Meta, 1)
#meta.overlay(BF_gb_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("BEAF RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 15,labels=bquote('n'==.(dim(Genes_BF_gb_up)[1])), cex=1.25)
title(main='Scaled gene bodies of body up \nregulated genes after BEAF RNAi', cex=1)

gb_scaled_matrix = collect.scaled_gbody(Genes_BF_pr_up, LZ_12_plus, LZ_12_minus, 300)
LZBF_pr_up_Meta = meta.subsample(Genes_BF_pr_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_pr_up_Meta, 1)
#meta.overlay(BF_pr_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("LacZ RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 15,labels=bquote('n'==.(dim(Genes_BF_pr_up)[1])), cex=1.25)
title(main='LacZ scaled gene bodies of promoter up \nregulated genes after BEAF RNAi', cex=1)

gb_scaled_matrix = collect.scaled_gbody(Genes_BF_gb_up, LZ_12_plus, LZ_12_minus, 300)
LZBF_gb_up_Meta = meta.subsample(Genes_BF_gb_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_gb_up_Meta, 1)
#meta.overlay(BF_gb_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("LacZ RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 8,labels=bquote('n'==.(dim(Genes_BF_gb_up)[1])), cex=1.25)
title(main='LacZ scaled gene bodies of body up \nregulated genes after BEAF RNAi', cex=1)


#This works for single plots:PB up genes, PB and LZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_PB_pr_up, PB_12_plus, PB_12_minus, 300)
PB_pr_up_Meta = meta.subsample(Genes_PB_pr_up, PB_12_plus, PB_12_minus, 300)
meta.plot(PB_pr_up_Meta, 1)
#meta.overlay(PB_pr_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("Pbro RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 30,labels=bquote('n'==.(dim(Genes_PB_pr_up)[1])), cex=1.25)
title(main='Scaled gene bodies of promoter up \nregulated genes after Pbro RNAi', cex=1)

gb_scaled_matrix = collect.scaled_gbody(Genes_PB_gb_up, PB_12_plus, PB_12_minus, 300)
PB_gb_up_Meta = meta.subsample(Genes_PB_gb_up, PB_12_plus, PB_12_minus, 300)
meta.plot(PB_gb_up_Meta, 1)
#meta.overlay(PB_gb_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("Pbro RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 500,labels=bquote('n'==.(dim(Genes_PB_gb_up)[1])), cex=1.25)
title(main='Scaled gene bodies of body up \nregulated genes after Pbro RNAi', cex=1)

gb_scaled_matrix = collect.scaled_gbody(Genes_PB_pr_up, LZ_12_plus, LZ_12_minus, 300)
LZPB_pr_up_Meta = meta.subsample(Genes_PB_pr_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_pr_up_Meta, 1)
#meta.overlay(PB_pr_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("LacZ RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 50,labels=bquote('n'==.(dim(Genes_PB_pr_up)[1])), cex=1.25)
title(main='LacZ scaled gene bodies of promoter up \nregulated genes after Pbro RNAi', cex=1)

gb_scaled_matrix = collect.scaled_gbody(Genes_PB_gb_up, LZ_12_plus, LZ_12_minus, 300)
LZPB_gb_up_Meta = meta.subsample(Genes_PB_gb_up, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_gb_up_Meta, 1)
#meta.overlay(PB_gb_up_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("LacZ RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 20,labels=bquote('n'==.(dim(Genes_PB_gb_up)[1])), cex=1.25)
title(main='LacZ scaled gene bodies of body up \nregulated genes after Pbro RNAi', cex=1)


#This works for single plots:BF down genes, BF and LZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_BF_pr_down, BF_12_plus, BF_12_minus, 300)
BF_pr_down_Meta = meta.subsample(Genes_BF_pr_down, BF_12_plus, BF_12_minus, 300)
meta.plot(BF_pr_down_Meta, 1)
#meta.overlay(BF_pr_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("BEAF RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 30,labels=bquote('n'==.(dim(Genes_BF_pr_down)[1])), cex=1.25)
title(main='Scaled gene bodies of promoter down \nregulated genes after BEAF RNAi', cex=1)

gb_scaled_matrix = collect.scaled_gbody(Genes_BF_gb_down, BF_12_plus, BF_12_minus, 300)
BF_gb_down_Meta = meta.subsample(Genes_BF_gb_down, BF_12_plus, BF_12_minus, 300)
meta.plot(BF_gb_down_Meta, 1)
#meta.overlay(BF_gb_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("BEAF RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 15,labels=bquote('n'==.(dim(Genes_BF_gb_down)[1])), cex=1.25)
title(main='Scaled gene bodies of body down \nregulated genes after BEAF RNAi', cex=1)

gb_scaled_matrix = collect.scaled_gbody(Genes_BF_pr_down, LZ_12_plus, LZ_12_minus, 300)
LZBF_pr_down_Meta = meta.subsample(Genes_BF_pr_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_pr_down_Meta, 1)
#meta.overlay(BF_pr_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("LacZ RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 15,labels=bquote('n'==.(dim(Genes_BF_pr_down)[1])), cex=1.25)
title(main='LacZ scaled gene bodies of promoter down \nregulated genes after BEAF RNAi', cex=1)

gb_scaled_matrix = collect.scaled_gbody(Genes_BF_gb_down, LZ_12_plus, LZ_12_minus, 300)
LZBF_gb_down_Meta = meta.subsample(Genes_BF_gb_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZBF_gb_down_Meta, 1)
#meta.overlay(BF_gb_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("LacZ RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 8,labels=bquote('n'==.(dim(Genes_BF_gb_down)[1])), cex=1.25)
title(main='LacZ scaled gene bodies of body down \nregulated genes after BEAF RNAi', cex=1)


#This works for single plots:PB down genes, PB and LZ RNAi
gb_scaled_matrix = collect.scaled_gbody(Genes_PB_pr_down, PB_12_plus, PB_12_minus, 300)
PB_pr_down_Meta = meta.subsample(Genes_PB_pr_down, PB_12_plus, PB_12_minus, 300)
meta.plot(PB_pr_down_Meta, 1)
#meta.overlay(PB_pr_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("Pbro RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 30,labels=bquote('n'==.(dim(Genes_PB_pr_down)[1])), cex=1.25)
title(main='Scaled gene bodies of promoter down \nregulated genes after Pbro RNAi', cex=1)

gb_scaled_matrix = collect.scaled_gbody(Genes_PB_gb_down, PB_12_plus, PB_12_minus, 300)
PB_gb_down_Meta = meta.subsample(Genes_PB_gb_down, PB_12_plus, PB_12_minus, 300)
meta.plot(PB_gb_down_Meta, 1)
#meta.overlay(PB_gb_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("Pbro RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 15,labels=bquote('n'==.(dim(Genes_PB_gb_down)[1])), cex=1.25)
title(main='Scaled gene bodies of body down \nregulated genes after Pbro RNAi', cex=1)

gb_scaled_matrix = collect.scaled_gbody(Genes_PB_pr_down, LZ_12_plus, LZ_12_minus, 300)
LZPB_pr_down_Meta = meta.subsample(Genes_PB_pr_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_pr_down_Meta, 1)
#meta.overlay(PB_pr_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("LacZ RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 15,labels=bquote('n'==.(dim(Genes_PB_pr_down)[1])), cex=1.25)
title(main='LacZ scaled gene bodies of promoter down \nregulated genes after Pbro RNAi', cex=1)

gb_scaled_matrix = collect.scaled_gbody(Genes_PB_gb_down, LZ_12_plus, LZ_12_minus, 300)
LZPB_gb_down_Meta = meta.subsample(Genes_PB_gb_down, LZ_12_plus, LZ_12_minus, 300)
meta.plot(LZPB_gb_down_Meta, 1)
#meta.overlay(PB_gb_down_sense_meta, BinSize, strand='+', col1=rgb(0, 0, 1, 0.5), col2=rgb(0, 0, 1, 0.25))
legend("topright", c("LacZ RNAi"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(rgb(1,0,0), rgb(0,0, 1)), y.intersp = 1.5)
text(-100, 8,labels=bquote('n'==.(dim(Genes_PB_gb_down)[1])), cex=1.25)
title(main='LacZ scaled gene bodies of body down \nregulated genes after Pbro RNAi', cex=1)

