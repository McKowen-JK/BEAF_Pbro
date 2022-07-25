#started on 150921, plot DESeq2 results for looking for changes in promoter read counts in the 50bp max pause window found with Greg's script
#Based off of /media/jacobmtome/HDD/JMT_Lis_Lab/2015_07_PRO_seq/Deseq2/Fixed_Minus_Strand/Compare_Plot_Results.R
#use the results from /media/jacobmtome/HDD/JMT_Lis_Lab/2015_07_PRO_seq/Deseq2/MaxPauseWindow/DESeq2_JMT_Promoter_MaxPauseWindow.R

setwd('/home/chart/PROseq2/DESeq/')
GeneList = read.csv("/home/chart/scripts/gene_lists/finalsetMJG_BedFormat_LacZ_50bpPause.csv", header=T)

# can change padj to be more or less stringent 0.1 to 0.001
# second padj line determines if gene is not changed in 2nd dataset
padj_cutoff = 0.05
padj_cutoff_other = 0.5   #to call a gene as not changed in one condition, require it to have a decently high p value
log2FoldChange_cutoff = 0

BF_LZ_Pause = read.csv('/home/chart/PROseq2/DESeq/DESeq_gene_lists/BF_LZ_50bpPause_DESeq2.csv')
PB_LZ_Pause = read.csv('/home/chart/PROseq2/DESeq/DESeq_gene_lists/PB_LZ_50bpPause_DESeq2.csv')

BF_PB = merge(BF_LZ_Pause, PB_LZ_Pause, by  ='X', suffixes = c('_BF', '_PB'))

BF_All = subset(BF_PB, padj_BF < padj_cutoff)
BF_All_genelist = merge(GeneList, BF_All, by.x = 'name', by.y = 'X')
BF_up = subset(BF_All, log2FoldChange_BF > log2FoldChange_cutoff)
BF_up_genelist = merge(GeneList, BF_up, by.x = 'name', by.y = 'X')
BF_dn = subset(BF_All, log2FoldChange_BF < log2FoldChange_cutoff)
BF_dn_genelist = merge(GeneList, BF_dn, by.x = 'name', by.y = 'X')

PB_All = subset(BF_PB, padj_PB < padj_cutoff)
PB_All_genelist = merge(GeneList, PB_All, by.x = 'name', by.y = 'X')
PB_up = subset(PB_All, log2FoldChange_PB > log2FoldChange_cutoff)
PB_up_genelist = merge(GeneList, PB_up, by.x = 'name', by.y = 'X')
PB_dn = subset(PB_All, log2FoldChange_PB < log2FoldChange_cutoff)
PB_dn_genelist = merge(GeneList, PB_dn, by.x = 'name', by.y = 'X')

#BF_and_PB = subset(BF_PB, padj_BF < padj_cutoff & padj_PB < padj_cutoff)
#BF_and_PB_genelist = merge(GeneList, BF_and_PB, by.x = 'name', by.y = 'X')

#write.csv(as.data.frame(BF_Only_genelist), file='BF_Only_SignificantlyChanged.csv')
#write.csv(as.data.frame(PB_Only_genelist), file='PB_Only_SignificantlyChanged.csv')

#plot(BF_PB[,2], BF_PB[,3], log='x', xlab='Normalized Mean Expression (BEAF)', ylab='log2 Fold Change (BEAF)', xlim=c(0.1,500000), cex.lab=1.25, cex.axis=1.25, ylim=c(-2.85, 2.85))#plot the log fold change and mean expression (MA plot)
#plot(BF_PB[,2], BF_PB[,3], pch=20, cex=0.5, col='dimgrey', log='x', xlab='Normalized Mean Expression (BEAF)', ylab='log2 Fold Change (BEAF)', xlim=c(0.1,500000), cex.lab=0.75, cex.axis=0.75, ylim=c(-2.85, 2.85), las=(1 & 2))#plot the log fold change and mean expression (MA plot)
plot(BF_PB[,2], BF_PB[,3], pch=20, cex=0.5, col='gray66', xaxt='n', yaxt='n', log='x', xlab='Average promoter reads', ylab='BEAF/LacZ promoter FC', xlim=c(0.1,500000), cex.lab=0.75, cex.axis=0.75, ylim=c(-2.85, 2.85), las=(1 & 2))#plot the log fold change and mean expression (MA plot)
axTicks<-seq(from=-1, to=5, by=1)
labels<-sapply(axTicks, function(i) as.expression(bquote(10^ .(i))))
axis(1, at=c(0.1, 1, 10, 100, 1000, 10000, 100000), labels=labels, cex.axis=0.75)
ayTicks<-seq(from=-3, to=3, by=1)
labels<-sapply(ayTicks, function(i) as.numeric(2^(i)))
labels=round(labels, digits=2)#Would trunc work as well as round?
axis(2, at=c(-3, -2, -1, 0, 1, 2, 3), labels=labels, las=2, cex.axis=0.75)
points(BF_up[,2], BF_up[,3], pch=20, cex=0.25, col='firebrick1', lwd=3)  #now add points for the significant genes from BEAF
points(BF_dn[,2], BF_dn[,3], pch=20, cex=0.25, col='dodgerblue', lwd=3)  #now add points for the significant genes from BEAF
#pch changes point size: 1: circle; 16: filled circle; 19: solid circle (larger); 20: bullet (smaller circle)
#cex changes scaling; cex.main changes title scaling; /n makes a line break; las=(1&2) rotates y-axis labels 90 degrees
title(main='Effect of BEAF RNAi on pause counts: padj < 0.05', cex.main=0.75)
text(.06,2.1,labels=bquote('up' ==.(dim(BF_up)[1])), col='firebrick1', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
text(.06,-2.6,labels=bquote('down' ==.(dim(BF_dn)[1])), col='dodgerblue', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it


plot(BF_PB[,8], BF_PB[,9], pch=20, cex=0.5, col='gray66', xaxt='n', yaxt='n', log='x', xlab='Average promoter reads', ylab='Pbro/LacZ promoter FC', xlim=c(0.1,500000), cex.lab=0.75, cex.axis=0.75, ylim=c(-3.75, 2.75), las=(1 & 2))#plot the log fold change and mean expression (MA plot)
axTicks<-seq(from=-1, to=5, by=1)
labels<-sapply(axTicks, function(i) as.expression(bquote(10^ .(i))))
axis(1, at=c(0.1, 1, 10, 100, 1000, 10000, 100000), labels=labels, cex.axis=0.75)
ayTicks<-seq(from=-4, to=3, by=1)
labels<-sapply(ayTicks, function(i) as.numeric(2^(i)))
labels=round(labels, digits=2)
axis(2, at=c(-4, -3, -2, -1, 0, 1, 2, 3), labels=labels, las=2, cex.axis=0.75)
points(PB_up[,8], PB_up[,9], pch=20, cex=0.5, col='magenta', lwd=3)  #now add points for the significant genes from BEAF
points(PB_dn[,8], PB_dn[,9], pch=20, cex=0.5, col='seagreen', lwd=3)  #now add points for the significant genes from BEAF
title(main='Effect of PBro RNAi on pause counts: padj < 0.05', cex.main=0.75)
text(.06,1.9,labels=bquote('up' ==.(dim(PB_up)[1])), col='magenta', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
text(.06,-3.3,labels=bquote('down' ==.(dim(PB_dn)[1])), col='seagreen', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
#abline(h=0, col = 'blue')
#look up colors in R for different colors


#Now make 2 color gene body MA plots
#Note that some temporary files have the same name as for the promoter region script
#So run separately (promoter first, then gene body) to avoid overwriting problems

setwd('/home/chart/PROseq2/DESeq/')
GeneList = read.table("/home/chart/scripts/gene_lists/finalsetMJG_BedFormat.dat", header=T)

# can change padj to be more or less stringent 0.1 to 0.001
# second padj line determines if gene is not changed in 2nd dataset
padj_cutoff = 0.05
padj_cutoff_other = 0.5   #to call a gene as not changed in one condition, require it to have a decently high p value
log2FoldChange_cutoff = 0

BF_LZ_GeneBody = read.csv('/home/chart/PROseq2/DESeq/DESeq_gene_lists/BF_LZ_GeneBody_DESeq2.csv')
PB_LZ_GeneBody = read.csv('/home/chart/PROseq2/DESeq/DESeq_gene_lists/PB_LZ_GeneBody_DESeq2.csv')

BF_PB = merge(BF_LZ_GeneBody, PB_LZ_GeneBody, by  ='X', suffixes = c('_BF', '_PB'))

#BF_All = subset(BF_PB, padj_BF < padj_cutoff)
#BF_All_genelist = merge(GeneList, BF_All, by.x = 'name', by.y = 'X')
BF_up = subset(BF_PB, log2FoldChange_BF > log2FoldChange_cutoff & padj_BF < padj_cutoff)
BF_up_genelist = merge(GeneList, BF_up, by.x = 'name', by.y = 'X')
BF_dn = subset(BF_PB, log2FoldChange_BF < log2FoldChange_cutoff & padj_BF < padj_cutoff)
BF_dn_genelist = merge(GeneList, BF_dn, by.x = 'name', by.y = 'X')

#PB_All = subset(BF_PB, padj_PB < padj_cutoff)
#PB_All_genelist = merge(GeneList, PB_All, by.x = 'name', by.y = 'X')
PB_up = subset(BF_PB, log2FoldChange_PB > log2FoldChange_cutoff & padj_PB < padj_cutoff)
PB_up_genelist = merge(GeneList, PB_up, by.x = 'name', by.y = 'X')
PB_dn = subset(BF_PB, log2FoldChange_PB < log2FoldChange_cutoff & padj_PB < padj_cutoff)
PB_dn_genelist = merge(GeneList, PB_dn, by.x = 'name', by.y = 'X')

#BF_and_PB = subset(BF_PB, padj_BF < padj_cutoff & padj_PB < padj_cutoff)
#BF_and_PB_genelist = merge(GeneList, BF_and_PB, by.x = 'name', by.y = 'X')

#write.csv(as.data.frame(BF_Only_genelist), file='BF_Only_SignificantlyChanged.csv')
#write.csv(as.data.frame(PB_Only_genelist), file='PB_Only_SignificantlyChanged.csv')

plot(BF_PB[,2], BF_PB[,3], pch=20, cex=0.5, col='gray66', xaxt='n', yaxt='n', log='x', xlab='Average gene body reads', ylab='BEAF/LacZ gene body FC', xlim=c(0.1,500000), cex.lab=0.75, cex.axis=0.75, ylim=c(-2.85, 2.85), las=(1 & 2))#plot the log fold change and mean expression (MA plot)
axTicks<-seq(from=-1, to=5, by=1)
labels<-sapply(axTicks, function(i) as.expression(bquote(10^ .(i))))
axis(1, at=c(0.1, 1, 10, 100, 1000, 10000, 100000), labels=labels, cex.axis=0.75)
ayTicks<-seq(from=-3, to=3, by=1)
labels<-sapply(ayTicks, function(i) as.numeric(2^(i)))
labels=round(labels, digits=2)
axis(2, at=c(-3, -2, -1, 0, 1, 2, 3), labels=labels, las=2, cex.axis=0.75)
points(BF_up[,2], BF_up[,3], pch=20, cex=0.25, col='firebrick1', lwd=3)  #now add points for the significant genes from BEAF
points(BF_dn[,2], BF_dn[,3], pch=20, cex=0.25, col='dodgerblue', lwd=3)  #now add points for the significant genes from BEAF
# col='goldenrod', col='darkorange', col='orangered', col='dimgrey'
title(main='Effect of BEAF RNAi on body counts: padj < 0.05', cex.main=0.75)
text(.06,2.1,labels=bquote('up' ==.(dim(BF_up)[1])), col='firebrick1', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
text(.06,-2.6,labels=bquote('down' ==.(dim(BF_dn)[1])), col='dodgerblue', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it


plot(BF_PB[,8], BF_PB[,9], pch=20, cex=0.5, col='gray66', xaxt='n', yaxt='n', log='x', xlab='Average gene body reads', ylab='Pbro/LacZ gene body FC', xlim=c(0.1,500000), cex.lab=0.75, cex.axis=0.75, ylim=c(-3.75, 2.75), las=(1 & 2))#plot the log fold change and mean expression (MA plot)
axTicks<-seq(from=-1, to=5, by=1)
labels<-sapply(axTicks, function(i) as.expression(bquote(10^ .(i))))
axis(1, at=c(0.1, 1, 10, 100, 1000, 10000, 100000), labels=labels, cex.axis=0.75)
ayTicks<-seq(from=-4, to=3, by=1)
labels<-sapply(ayTicks, function(i) as.numeric(2^(i)))
labels=round(labels, digits=2)
axis(2, at=c(-4, -3, -2, -1, 0, 1, 2, 3), labels=labels, las=2, cex.axis=0.75)
points(PB_up[,8], PB_up[,9], pch=20, cex=0.25, col='magenta', lwd=3)  #now add points for the significant genes from BEAF
points(PB_dn[,8], PB_dn[,9], pch=20, cex=0.25, col='seagreen', lwd=3)  #now add points for the significant genes from BEAF
title(main='Effect of PBro RNAi on body counts: padj < 0.05', cex.main=0.75)
text(.06,1.9,labels=bquote('up' ==.(dim(PB_up)[1])), col='magenta', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
text(.06,-3.3,labels=bquote('down' ==.(dim(PB_dn)[1])), col='seagreen', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
#abline(h=0, col = 'blue')
#look up colors in R for different colors

