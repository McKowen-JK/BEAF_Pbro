#started on 150921, plot DESeq2 results for looking for changes in promoter read counts in the 50bp max pause window found with Greg's script
#Based off of /media/jacobmtome/HDD/JMT_Lis_Lab/2015_07_PRO_seq/Deseq2/Fixed_Minus_Strand/Compare_Plot_Results.R
#use the results from /media/jacobmtome/HDD/JMT_Lis_Lab/2015_07_PRO_seq/Deseq2/MaxPauseWindow/DESeq2_JMT_Promoter_MaxPauseWindow.R


setwd('/home/chart/PROseq2/DESeq/')
GeneList = read.table("/home/chart/scripts/gene_lists/finalsetMJG_BedFormat.dat", header=T)


# can change padj to be more or less stringent 0.1 to 0.001
# second padj line determines if gene is not changed in 2nd dataset
padj_cutoff = 0.05
padj_cutoff_other = 0.5   #to call a gene as not changed in one condition, require it to have a decently high p value

BF_LZ_GeneBody = read.csv('/home/chart/PROseq2/DESeq/DESeq_gene_lists/BF_LZ_GeneBody_DESeq2.csv')
PB_LZ_GeneBody = read.csv('/home/chart/PROseq2/DESeq/DESeq_gene_lists/PB_LZ_GeneBody_DESeq2.csv')


BF_PB = merge(BF_LZ_GeneBody, PB_LZ_GeneBody, by  ='X', suffixes = c('_BF', '_PB'))

BF_All = subset(BF_PB, padj_BF < padj_cutoff)
BF_All_genelist = merge(GeneList, BF_All, by.x = 'name', by.y = 'X')
BF_Only = subset(BF_PB, padj_BF < padj_cutoff & padj_PB > padj_cutoff_other)
BF_Only_genelist = merge(GeneList, BF_Only, by.x = 'name', by.y = 'X')

PB_All = subset(BF_PB, padj_PB < padj_cutoff)
PB_All_genelist = merge(GeneList, PB_All, by.x = 'name', by.y = 'X')
PB_Only = subset(BF_PB, padj_PB < padj_cutoff & padj_BF > padj_cutoff_other)
PB_Only_genelist = merge(GeneList, PB_Only, by.x = 'name', by.y = 'X')

BF_and_PB = subset(BF_PB, padj_BF < padj_cutoff & padj_PB < padj_cutoff)
BF_and_PB_genelist = merge(GeneList, BF_and_PB, by.x = 'name', by.y = 'X')

#write.csv(as.data.frame(BF_Only_genelist), file='BF_Only_SignificantlyChanged.csv')
#write.csv(as.data.frame(PB_Only_genelist), file='PB_Only_SignificantlyChanged.csv')

plot(BF_PB[,2], BF_PB[,3], log='x', xlab='Normalized Mean Expression (BEAF)', ylab='log2 Fold Change (BEAF)', xlim=c(0.1,500000), cex.lab=1.25, cex.axis=1.25, ylim=c(-2.75, 2.75))#plot the log fold change and mean expression (MA plot)
points(BF_All[,2], BF_All[,3], col='red', lwd=3)  #now add points for the significant genes from BEAF
#points(WT_Only_Promoter[,2], WT_Only_Promoter[,3], col='green', lwd=3)  #now add points for the significant genes from the aptamer
#points(Apt_WT_Promoter[,2], Apt_WT_Promoter[,3], col='goldenrod', lwd=3)  #now add points for the significant genes from the aptamer
title(main='Effect of BEAF RNAi on \ngene body read counts: padj < 0.05', cex=1.5)
text(.06,-2.4,labels=bquote('BEAF All' ==.(dim(BF_All)[1])), col='red', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it

plot(BF_PB[,2], BF_PB[,3], log='x', xlab='Normalized Mean Expression (BEAF)', ylab='log2 Fold Change (BEAF)', xlim=c(0.1,500000), cex.lab=1.5, cex.axis=1.25, ylim=c(-2.75, 2.75))#plot the log fold change and mean expression (MA plot)
points(BF_Only[,2], BF_Only[,3], col='red', lwd=3)  #now add points for the significant genes from BEAF
title(main='Effect of BEAF RNAi (minus PB) on \ngene body read counts: padj < 0.05', cex=1.5)
text(.06,-2.4,labels=bquote('BEAF Only' ==.(dim(BF_Only)[1])), col='red', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
abline(h=0, col = 'blue')

plot(BF_PB[,8], BF_PB[,9], log='x', xlab='Normalized Mean Expression (PBro)', ylab='log2 Fold Change (PBro)', xlim=c(0.1,500000), cex.lab=1.25, cex.axis=1.25, ylim=c(-2.75, 2.75))#plot the log fold change and mean expression (MA plot)
points(PB_All[,8], PB_All[,9], col='red', lwd=3)  #now add points for the significant genes from BEAF
title(main='Effect of PBro RNAi on \ngene body read counts: padj < 0.05', cex=1.5)
text(.06,-2.4,labels=bquote('PBro All' ==.(dim(PB_All)[1])), col='red', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
abline(h=0, col = 'blue')

plot(BF_PB[,8], BF_PB[,9], log='x', xlab='Normalized Mean Expression (PBro)', ylab='log2 Fold Change (PBro)', xlim=c(0.1,500000), cex.lab=1.25, cex.axis=1.25, ylim=c(-2.75, 2.75))#plot the log fold change and mean expression (MA plot)
points(PB_Only[,8], PB_Only[,9], col='red', lwd=3)  #now add points for the significant genes from BEAF
title(main='Effect of PBro RNAi (minus BF) on \ngene body read counts: padj < 0.05', cex=1.5)
text(.06,-2.4,labels=bquote('PBro All' ==.(dim(PB_Only)[1])), col='red', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
abline(h=0, col = 'blue')

plot(BF_PB[,2], BF_PB[,3], log='x', xlab='Normalized Mean Expression (BEAF)', ylab='log2 Fold Change (BEAF)', xlim=c(0.1,500000), cex.lab=1.25, cex.axis=1.25, ylim=c(-2.75, 2.75))#plot the log fold change and mean expression (MA plot)
points(BF_and_PB[,2], BF_and_PB[,3], col='red', lwd=3)  #now add points for the significant genes from BEAF
title(main='Effect of BEAF RNAi (with PB) on \ngene body read counts: padj < 0.05', cex=1.5)
text(.06,-2.4,labels=bquote('BEAF and PBro' ==.(dim(BF_and_PB)[1])), col='red', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
abline(h=0, col = 'blue')

plot(BF_PB[,8], BF_PB[,9], log='x', xlab='Normalized Mean Expression (PBro)', ylab='log2 Fold Change (PBro)', xlim=c(0.1,500000), cex.lab=1.25, cex.axis=1.25, ylim=c(-2.75, 2.75))#plot the log fold change and mean expression (MA plot)
points(BF_and_PB[,8], BF_and_PB[,9], col='red', lwd=3)  #now add points for the significant genes from BEAF
title(main='Effect of PBro RNAi (with BF) on \ngene body read counts: padj < 0.05', cex=1.5)
text(.06,-2.4,labels=bquote('PBro and BEAF' ==.(dim(BF_and_PB)[1])), col='red', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
abline(h=0, col = 'blue')

#Now, use the data in the same tables to plot the same data, but with the counts and fold change in S2 cells (so, just change the columns used for x and y)
#151014_MA_forSeminar_WTversion_1
plot(Apt_S2_Promoter[,8], Apt_S2_Promoter[,9], log='x', xlab='Normalized Mean Expression (WT)', ylab='log2 Fold Change (WT)', xlim=c(0.1,10000), cex.lab=1.5, cex.axis=1.25, ylim=c(-2.5, 2.5))#plot the log fold change and mean expression (MA plot)
points(Apt_Only_Promoter[,8], Apt_Only_Promoter[,9], col='red', lwd=3)  #now add points for the significant genes from the aptamer
points(WT_Only_Promoter[,8], WT_Only_Promoter[,9], col='green', lwd=3)  #now add points for the significant genes from the aptamer

points(Apt_WT_Promoter[,8], Apt_WT_Promoter[,9], col='goldenrod', lwd=3)  #now add points for the significant genes from the aptamer
title(main='Effect of Copper Induction on Promoter Read Counts\nin WT Cells: padj < 0.01', cex=1.5)
text(.06,-2.4,labels=bquote('Aptamer Only' ==.(dim(Apt_Only_Promoter)[1])), col='red', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
text(.06,-2,labels=bquote('Changed in Both' ==.(dim(Apt_WT_Promoter)[1])), col='goldenrod', pos=4)
text(.06,-1.6,labels=bquote('WT Only' ==.(dim(WT_Only_Promoter)[1])), col='darkgreen', pos=4)
abline(h=0, col = 'blue')


Apt_Promoter_All = subset(Apt_S2_Promoter, padj_A < padj_cutoff)
WT_Promoter = subset(Apt_S2_Promoter, padj_WT < padj_cutoff)


plot(Apt_S2_Promoter[,2], Apt_S2_Promoter[,3], log='x', xlab='Normalized Mean Expression (Apt)', ylab='log2 Fold Change (Apt)', xlim=c(0.1,10000), cex.lab=1.5, cex.axis=1.25, ylim=c(-2.5, 2.5))#plot the log fold change and mean expression (MA plot)
points(Apt_Promoter_All[,2], Apt_Promoter_All[,3], col='red', lwd=3)  #now add points for the significant genes from the aptamer
points(WT_Promoter[,2], WT_Promoter[,3], col='green', lwd=3)  #now add points for the significant genes from the aptamer
points(Apt_WT_Promoter[,2], Apt_WT_Promoter[,3], col='goldenrod', lwd=3)  #now add points for the significant genes from the aptamer
title(main='Effect of Copper Induction on Promoter Read Counts\nin Aptamer Containing Cells: padj <0.01', cex=1.5)
text(.06,-2.4,labels=bquote('Changed With Aptamer' ==.(dim(Apt_Promoter_All)[1])), col='red', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
text(.06,-2,labels=bquote('Changed in Both' ==.(dim(Apt_WT_Promoter)[1])), col='goldenrod', pos=4)
text(.06,-1.6,labels=bquote('Changed in WT' ==.(dim(WT_Promoter)[1])), col='green', pos=4)
abline(h=0, col = 'blue')



#################################################  Now plot the genes that change in gene body vs promoter in the aptamer cells

#I GET A DIFFERENT NUMBER OF DIFFERENTIALLY EXPRESSED GENES IN THE PROMOTER BECAUSE I AM ONLY INCLUDING GENES THAT DO NOT GIVE A PADJ OF NA FOR BOTH
#A gene will give a padj of NA (and only NA there) if it did not pass the automatic independent filtering due to low normalized mean expression
#So, I'm not working with all of the genes with a promoter padj < 0.1 here... those that did not have enough data in the gene body to determine if they are differentially expressed are excluded
#the same is true above

Apt_Combined = merge(Apt_Promoter, Apt_GeneBody, by  ='X', suffixes = c('_P', '_GB'))

Apt_Promoter_Any = subset(Apt_Combined, padj_P < padj_cutoff)


Apt_Promoter_Only = subset(Apt_Combined, padj_P < padj_cutoff & padj_GB > padj_cutoff_other)
Apt_GB_Only = subset(Apt_Combined, padj_P > padj_cutoff_other & padj_GB < padj_cutoff)
Apt_GBandPromoter = subset(Apt_Combined, padj_P < padj_cutoff & padj_GB < padj_cutoff)

plot(Apt_Combined[,2], Apt_Combined[,3], log='x', xlab='Normalized Mean Expression (Apt)', ylab='log2 Fold Change (Apt)', xlim=c(0.1,10000), cex.lab=1.25, cex.axis=1.25, ylim=c(-2.5, 2.5))#plot the log fold change and mean expression (MA plot)
points(Apt_Promoter_Only[,2], Apt_Promoter_Only[,3], col='red', lwd=3)  #now add points for the significant genes from the aptamer
points(Apt_GB_Only[,2], Apt_GB_Only[,3], col='green', lwd=3)  #now add points for the significant genes from the aptamer
points(Apt_GBandPromoter[,2], Apt_GBandPromoter[,3], col='goldenrod', lwd=3)  #now add points for the significant genes from the aptamer
title(main='Effect of Copper Induction on Promoter Read Counts\nin Aptamer Containing Cells', cex=1.5)
text(.06,-2.4,labels=bquote('Changed in Promoter Only' ==.(dim(Apt_Promoter_Only)[1])), col='red', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
text(.06,-2,labels=bquote('Changed in Both' ==.(dim(Apt_GBandPromoter)[1])), col='goldenrod', pos=4)
text(.06,-1.6,labels=bquote('Changed in Gene Body Only' ==.(dim(Apt_GB_Only)[1])), col='green', pos=4)
abline(h=0, col = 'blue')




#plot the same differential expression data with the gene body MA data
plot(Apt_Combined[,8], Apt_Combined[,9], log='x', xlab='Normalized Mean Expression (Apt)', ylab='log2 Fold Change (Apt)', xlim=c(0.1,10000), cex.lab=1.25, cex.axis=1.25, ylim=c(-2.5, 2.5))#plot the log fold change and mean expression (MA plot)
points(Apt_Promoter_Only[,8], Apt_Promoter_Only[,9], col='red', lwd=3)  #now add points for the significant genes from the aptamer
points(Apt_GB_Only[,8], Apt_GB_Only[,9], col='green', lwd=3)  #now add points for the significant genes from the aptamer
points(Apt_GBandPromoter[,8], Apt_GBandPromoter[,9], col='goldenrod', lwd=3)  #now add points for the significant genes from the aptamer
title(main='Effect of Copper Induction on Gene Body Read Counts\nin Aptamer Containing Cells: padj < 0.01', cex=1.5)
text(.06,-2.4,labels=bquote('Changed in Promoter Only' ==.(dim(Apt_Promoter_Only)[1])), col='red', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
text(.06,-2,labels=bquote('Changed in Both' ==.(dim(Apt_GBandPromoter)[1])), col='goldenrod', pos=4)
text(.06,-1.6,labels=bquote('Changed in Gene Body Only' ==.(dim(Apt_GB_Only)[1])), col='green', pos=4)
abline(h=0, col = 'blue')



#now use all significantly differentially expressed genes, regardless of whether the other set was able to be called as sig or not
Apt_Promoter_All = subset(Apt_Combined, padj_P < padj_cutoff)
Apt_GB_All = subset(Apt_Combined, padj_GB < padj_cutoff)


plot(Apt_Combined[,8], Apt_Combined[,9], log='x', xlab='Normalized Mean Expression (Apt)', ylab='log2 Fold Change (Apt)', xlim=c(0.1,10000), cex.lab=1.25, cex.axis=1.25, ylim=c(-2.5, 2.5))#plot the log fold change and mean expression (MA plot)
points(Apt_Promoter_All[,8], Apt_Promoter_All[,9], col='red', lwd=3)  #now add points for the significant genes from the aptamer
points(Apt_GB_All[,8], Apt_GB_All[,9], col='green', lwd=3)  #now add points for the significant genes from the aptamer
points(Apt_GBandPromoter[,8], Apt_GBandPromoter[,9], col='goldenrod', lwd=3)  #now add points for the significant genes from the aptamer
title(main='Effect of Copper Induction on Gene Body Read Counts\nin Aptamer Containing Cells: padj < 0.01', cex=1.5)
text(.1,2.35,labels=bquote('Changed in Promoter' ==.(dim(Apt_Promoter_All)[1])), col='red', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
text(.1,2.1,labels=bquote('Changed in Both' ==.(dim(Apt_GBandPromoter)[1])), col='goldenrod', pos=4)
text(.1,1.85,labels=bquote('Changed in Gene Body' ==.(dim(Apt_GB_All)[1])), col='green', pos=4)

plot(Apt_Combined[,2], Apt_Combined[,3], log='x', xlab='Normalized Mean Expression (Apt)', ylab='log2 Fold Change (Apt)', xlim=c(0.1,10000), cex.lab=1.25, cex.axis=1.25, ylim=c(-2.5, 2.5))#plot the log fold change and mean expression (MA plot)
points(Apt_Promoter_All[,2], Apt_Promoter_All[,3], col='red', lwd=3)  #now add points for the significant genes from the aptamer
points(Apt_GB_All[,2], Apt_GB_All[,3], col='green', lwd=3)  #now add points for the significant genes from the aptamer
points(Apt_GBandPromoter[,2], Apt_GBandPromoter[,3], col='goldenrod', lwd=3)  #now add points for the significant genes from the aptamer
title(main='Effect of Copper Induction on Promoter Read Counts\nin Aptamer Containing Cells', cex=1.5)
text(.1,2.35,labels=bquote('Changed in Promoter' ==.(dim(Apt_Promoter_All)[1])), col='red', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
text(.1,2.1,labels=bquote('Changed in Both' ==.(dim(Apt_GBandPromoter)[1])), col='goldenrod', pos=4)
text(.1,1.85,labels=bquote('Changed in Gene Body' ==.(dim(Apt_GB_All)[1])), col='green', pos=4)



#################################################  Now plot the genes that change in gene body vs promoter in the wild type cells

S2_Promoter_GeneBody = merge(S2_Promoter, S2_GeneBody, by  ='X', suffixes = c('_P', '_GB'))


S2_Promoter_all = subset(S2_Promoter_GeneBody, padj_P < padj_cutoff)
S2_GeneBody_all = subset(S2_Promoter_GeneBody, padj_GB < padj_cutoff)
S2_Promoter_AndGeneBody = subset(S2_Promoter_GeneBody, padj_P < padj_cutoff & padj_GB < padj_cutoff)

plot(S2_Promoter_GeneBody[,2], S2_Promoter_GeneBody[,3], log='x', xlab='Normalized Mean Expression (Apt)', ylab='log2 Fold Change (Apt)', xlim=c(0.1,10000), cex.lab=1.25, cex.axis=1.25, ylim=c(-2.5, 2.5))#plot the log fold change and mean expression (MA plot)
points(S2_Promoter_all[,2], S2_Promoter_all[,3], col='red', lwd=3)  #now add points for the significant genes from the aptamer
points(S2_GeneBody_all[,2], S2_GeneBody_all[,3], col='green', lwd=3)  #now add points for the significant genes from the aptamer
points(S2_Promoter_AndGeneBody[,2], S2_Promoter_AndGeneBody[,3], col='goldenrod', lwd=3)  #now add points for the significant genes from the aptamer
title(main='Effect of Copper Induction on Promoter Read Counts\nin WT Cells', cex=1.5)
text(.1,2.35,labels=bquote('Changed Promoter' ==.(dim(S2_Promoter_all)[1])), col='red', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
text(.1,2.1,labels=bquote('Changed in Both' ==.(dim(S2_Promoter_AndGeneBody)[1])), col='goldenrod', pos=4)
text(.1,1.85,labels=bquote('Changed in Gene Body' ==.(dim(S2_GeneBody_all)[1])), col='green', pos=4)


#################################################  Make an MA Plot of the up and down regulated genes in Apt containing cells, that are only up or down due to the aptamer

Apt_Changed_Promoter = subset(Apt_S2_Promoter, padj_A < padj_cutoff & padj_WT > padj_cutoff_other)

#write.csv(as.data.frame(Apt_Changed_Promoter), file='Apt_SignificantlyChanged_Promoter_NoChangeWT.csv')

Apt_Only_Promoter_Up = subset(Apt_Changed_Promoter, log2FoldChange_A > 0)
Apt_Only_Promoter_Down = subset(Apt_Changed_Promoter, log2FoldChange_A < 0)

plot(Apt_S2_Promoter[,2], Apt_S2_Promoter[,3], log='x', xlab='Normalized Mean Expression (Apt)', ylab='log2 Fold Change (Apt)', xlim=c(0.1,10000), cex.lab=1.25, cex.axis=1.25, ylim=c(-2.5, 2.5))#plot the log fold change and mean expression (MA plot)
points(Apt_Only_Promoter_Up[,2], Apt_Only_Promoter_Up[,3], col='red', lwd=3)  #now add points for the significant genes from the aptamer

points(Apt_WT_Promoter[,2], Apt_WT_Promoter[,3], col='darkgoldenrod1', lwd=3)  #now add points for the significant genes from the aptamer

points(Apt_Only_Promoter_Down[,2], Apt_Only_Promoter_Down[,3], col='green', lwd=3)  #now add points for the significant genes from the aptamer
title(main='Genes with Significantly Changed Promoter Read Counts\nin Aptamer Containing Cells Only: padj > 0.01', cex=1.5)
text(.06,2.4,labels=bquote('Increased Promoter Counts' ==.(dim(Apt_Only_Promoter_Up)[1])), col='red', pos=4) #Pos=4 will make it place text to the right of the coordinate I gave, rather than centering on it
text(.06,2.0,labels=bquote('Changed in WT Too' ==.(dim(Apt_WT_Promoter)[1])), col='darkgoldenrod1', pos=4)
text(.06,1.6,labels=bquote('Decreased Promoter Counts' ==.(dim(Apt_Only_Promoter_Down)[1])), col='green', pos=4)
abline(h=0, col = 'blue')





###################See which genes are changed by looking at the max pause window like this vs just around the promoter
Old_AptOnly = read.csv('/media/jacobmtome/HDD/JMT_Lis_Lab/2015_07_PRO_seq/Deseq2/Fixed_Minus_Strand/AptOnly_SignificantlyChanged_Promoter.csv')
Old_New = merge(Old_AptOnly, Apt_Only_Promoter, by.x = 'X.1', by.y = 'X')
#OF THE 128 GENES THAT CHANGED USING THE LARGE WINDOW, 120 WERE CHANGED IN THE MAX PAUSE WINDOW.  SO, FOR THE MOST PART, i HAVE JUST ADDED GENES BY REDUCING NOISE
