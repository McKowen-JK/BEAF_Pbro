#!/bin/bash
## Script to compute matrix and plot a heatmap with deeptools
##to run this script all I did was copy and paste it into the terminal

##This works, but in the original .csv I deleted some columns (12 max), rearranged columns, and converted to tab-delimited by saving in LibreOffice calc as .csv, checking 'Edit filter settings', click Save, select 'Field delimiter' {Tab}, then Save; also had to start the header row with # (#Chrom) 


##This is a reference for writing the script, from http://deeptools.readthedocs.io/en/latest/content/tools/computeMatrix.html example 1
computeMatrix reference-point \ # choose the mode
       --referencePoint TSS \ # alternatives: TES, center (TSS is default)
       -b 1000 -a 1000 \ # define the region you are interested in
       -R /home/chart/modENCODE_bigWigs/oli_bf_chipSeq_bedgraph/SRR1042411.bed.gz \
       -S /home/chart/modENCODE_bigWigs/BEAF32_ChIPseq.GSM1278639.bigWig  \
       --skipZeros \
       -o bfMatrix.mat.g \ # to be used with plotHeatmap and plotProfile
       --outFileSortedRegions regions1_H3K4me3_l2r_genes.bed

plotHeatmap -m bfMatrix.mat.gz \
      -out ExampleHeatmap1.png \

## Set some enviroment variables.
##the bedgraph directory is where you want the bedgraph files to go
##example: bedgraphDIR=/home/chart/modENCODE_bigWigs/modEncode_bedgraph
##the macs2 directory is where you want the peak list files to go
##example: macs2DIR=/home/chart/modENCODE_bigWigs/modEncode_macs2
##starting direction is where the bigwig files are
##example: STARTDIR=/home/chart/modENCODE_bigWigs


## Script using shortened names to compute matrix and plot a heatmap with deeptools
##to run this script all I did was copy and paste it into the terminal

##This works, but in the original .csv I deleted some columns (12 max), rearranged columns, and converted to tab-delimited by saving in LibreOffice calc as .csv, checking 'Edit filter settings', click Save, select 'Field delimiter' {Tab}, then Save; also had to start the header row with # (#Chrom) 


##modified for PRO-seq data (LacZ RNAi active genes arranged by BEAF ChIP-seq)


#!/bin/bash
##give gene lists (tab-delimited csv and bed) and bigwig files shortened names
bf_chip=/home/chart/modENCODE_bigWigs/oli_bf_chipSeq_bedgraph/test_oli/beaf_chSeq_srr1042411_treat_pileup.bigWig
bf_peaks=/home/chart/modENCODE_bigWigs/oli_bf_chipSeq_bedgraph/beaf_SRR1042411_macs2_peaks4.csv
LZact=/home/chart/PROseq2/pausing/activePausedGeneLists/2021apr_Genes_LZ_active4b.csv
LZplus=/home/chart/PROseq2/bigwig/LacZ_12_normed_plus_noMito.bw
LZminus=/home/chart/PROseq2/bigwig/LacZ_12_normed_minus_noMito.bw

##make a matrix of LacZ RNAi active genes (6069) based on oli BEAF ChIP-seq signal and output a bed file sorted by Bf ChIP-seq intensity
computeMatrix reference-point -b 300 -a 300 -R $LZact -S $bf_chip --missingDataAsZero --skipZeros -o /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_BfChIPseq300.mat.gz

##make a red-blue heatmap of LacZ RNAi active genes (6069) sorted by Bf intensity
plotHeatmap -m /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_BfChIPseq300.mat.gz --colorMap RdBu --outFileSortedRegions /home/chart/MNase_seq2/MNase_120_180_bw/LZact_sorted_BfChIPseq300.bed -out /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_BfChIPseq300.png 
LZact_ChIP=/home/chart/MNase_seq2/MNase_120_180_bw/LZact_sorted_BfChIPseq300.bed


##make separate matrices of LZact_ChIP genes sorted by PRO-seq LacZ RNAi signal for plus and minus bigwigs
computeMatrix reference-point -b 300 -a 300 -R $LZact_ChIP -S $LZplus --missingDataAsZero --skipZeros -o /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300p.mat.gz

computeMatrix reference-point -b 300 -a 300 -R $LZact_ChIP -S $LZminus --missingDataAsZero --skipZeros -o /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300m.mat.gz

##maybe I only need the plus heatmap: make a red-blue or red heatmap of LacZ RNAi active genes (6069) plus-strand PRO-seq signal sorted by Bf intensity
plotHeatmap -m /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300p.mat.gz --sortRegions no --colorMap RdBu -out /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300p_RdBu.png 
plotHeatmap -m /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300p.mat.gz --sortRegions no --colorMap Reds -out /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300p_red.png 

##maybe I only need the plus heatmap: make a red-blue or red heatmap of LacZ RNAi active genes (6069) plus-strand PRO-seq signal sorted by Bf intensity, default sorting (--zMin and --zMax to set scale; --colorMap Reds_r to reverse the scale)
plotHeatmap -m /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300p.mat.gz --colorMap RdBu -out /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300p_RdBu2.png 
plotHeatmap -m /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300p.mat.gz --colorMap Reds -out /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300p_red2.png 

plotHeatmap -m /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300p.mat.gz --zMin 0 --zMax 0.7 --colorMap RdBu_r -out /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300p_RdBu2b.png 
plotHeatmap -m /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300p.mat.gz --zMin 0 --zMax 0.7 --colorMap Reds_r -out /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300p_red2b.png 


##combine plus and minus matrices using rbind
##Doesn't work, plus is top of matrix, minus is bottom, instead of combining by genes; heatmap is wrong
computeMatrixOperations rbind -m /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300p.mat.gz /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300m.mat.gz -o /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300.mat.gz

##make a red-blue or red heatmap of LacZ RNAi active genes (6069) PRO-seq signal sorted by Bf intensity
##no good, the matrix is not right
plotHeatmap -m /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300.mat.gz --sortRegions no --colorMap RdBu -out /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300.png 
plotHeatmap -m /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300.mat.gz --sortRegions no --colorMap Reds -out /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300b.png 


##make a matrix of LZact_ChIP genes sorted by PRO-seq LacZ RNAi signal (is this right or do I need to make plus and minus separately and then do something like rbind?)
##Doesn't work, plus and minus are plotted as separate heatmaps
computeMatrix reference-point -b 300 -a 300 -R $LZact_ChIP -S $LZplus $LZminus --missingDataAsZero --skipZeros -o /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300pm.mat.gz

##make a red-blue heatmap of LacZ RNAi active genes (6069) PRO-seq signal sorted by Bf intensity
##Doesn't work, plus and minus are plotted as separate heatmaps
plotHeatmap -m /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300pm.mat.gz --sortRegions no --colorMap RdBu -out /home/chart/MNase_seq2/MNase_120_180_meta/2021may_LZact_PROseq_ChIPseq300pm.png 



