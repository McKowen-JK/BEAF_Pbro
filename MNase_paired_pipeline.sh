#!/bin/bash
##to run this script all I did was copy and paste it into the terminal

## Set some enviroment variables.
##the working directory is where you want the finished files to go
WORKINGDIR=/media/macproadmin/miescher/MNase_seq/test_run
##starting direction is where the fastq files are
STARTDIR=/media/macproadmin/miescher/MNase_seq
#genome direction is where the bowtie2 genome indexes are
GENOME=/media/macproadmin/miescher/genomes_etc/bowtie2_refGenome/dm3
bgToBigWig="/home/macproadmin/users/GROseq_Utils/mkbigWig/bedGraphToBigWig"

## Copy data files to the working directory.
#cp $GENOME/* $WORKINGDIR/


##put the common part of the names of zipped paired-end fastq files you want to run
data=(BF1 BF2 PB1 PB2 LZ1 LZ2)
##also make a list with the names without the replicate designation, so you can combined replicates (see below)
dataCombined=(BF PB LZ)
#####################################################################################
## Align reads.
for i in ${data[*]}
do
  echo ${i}
## Remove reads that don't pass filters.
  echo Filter out poor reads
  grep -A 3 '^@.* [^:]*:N:[^:]*:' ${i}_R1.fastq | sed '/^--$/d' | gzip > ${i}_R1.filt.fastq.gz 
  grep -A 3 '^@.* [^:]*:N:[^:]*:' ${i}_R2.fastq | sed '/^--$/d' | gzip > ${i}_R2.filt.fastq.gz
## Align using bowtie2.
  echo Bowtie2
  /home/macproadmin/Downloads/bowtie2-2.2.2/bowtie2 --no-mixed --no-discordant -p 6 -x $GENOME/Dm3.chrom.sizes -1 ${i}_R1.filt.fastq.gz -2 ${i}_R2.filt.fastq.gz -S ${i}.sam
  echo Sam to BAM
  samtools view -b -S ${i}.sam > ${i}.bam  ## Don't sort here ... removes mate-pair order.
  echo BAM to sorted BED
    bedtools bamtobed -bedpe -i ${i}.bam | awk 'BEGIN{OFS="\t"} ($8 > 0){print $1,$2,$6}' | sed 's/dmel_mitochondrion_genome/M/' | sort -k1,1 -k2,2n > ${i}_sorted.bed
  wc -l ${i}_sorted.bed
  echo Select 120-180bp reads
  awk 'BEGIN{OFS="\t"} ($3-$2 >120 && $3-$2 < 181) {print "chr"$1,$2,$3,"N",0,"+"}' ${i}_sorted.bed >${i}_120-181_sorted.bed
#if using sorted.bed files to make bedGaph files etc, uncomment the for and do lines and start there
# (after running the bin/bash, environment variables, data and dataCombined lines)
#for i in ${data[*]}
#do
  echo ${i}
  wc -l ${i}_120-181_sorted.bed
  echo bed to bedGraph
  genomeCoverageBed -bg -i ${i}_120-181_sorted.bed -g $GENOME/Dm3.chrom.sizes > ${i}_120-181.bedGraph
  echo bedGraph to bigWig
  $bgToBigWig ${i}_120-181.bedGraph $GENOME/Dm3.chrom.sizes ${i}_120-181.bw
done
STOPPED HERE

_____________________________

##convert the reads to map just the center of the MNase-seq read (important for running edgeR)
for i in ${data[*]}
do
  echo ${i}
## calculate the center of the read and make a new bed file
  echo make bed file with center points
  awk '{AVE=int(($2+$3)/2+0.5); $2=AVE; $3=AVE;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ${i}_120-181_sorted.bed >${i}_120-181_sorted_centers.bed
  echo bed to bedGraph
  genomeCoverageBed -bg -i ${i}_120-181_sorted_centers.bed -g $GENOME/Dm3.chrom.sizes > ${i}_120-181_sorted_centers.bedGraph
  echo bedGraph to bigWig
  $bgToBigWig ${i}_120-181_sorted_centers.bedGraph $GENOME/Dm3.chrom.sizes ${i}_120-181_sorted_centers.bw
done

##make files for the combined replicates
##common part of the replicate file names 

for i in ${dataCombined[*]}
do
  echo ${i}
  echo Combine replicates	
  cat ${i}1_120-181_sorted.bed ${i}2_120-181_sorted.bed | sort -k1,1 -k2,2n > ${i}_120-181.combined.bed
  wc -l ${i}_120-181.combined.bed
  echo bed to bedGraph
  genomeCoverageBed -bg -i ${i}_120-181.combined.bed -g $GENOME/Dm3.chrom.sizes > ${i}_120-181.combined.bedGraph
  echo bedGraph to bigWig
  $bgToBigWig ${i}_120-181.combined.bedGraph $GENOME/Dm3.chrom.sizes ${i}_120-181.combined.bw
done

##to make normalized bedGraph files (in reads per million total reads), we need total combined reads from each treatment
##use wc -l (see above) to count the number of lines in the combined.bed files
##NEED TO FIND A WAY TO AUTOMATE THIS!!!!!!!!!!!
##BF: 72.365747 e6 reads
##PB: 66.339014 e6 reads
##LZ: 62.516925 e6 reads
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4/72.365747}' BF_120-181.combined.bedGraph >BF_120-181.combined_normalized.bedGraph
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4/66.339014}' PB_120-181.combined.bedGraph >PB_120-181.combined_normalized.bedGraph
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4/62.516925}' LZ_120-181.combined.bedGraph >LZ_120-181.combined_normalized.bedGraph

for i in ${dataCombined[*]}
do
  echo ${i}
  $bgToBigWig ${i}_120-181.combined_normalized.bedGraph $GENOME/Dm3.chrom.sizes ${i}_120-181.combined_normalized.bw
done

##adding headers to make the bedgraphs usable in UCSC browser
sed -i '1s/^/track type=bedGraph name="BF_120-180" description="combined_duplicates_S2_size_selected_200-400bp_MNase-seq_120-180bp_reads" visibility=2 graphType=bar color=0,0,200\n/' BF_120-181.combined_normalized.bedGraph

sed -i '1s/^/track type=bedGraph name="PB_120-180" description="combined_duplicates_S2_size_selected_200-400bp_MNase-seq_120-180bp_reads" visibility=2 graphType=bar color=255,75,75\n/' PB_120-181.combined_normalized.bedGraph

sed -i '1s/^/track type=bedGraph name="LZ_120-180" description="combined_duplicates_S2_size_selected_200-400bp_MNase-seq_120-180bp_reads" visibility=2 graphType=bar color=200,0,0\n/' LZ_120-181.combined_normalized.bedGraph

##clean up the computer space by removing non-essential files!
#rm *.bt2
#gzip *.fastq
#gzip *.bedGraph
#gzip *.bed
#gzip *.bw
gzip *.bam
gzip *.sam
####remember to save essential files () to the back-up server!!! 

