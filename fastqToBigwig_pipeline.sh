#run the pipeline on all of the combin

#================================================================================
#PRO-seq shell script
#Ver. 1.1
#Created by Kevin Lee(hwonlee5@gmail.com), Leighton Core(leightoncore@gmail.com)
#Revised by Ravi Patel(ravipatel4@gmail.com)	     
#================================================================================
#===========================================================================================================================================
#Ingredients					      
#1. Your FASTA, rDNA (ebwt), genome (ebwt) *ebwt can be generated by "bowtie-build infile outfile"*  				      
#2. One python script: passfilter.py
#3. fastx_tools					      
#4. bowtie					      
#5. Samtools
#6. BEDTools				              
#===========================================================================================================================================
#Usage
#1. Input file path of necessary files and scripts
#2. Run the script by "sh SCRIPTPATH"
#(optional) if no rDNA, 
#hashtag or remove all the lines in between "++++++++++++++" 
#AND hashtag or remove one "done" in the end of the script
#===========================================================================================================================================
#Input

#try first without the whole path
fastq='LZ1
LZ2
LacZ_12
BF1
BF2
BEAF_12
PB1
PB2
PBro_12'

rDNA="/home/macproadmin/users/GROseq_Utils/refGenomes/dm3/rRNA_DM/rRNA"

genome="/home/macproadmin/users/GROseq_Utils/refGenomes/dm3/genome_dm3/dm_5.22"

passfilter="/home/macproadmin/users/GROseq_Utils/Processing_Scripts/PROSEQ/pipeline/passfilter.py"

chrominfo="/home/macproadmin/users/Fabiana/Dm3.chrom.sizes"

#chrominfo format example = chr1 52949123 (with no title)

bgToBigWig="/home/macproadmin/users/GROseq_Utils/mkbigWig/bedGraphToBigWig"

phix="/home/macproadmin/users/GROseq_Utils/refGenomes/phiX/phiX"

nmer="26"
#===========================================================================================================================================


for file in ${fastq}
do
for pf in ${passfilter}
do
for nmer in ${nmer}
do

total=$(wc -l ${file}.fastq | cut -f1 -d" ")

#echo -e "\n\n####################\n"${file} passfiltering...

#Takes the fastq file and makes _q.tmp file with reads that passed the quality filter which is all the "N" and discarding all the "Y" quality reads.
#python ${pf} ${file}.fastq _q.tmp
#wc -l _q.tmp | awk '{total="'$total'"; printf "Number of reads after passfilter: %s (%s%s)\n", $1/4, $1/total*100, "%"}' >> ${file}.runStat

echo ${file} clipping...

#Takes _q.temp file and clips any residual 3'adapter sequence which could be all or part of TGGAATTCTCGGGTGCCAAGG and generates another temporary file called _c.tmp
fastx_clipper -i ${file}.fastq -o _c.tmp -a TGGAATTCTCGGGTGCCAAGG -l 15 -Q33  
wc -l _c.tmp | awk '{total="'$total'"; printf "Number of reads after clipping: %s (%s%s)\n", $1/4, $1/total*100, "%"}' >> ${file}.runStat

echo ${file} trimming...

#Takes _c.tmp file and trims sequences to set number of bases.
fastx_trimmer -l ${nmer} -i _c.tmp -o _c1.tmp -Q33

echo ${file} reverse complementing...

#Takes the _c1.tmp file and flips the read because the illumina sequences the DNA from 5'end so the adapters in PRO-seq is reversed.
fastx_reverse_complement -i _c1.tmp -o _3.tmp -Q33

#for phix in ${phi}
#do
echo ${phix} aligning to phix...    
echo -e "\n"${phix} aligning to phix...   >> ${file}.runStat

bowtie -p9 -v2 -M1 -q --sam --un ${file}_phixUnalign ${phix} _3.tmp ${file}_phixAlign 2>&1| tee ${file}_phix.log   #p9, number of threads, -v2 allow two mismatches, -M1 for multiple (the 1 means greater than 1 alignment) alignment, randomly assigns one of them (used for ribosomal and PhiX so that they can be discarded later), -q it is in fastq format, --sam output is in sam format, --un outputs the unaligned reads to the specified file (here, that is what does not align to PhiX... we will use that as the input for the next step), ${phix} is the ebwt file from bowtie (the reference genome)
cat ${file}_phix.log >> ${file}.runStat

#+++++++++++++++++++++++++++++++++++++++++++++++++++++HASHTAG IF NO rDNA+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for rib in ${rDNA}
do


echo ${file} aligning to ribosomal DNA...
echo -e "\n"${file} aligning to ribosomal DNA...  >> ${file}.runStat

#Generate 3 files: ribUnalign (not aligned to rDNA sequences), ribAlign (aligned to rDNA), _rib.log (bowtie stats)
bowtie -p9 -v2 -M1 -q --sam --un ${file}_ribUnalign ${rib} ${file}_phixUnalign ${file}_ribAlign 2>&1| tee ${file}_rib.log
cat ${file}_rib.log >> ${file}.runStat
##bowtie -p3 -k2 -m1 -q --sam --un ${file}_ribUnalign ${rib} _3.tmp ${file}_ribAlign 2>&1| tee ${file}_rib.log      
   
#+++++++++++++++++++++++++++++++++++++++++++++++++++++HASHTAG IF NO rDNA+++++++++++++++++++++++++++++++++++++++++++++++++++++++++


for gen in ${genome}
do

echo ${file} aligning to genome...
echo -e "\n"${file} aligning to genome...  >> ${file}.runStat

#Generate 3 files: Unalign (not aligned to genome), Align (aligned to genome), _gen.log (bowtie stats)
bowtie -p9 -v2 -m1 -q --sam --un ${file}_Unalign ${gen} ${file}_ribUnalign ${file}_Align.sam 2>&1| tee ${file}_gen.log   #lowercase -m1 discards reads that have multiple alignments  
cat ${file}_gen.log >> ${file}.runStat

echo ${file} converting SAM to BED...

#Converts SAM to BED
samtools view ${file}_Align.sam -Sb | bamToBed -i stdin > ${file}.bed

for chr in ${chrominfo}
do

#Converts BED to Sorted Bed using command line arguments
echo Sorting BED...
sort -k1,1 -k2,2n ${file}.bed > ${file}_sorted.bed

#add 'chr' to chromosome names
echo adding chr to chrom names
cat ${file}_sorted.bed | awk '{printf "%s%s\t%s\t%s\t%s\t%s\t%s\n", "chr",$1, $2, $3, $4, $5, $6}' > ${file}_sorted_chr.bed

#Generate non-nomalized Bedgraphs 
echo Generating non-normalized Bedgraphs...
awk '$6 == "+"' ${file}_sorted_chr.bed | genomeCoverageBed -i stdin -3 -bg -g ${chr} > ${file}_plus.bedgraph    #-3 uses only the 3' position 
awk '$6 == "-"' ${file}_sorted_chr.bed | genomeCoverageBed -i stdin -3 -bg -g ${chr} > ${file}_m.bedgraph       #name m because this is temporary before making the file where the read depth is multiplied by -1
awk '{$4=$4*-1; print}' ${file}_m.bedgraph > ${file}_minus.bedgraph

#move the mitoremover thing here, so that I can make bigwigs for both later
python MitoRemoverIndividual.py ${file}_plus.bedgraph     #make a file with no mitochondrial reads before making the bigWig 
python MitoRemoverIndividual.py ${file}_minus.bedgraph


#Generate normalized Bedgraphs
a=$(awk '{ sum += $4 } END { print sum }' ${file}_plus_noMito.bedgraph )
b=$(awk '{ sum += $4 } END { print sum }' ${file}_minus_noMito.bedgraph )
d=$(($b*-1))
c=$(expr $a + $d)

echo track type=bedGraph name=${file}.norm.pl description=plus visibility=full color=190,0,30 autoScale=on priority=1 > ${file}_normed_plus_noMito.bedgraph    #header needed for the genome browser
echo track type=bedGraph name=${file}.norm.mn description=minus visibility=full color=40,10,180 autoScale=on priority=2 > ${file}_normed_minus_noMito.bedgraph
echo $c | awk '{ c="'$c'"; printf "%s\t%s\t%s\t%s\n", $1, $2, $3, ($4*1000000)/c}' ${file}_plus_noMito.bedgraph  >> ${file}_normed_plus_noMito.bedgraph       #makes the reads per million mapped calculation
echo $c | awk '{ c="'$c'"; printf "%s\t%s\t%s\t%s\n", $1, $2, $3, ($4*1000000)/c}' ${file}_minus_noMito.bedgraph >> ${file}_normed_minus_noMito.bedgraph

for bigwig in ${bgToBigWig}
do

#Make BigWig
echo Making BigWigs...          #crashed because I needed to include normed in the name

${bigwig} ${file}_normed_plus_noMito.bedgraph ${chr} ${file}_normed_plus_noMito.bw
${bigwig} ${file}_normed_minus_noMito.bedgraph ${chr} ${file}_normed_minus_noMito.bw

${bigwig} ${file}_plus_noMito.bedgraph ${chr} ${file}_plus_noMito.bw    #make bigwigs for the normalized and non-normalized data (could need non norm for Pausing Index, DeSeq)
${bigwig} ${file}_minus_noMito.bedgraph ${chr} ${file}_minus_noMito.bw

rm *.tmp
rm ${file}_m.bedgraph
rm *.sam
rm *Align
rm *Unalign*

#hashtag or remove one "done" if no rDNA
#done
done
done
done
done
done
done
done


## check for existence of subdirectories, make if non-existant 
if [ ! -d "./bed" ]; then 
mkdir "./bed"
fi
if [ ! -d "./bedgraph" ]; then 
mkdir "./bedgraph"
fi
if [ ! -d "./logs" ]; then 
mkdir "./logs"
fi
if [ ! -d "./bw" ]; then 
mkdir "./bw"
fi
if [ ! -d "./bw" ]; then 
mkdir "./bw"
fi


## Move aligned files to corresponding directories
mv *.bed ./bed 
mv *.bw ./bw
mv *.bedgraph ./bedgraph 
mv *.log ./logs
#mv *.txt ./logs 
#mv *.bam* ./bam
mv *run* ./logs


