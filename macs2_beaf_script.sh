#!/bin/bash
##to run this script all I did was copy and paste it into the terminal

## Set some enviroment variables. (Actually not used; 
##if variable is used, path starts with $ ; if variable is not used, delete $)
##genome direction is where the bowtie2 genome indexes are (not needed, bam already has genome associated?)
##GENOME=/home/chart/genomes/bowtie2_refGenome/dm3

cd /home/chart/scripts/macs2/MACS2
macs2 callpeak -t '/home/chart/modENCODE_bigWigs/oli_bf_chipSeq_bedgraph/SRR1042411_Align.bam' --outdir '/home/chart/modENCODE_bigWigs/oli_bf_chipSeq_bedgraph' -n bf_SRR1042411_2016sep27 -B --call-summits




