cat rep.txt | while read line
do
#awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$3-$2}' $line.s.bed | awk -F '\t' '$4 >= 120 && $4 <= 180' | awk -F '\t' '{print $1"\t"$2"\t"$3}' > $line.120_180.bed
#awk -F '\t' '{X=25; mid=(int($2)+int($3))/2;printf("%s\t%d\t%d\n",$1,(mid-X<0?0:mid-X),mid+X);}' $line.120_180.bed | awk -F '\t' 'BEGIN{OFS="\t";} $4=(FNR FS $4)' > $line.x.frag50.bed
#bedtools sort -i $line.x.frag50.bed > $line.frag50.bed
#rm $line.x.frag50.bed
#bedToBam -i $line.frag50.bed -g dm3.genome > $line.filter.bam
#samtools index $line.filter.bam
bamCoverage -b $line.filter.bam -bs 1 -p 6 --normalizeUsing RPKM -o $line.f50.RPKM.bw
done


