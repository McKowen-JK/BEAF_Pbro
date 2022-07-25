#Edit the script used to see if the maximum pause window moves 
    #to make a gene list of where the maximum 50 bp window is in the untreated cells
    #I will then use that gene list go get counts for DESeq2
    #Jay brought it up in my group meeting that this could lower noise by focusing on the pause region
    #I see this as a bridge to measuring shifts in distribution within pause: 
      #Next, I want to use this genelist to measure changes in the distribution of reads within this window (i.e. fraction in each quarter... Greg's shape stuff)



library(bigWig)

finalset =  read.table('/home/chart/scripts/finalsetMJG_BedFormat.dat', header=TRUE)
basepath = "/home/chart/PROseq2/bigwig/"
map=load.bigWig("/home/chart/genomes/dm3_single_base_mappability.bigWig")

#Greg's function to scan promoter and find the 50bp region with highest number of counts
getCounts.sliding.pr <- function(wig.p, wig.m, map, genes, off1 = -50, off2 = 150, windowsize = 50) {
  N = dim(genes)[1]
  plusStrand = genes[,6] == '+'
  
  result = vector(mode="integer", length=N)
  mapresult = vector(mode="integer", length=N)
  pausing_region_start = vector(mode="integer", length=N)
  pausing_region_end = vector(mode="integer", length=N)
  
  scanlen = off2 - off1 - windowsize
  
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 2])
    end = as.integer(genes[i, 3])
    
    if (plusStrand[i]) {
      wig = wig.p
      
      scanResults = c()
      for (n in c(1:scanlen)){
        
        qStart = start + off1 + n
        qEnd = start + off1 + n + windowsize
        
        data = query.bigWig(wig, chrom, qStart, qEnd)
        scanResults = append(scanResults, abs(sum(data[,3])))
        
      }
      Offset = which.max(scanResults)  #grabs the index of the first instance of the maximum value ( = bp offset during scanning)
      result[i] = max(scanResults)
      pausing_region_start[i] = start + off1 + Offset
      pausing_region_end[i] = start + off1 + Offset + windowsize
      qMapStart = start + off1 + Offset + 25
      qMapEnd  = start + off1 + Offset + windowsize + 25
      mapdata = query.bigWig(map, chrom, qMapStart, qMapEnd)
      mapresult[i] = abs(sum(mapdata[,3]))
    } 
    else {
      scanResults = c()
      
      wig = wig.m
      if (is.null(wig.m))
        wig = wig.p
      
      for (n in c(1: scanlen)){
        qStart = end - off1 - n - windowsize
        qEnd = end - off1 - n
        
        data = query.bigWig(wig, chrom, qStart, qEnd)
        scanResults = append(scanResults, abs(sum(data[,3])))
      }
      Offset = which.max(scanResults)  #grabs the index of the first instance of the maximum value ( = bp offset during scanning)
      result[i] = max(scanResults)
      pausing_region_start[i] = end - off1 - Offset - windowsize
      pausing_region_end[i] = end - off1 - Offset
      qMapStart = end - off1 - Offset - windowsize
      qMapEnd = end - off1 - Offset 
      mapdata = query.bigWig(map, chrom, qMapStart, qMapEnd)
      mapresult[i] = abs(sum(mapdata[,3]))
    }
  }
  return(cbind(result,mapresult,pausing_region_start,pausing_region_end))
}

LacZ_minus = load.bigWig('/home/chart/PROseq2/bigwig/LacZ_12_minus_noMito.bw')
LacZ_plus = load.bigWig('/home/chart/PROseq2/bigwig/LacZ_12_plus_noMito.bw')

LacZ_MaxPause = getCounts.sliding.pr(LacZ_plus, LacZ_minus, map, finalset, -50, 150, 50)


MaxPauseGenelist = cbind(finalset, LacZ_MaxPause[,3:4])


MaxPauseGenelistFinal = MaxPauseGenelist[,c(1, 7:8, 4:6)]

write.csv(MaxPauseGenelistFinal, '/home/chart/scripts/finalsetMJG_BedFormat_LacZ_50bpPause.csv', row.names = FALSE)


