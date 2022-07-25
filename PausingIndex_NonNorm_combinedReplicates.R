#The goal here is to make a list of all genes that are paused, to use as a filter for my CDF analysis (I'm only interested in paused genes)
  #To test if a gene is paused, we use a Fisher's exact test to see if the density of mapable reads in the promoter is different than that across the gene body
  #I edited Fabiana's script for this.  She use the max pause window from here control to look in the otheres
  #I want to look for pausing in all of my conditions, and believe that the window could shift, so I used the max pause window determined individually for each treatment
  #in the end, I made a gene list that has all genes that are paused (with fe test p value < .01) in all four of my conditions
  #I have the least number of paused genes (by a lot... 1200 fewer than next highest) for my aptamer Cu cells
  #but, I think that this is because it has the lowest read depth (after looking at some genes with marginal p value), so more reads could make this list larger

#Use non-normalized, combined replicate bigwigs (not normalized) because the fisher exact test needs count data so it will not work with normalized data

#From this analysis, I also have pausing information for all genes. When I look at just pausing indices, I see something interesting:
  #There is a global decrease in pausing (slight), in the induced aptamer containing cells
  #this supports what I saw with the CDF of promoter read counts, 
  #but I like this analysis more as it is a ratio and thus less ambiguous than looking simply at a reduciton in normalized read count

library(lattice)
library(grid)
library(bigWig)
#library(latticeExtra)

finalset = read.csv("/home/chart/scripts/gene_lists/finalsetMJG_CG.csv", header=T)
basepath = "/home/chart/PROseq2/bigwig/"
map=load.bigWig("/home/chart/genomes/dm3_single_base_mappability.bigWig")

#Greg's function to scan promoter and find the 50bp region with highest number of countws
getCounts.sliding.pr <- function(wig.p, wig.m, map, genes, off1 = -50, off2 = 150, windowsize = 50) {
  N = dim(genes)[1]
  plusStrand = genes[,2] == '+'
  
  result = vector(mode="integer", length=N)
  mapresult = vector(mode="integer", length=N)
  pausing_region_start = vector(mode="integer", length=N)
  pausing_region_end = vector(mode="integer", length=N)

  scanlen = off2 - off1 - windowsize
  
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 3])
    end = as.integer(genes[i, 4])
    
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
      qMapStart = start + off1 + Offset + 25    #offset by 25 here only because mappability is calculated using 26mers, so on the + strand a bp's mappablility depends on the 25 bp leading up to it
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

getCounts.gb <- function(wig.p, wig.m, genes, off1) {
  N = dim(genes)[1]
  plusStrand = genes[,2] == '+'
  
  result = vector(mode="integer", length=N)
  
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 3])
    end = as.integer(genes[i, 4])

    if (plusStrand[i]) {
      wig = wig.p

      qStart = start + 200
      qEnd = end - 200
            
      if (qStart > qEnd) {
        result[i] = NA
        next
      }

      data = query.bigWig(wig, chrom, qStart, qEnd)

      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    } else {
      wig = wig.m
      if (is.null(wig.m))
        wig = wig.p

      qStart = start + 200
      qEnd = end - 200

      if (qStart > qEnd) {
        result[i] = NA
        next
      }

      data = query.bigWig(wig, chrom, qStart, qEnd)

      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    }
  }

  return(cbind(result))
}

getCounts.gb.mappability <- function(genes, off1) {
  N = dim(genes)[1]
  plusStrand = genes[,2] == '+'
  
  result = vector(mode="integer", length=N)
  
  for (i in 1:N) {
    chrom = as.character(genes[i, 1])
    start = as.integer(genes[i, 3])
    end = as.integer(genes[i, 4])

    if (plusStrand[i]) {

      qStart = start + 225   #offset by 25 here only because mappability is calculated using 26mers, so on the + strand a bp's mappablility depends on the 25 bp leading up to it
      qEnd = end - 175
          
      if (qStart > qEnd) {
        result[i] = NA
        next
      }

      data = query.bigWig(map, chrom, qStart, qEnd)
    
      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    } else {
      
      qStart = start + 200
      qEnd = end - 200

      if (qStart > qEnd) {
        result[i] = NA
        next
      }
       
      data = query.bigWig(map, chrom, qStart, qEnd)

      if (!is.null(data))
        result[i] = abs(sum(data[,3]))
    }
  }

  return(cbind(result))
}

approx.ratio.CI <- function(x1, x2, alpha=0.05) {
  t = qnorm(1 - alpha/2)
  n = x1 + x2
  zp = (t^2/2 + x1 + 1/2)^2 - ((n + t^2)/n) * (x1 + 1/2)^2
  zn = (t^2/2 + x1 - 1/2)^2 - ((n + t^2)/n) * (x1 - 1/2)^2
  
  a = (t^2/2 + x1 + 1/2 + sqrt(zp)) / (n + t^2/2 - x1 - 1/2 - sqrt(zp))
  b = (t^2/2 + x1 - 1/2 - sqrt(zn)) / (n + t^2/2 - x1 + 1/2 + sqrt(zn))

  return(c(b, a))
}

approx.ratios.CI <- function(num.counts, denom.counts, alpha=0.05) {
  stopifnot(length(num.counts) == length(denom.counts))
  N = length(num.counts)

  result = matrix(data=0, nrow=N, ncol=2)

  for (i in 1:N)
    result[i, ] = approx.ratio.CI(num.counts[i], denom.counts[i], alpha)

  return(result)
}

load.wigset <- function(row) {
  file = wig.table[row, 1]
  wig.p = NULL
  if (file != "")
    wig.p = load.bigWig(paste(basepath, file, sep=''))
  file = wig.table[row, 2]
  wig.m = NULL
  if (file != "")
    wig.m = load.bigWig(paste(basepath, file, sep=''))
  return(list(wig.p, wig.m, wig.table[row, 3]))
}

unload.wigset <- function(set) {
  if (!is.null(set[[1]]))
    unload.bigWig(set[[1]])
  if (!is.null(set[[2]]))
    unload.bigWig(set[[2]])
}

pause.index <- function(wig.p, wig.m, genes, name, alpha=0.05) {
  pr.info = getCounts.sliding.pr(wig.p, wig.m, map, genes, -50, 150, 50)
  
  pr.counts = pr.info[,1]
  gb.counts = getCounts.gb(wig.p, wig.m, genes, 200)

  pr.map = pr.info[,2]
  gb.map = getCounts.gb.mappability(genes, 200) 

  pr.density = pr.counts / pr.map
  pr.density[pr.map == 0] = NA
  gb.density = gb.counts[,1] / gb.map[,1]
  gb.density[gb.map[,1] == 0] = NA

  pause.index = pr.counts / gb.counts[,1]
  pause.index[gb.counts[,1] == 0] = NA

  cis = approx.ratios.CI(pr.counts, gb.counts[,1], alpha)

  pause.index.scaled = pause.index * (gb.map[,1]/pr.map)
  pause.index.scaled[gb.map[,1] == 0] = NA
  pause.index.scaled[pr.map == 0] = NA
  ci.low.scaled = cis[,1] * (gb.map[,1]/pr.map)
  ci.low.scaled[gb.map[,1] == 0] = NA
  ci.low.scaled[pr.map == 0] = NA
  ci.high.scaled = cis[,2] * (gb.map[,1]/pr.map)
  ci.high.scaled[gb.map[,1] == 0] = NA
  ci.high.scaled[pr.map == 0] = NA

  result = cbind(pr.info[,3], pr.info[,4], pr.counts, gb.counts[,1], pr.map, gb.map[,1], pr.density, gb.density, pause.index, cis[,1], cis[,2], pause.index.scaled, ci.low.scaled, ci.high.scaled)
  colnames(result) <- c(paste(name,"_pausing_region_start", sep=""), paste(name,"_pausing_region_end", sep=""), paste(name,"_pr_counts", sep=""), paste(name,"_gb_counts", sep=""), paste(name,"_pr_map", sep=""), paste(name,"_gb_map", sep=""), paste(name,"_pr_density", sep=""), paste(name,"_gb_density", sep=""), paste(name,"_pause_index", sep=""), paste(name,"_pause_index_CI_low", sep=""), paste(name,"_pause_index_CI_high", sep=""), paste(name,"_pause_index_scaled", sep=""), paste(name,"_pause_index_CI_low_scaled", sep=""), paste(name,"_pause_index_CI_high_scaled", sep=""))
  
  return(result)
}


################################Find the pausing index for each in a new max pause window.  
  #I want to use this to filter for pausing here, so I think I want to allow it to use the max window for each conditions (could change when it is shifted, though it is still paused)


wig.table = rbind(
  c("LacZ_12_plus_noMito.bw", "LacZ_12_minus_noMito.bw", "LacZ12"),
  c("BEAF_12_plus_noMito.bw", "BEAF_12_minus_noMito.bw", "BEAF12"),
  c("PBro_12_plus_noMito.bw", "PBro_12_minus_noMito.bw", "PBro12"))



N = dim(wig.table)[1]
pi.res = vector(mode="list", length=N)
for (i in 1:N) {
  cat("* loading", i, "\n")
  wigs = load.wigset(i)#this should give all the bigwig files
  cat("* computing ...\n")
  pi.res[[i]] = pause.index(wigs[[1]], wigs[[2]], finalset, wigs[[3]])
  cat("* unloading.\n")
  unload.wigset(wigs)
}


#keep separate dataframes for each in addition to the one big one

pi.output <- cbind(finalset, pi.res[[1]][,1:14],  pi.res[[2]][,1:14],  pi.res[[3]][,1:14])
#write.csv(pi.output, '/media/jacobmtome/HDD/JMT_Lis_Lab/2015_07_PRO_seq/PausingIndex/151028_PausingInformation_all4.csv')
write.csv(pi.output, '/home/chart/PROseq2/pausing/2020oct06_pauseInfo_LzBfPb12.csv')
#Above line makes pi file for all files, I only need control (LacZ)

LacZ12_pi= cbind(finalset, pi.res[[1]][,1:14])
BEAF12_pi= cbind(finalset, pi.res[[2]][,1:14])
PBro12_pi= cbind(finalset, pi.res[[3]][,1:14])



plot(pi.output$LacZ12_pause_index, pi.output$BEAF12_pause_index, log='xy')
abline(a=0, b=1, col='red', lwd=3)

plot(pi.output$LacZ12_pause_index, pi.output$PBro12_pause_index, log='xy')
abline(a=0, b=1, col='red', lwd=3)


wilcox.test(pi.output$LacZ12_pause_index, pi.output$BEAF12_pause_index, paired = TRUE)
#Wilcoxon signed rank test with continuity correction

#data:  pi.output$LacZ_pause_index and pi.output$BEAF_pause_index
#V = 20915000, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(pi.output$LacZ12_pause_index, pi.output$PBro12_pause_index, paired = TRUE)
#Wilcoxon signed rank test with continuity correction

#data:  pi.output$LacZ_pause_index and pi.output$PBro_pause_index
#V = 22366000, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0




##adapted from Leighton

###The following tests for pausing... paused genes pass the fischer exact test
  #test if the read density in the promoter is significantly different than the read density across the whole gene (promoter plus gene body)
  #Fabs just did this for LacZ... i want to look at genes that are paused in both (or all conditions)
  #So, make this into a loop that makes a new genelist of paused genes only for each condition
  #because R loops are bitches, use lapply: http://stackoverflow.com/questions/14429686/change-multiple-dataframes-in-a-loop

fe.test= function(data, prReads, prmappL, gbReads, gbmappL, colName='')
{
  newtable = data.frame(data)
  rows = nrow(data)
  for(ii in 1:nrow(data))
  {
    #print(data[ii,])
    fmat = matrix(c((data[ii, prReads]),(data[ii, prmappL]),(data[ii, gbReads]),(data[ii, gbmappL])), nrow=2) 
    #print(fmat)
    f = fisher.test(fmat,alternative='greater')#it has to be a one-tailed fisher-exact test, or it will give you the de-enriched genes as well
    newtable[ii, colName] = p.adjust(f$p.value, method = 'bonferroni', n = rows)
  }
  return(newtable)
}

MyDataSets = list(LacZ_pi, BEAF_pi, PBro_pi)

MyDataSets <- lapply(MyDataSets,function(df){
  df$EP<-as.integer((df[,15]+df[,16])*((df[,17]/(df[,17]+df[,18]))))
  
  df<-df[!is.na(df$EP),]
  
  df$EG<-as.integer(df[,15]+df[,16]-df$EP)
  
  df<-df[!is.na(df$EG),]
  df <- fe.test(df, 15, 27, 16, 28, colName="pVal_Alt")
  df
  })

LacZ_pi_ratios = MyDataSets[[1]]
BEAF_pi_ratios = MyDataSets[[2]]
PBro_pi_ratios = MyDataSets[[3]]


LacZ_pi_ratios_sig = LacZ_pi_ratios[LacZ_pi_ratios$pVal_Alt <= 0.05,]
BEAF_pi_ratios_sig = BEAF_pi_ratios[BEAF_pi_ratios$pVal_Alt <= 0.05,]
PBro_pi_ratios_sig = PBro_pi_ratios[PBro_pi_ratios$pVal_Alt <= 0.05,]

#I could not find a more elegant way to merge more than one dataframe
#PausedGenesForFilter = merge(merge(merge(LacZ_pi_ratios_sig, BEAF_pi_ratios_sig, by='name'), PBro_pi_ratios_sig, by='name'), by='name')

LacZ_pi_ratios_GeneList = LacZ_pi_ratios_sig[,c(1, 3:4, 13:14, 2, 7:9, 12, 15:29)]
write.csv(LacZ_pi_ratios_GeneList, '/home/chart/PROseq2/pausing/2020aug14_LacZ_PauseIndex_GeneList.csv', row.names = FALSE)
LacZ_pi_all_GeneList = LacZ_pi_ratios[,c(1, 3:4, 13:14, 2, 7:9, 12, 15:29)]
write.csv(LacZ_pi_all_GeneList, '/home/chart/PROseq2/pausing/2020aug14_LacZ_PauseIndex_all_GeneList.csv', row.names = FALSE)

BEAF_pi_ratios_GeneList = BEAF_pi_ratios_sig[,c(1, 3:4, 13:14, 2, 7:9, 12, 15:29)]
write.csv(BEAF_pi_ratios_GeneList, '/home/chart/PROseq2/pausing/2020aug14_BEAF_PauseIndex_GeneList.csv', row.names = FALSE)
BEAF_pi_all_GeneList = BEAF_pi_ratios[,c(1, 3:4, 13:14, 2, 7:9, 12, 15:29)]
write.csv(BEAF_pi_all_GeneList, '/home/chart/PROseq2/pausing/2020aug14_BEAF_PauseIndex_all_GeneList.csv', row.names = FALSE)

PBro_pi_ratios_GeneList = PBro_pi_ratios_sig[,c(1, 3:4, 13:14, 2, 7:9, 12, 15:29)]
write.csv(PBro_pi_ratios_GeneList, '/home/chart/PROseq2/pausing/2020aug14_PBro_PauseIndex_GeneList.csv', row.names = FALSE)
PBro_pi_all_GeneList = PBro_pi_ratios[,c(1, 3:4, 13:14, 2, 7:9, 12, 15:29)]
write.csv(PBro_pi_all_GeneList, '/home/chart/PROseq2/pausing/2020aug14_PBro_PauseIndex_all_GeneList.csv', row.names = FALSE)

#save.image("/media/jacobmtome/HDD/JMT_Lis_Lab/2015_07_PRO_seq/PausingIndex/CalculatePausingIndex_NewWindow_Each_NonNormalized.RData")




###########################################################Plots to examine changes in pausing index

LacZ_color  ="#B22222" #firebrick red
BEAF_color ="#483D8B" #blue hex color
PBro_color="#009933" #green

#make vectors of the pausing indices of the treatments


#make vectors that have all genes with pausing information
LacZ_PauseIndices_all = log10(LacZ_pi_ratios$LacZ_pause_index[!is.na(LacZ_pi_ratios$LacZ_pause_index) & !is.infinite(log10(LacZ_pi_ratios$LacZ_pause_index))])
BEAF_PauseIndices_all = log10(BEAF_pi_ratios$BEAF_pause_index[!is.na(BEAF_pi_ratios$BEAF_pause_index) & !is.infinite(log10(BEAF_pi_ratios$BEAF_pause_index))])
PBro_PauseIndices_all = log10(PBro_pi_ratios$PBro_pause_index[!is.na(PBro_pi_ratios$PBro_pause_index) & !is.infinite(log10(PBro_pi_ratios$PBro_pause_index))])


plot(density(PBro_PauseIndices_all))
lines(density(LacZ_PauseIndices_all), col='red')

plot(density(BEAF_PauseIndices_all))
lines(density(LacZ_PauseIndices_all), col='red')

plot(density(PBro_PauseIndices_all))
lines(density(BEAF_PauseIndices_all), col='purple')
lines(density(LacZ_PauseIndices_all), col='red')

#now show changes in pausing index with vioplot
library(vioplot)

#Plot the normalized reads... I think I need to filter out the lowly expressed guys

#png('151027_NormMeanExpression_MaxPauseWindow.png')
plot(1,1,xlim=c(0.5,3.5),ylim=range(c(LacZ_PauseIndices_all,BEAF_PauseIndices_all,PBro_PauseIndices_all)),type="n", xlab="",ylab="log(DESeq2 Norm. Mean Expression)",axes=FALSE, cex.lab=1.25)
## bottom axis, with user-specified labels
title(main='Normalized Mean Expression \nMeasured by DESeq2', cex=1.25)
axis(side=1,at=1:3,labels=c("LacZ", "BEAF", "PBro"), cex.axis=1.25)
axis(side=2, cex.axis=1.25)
vioplot(LacZ_PauseIndices_all,at=1,col=LacZ_color,add=TRUE)
vioplot(BEAF_PauseIndices_all,at=2,col=BEAF_color,add=TRUE)
vioplot(PBro_PauseIndices_all,at=3,col=PBro_color,add=TRUE)
grid()
abline(h=0, lty=3)
#dev.off()


################################If I want to do statistics, paired tests are best (wilcoxon), so do this with only the set of genes with good pausing index information for all (just has to be not na of inf... so no zero reads... I could maybe make this better with a minimum read count)
GoodData = !is.na(pi.output$LacZ_pause_index) & !is.infinite(log10(pi.output$LacZ_pause_index)) & !is.na(pi.output$BEAF_pause_index) & !is.infinite(log10(pi.output$BEAF_pause_index)) & !is.na(pi.output$PBro_pause_index) & !is.infinite(log10(pi.output$PBro_pause_index))

Good_pi.output = pi.output[GoodData,]


#####Test if the pausing indices change between samples.  With this more filtered data set, I now see a change in the wt cells with copper added
wilcox.test(Good_pi.output$PBro_pause_index, Good_pi.output$LacZ_pause_index, paired = TRUE)
#Wilcoxon signed rank test with continuity correction

#data:  Good_pi.output$PBro_pause_index and Good_pi.output$LacZ_pause_index
#V = 1867090, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(Good_pi.output$BEAF_pause_index, Good_pi.output$LacZ_pause_index, paired = TRUE)
#Wilcoxon signed rank test with continuity correction

#data:  Good_pi.output$BEAF_pause_index and Good_pi.output$LacZ_pause_index
#V = 3080720, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0


plot(density(log10(Good_pi.output$PBro_pause_index)), lwd=4, col = PBro_color, main='log10(Pausing Index) for All Genes\nwith Pausing Indices', cex.main=1.25, cex.lab=1.25, cex.axis=1.25, xlab='log10(Pausing Index)')
lines(density(log10(Good_pi.output$LacZ_pause_index)), col=LacZ_color, lwd=4)

lines(density(log10(Good_pi.output$BEAF_pause_index)), col=BEAF_color, lwd=4)
legend("topleft", c("LacZ", "BEAF", "PBro"), cex=1.125, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(BEAF_color, LacZ_color, PBro_color), y.intersp = 1.5)
text(2.2, .4, paste('N =', dim(Good_pi.output)[1]), cex = 1.25)





plot(ecdf(log10(Good_pi.output$BEAF_pause_index)), xlim=c(-2.5,2.25), col=BEAF_color, lwd = 4, xlab='log10(Pausing Index)', ylab='Cumulative Density', cex.axis=1.25, cex.lab=1.25, main='log10(Pausing Index) for All Genes \nWith > 0 Counts', cex.main=1.25)
lines(ecdf(log10(Good_pi.output$LacZ_pause_index)), col=LacZ_color, lwd = 4)
lines(ecdf(log10(Good_pi.output$PBro_pause_index)), col=PBro_color, lwd = 4, lty = 3)
legend("topleft", c("LacZ", "BEAF", "PBro"), cex=1.125, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(BEAF_color, LacZ_color, PBro_color), y.intersp = 1.5)
grid()
abline(v=mean(log10(Good_pi.output$PBro_pause_index)), col=PBro_color, lty=3, lwd=3)
abline(v=mean(log10(Good_pi.output$LacZ_pause_index)), col=LacZ_color, lty=3, lwd=3)
text(1.5, .1, paste('N =', length(Good_pi.output$PBro_pause_index)), cex = 1.25)



#now show changes in pausing index with vioplot

#Plot the normalized reads... I think I need to filter out the lowly expressed guys

meanlog10PI = mean(c(log10(Good_pi.output$PBro_pause_index),log10(Good_pi.output$LacZ_pause_index),log10(Good_pi.output$BEAF_pause_index)))

plot(1,1,xlim=c(0.5,3.5),ylim=range(c(log10(Good_pi.output$PBro_pause_index),log10(Good_pi.output$LacZ_pause_index),log10(Good_pi.output$BEAF_pause_index))),type="n", xlab="",ylab="log10(Pausing Index)",axes=FALSE, cex.lab=1.25)
## bottom axis, with user-specified labels
title(main='log10(Pausing Index) for Genes \nwith > 0 counts', cex=1.25)
axis(side=1,at=1:3,labels=c( "LacZ", "BEAF", "PBro"), cex.axis=1.25)
axis(side=2, cex.axis=1.25)
vioplot(log10(Good_pi.output$LacZ_pause_index),at=1,col=LacZ_color,add=TRUE)
vioplot(log10(Good_pi.output$BEAF_pause_index),at=2,col=BEAF_color,add=TRUE)
vioplot(log10(Good_pi.output$PBro_pause_index),at=3,col=PBro_color,add=TRUE)
abline(h=meanlog10PI, lty=3)
text(4, 4, paste('N =', length(Good_pi.output$PBro_pause_index)), cex = 1.25)
#grid()




################################If I want to do statistics, paired tests are best (wilcoxon), so do this with only the set of genes with good pausing index information for all (just has to be not na of inf... so no zero reads... I could maybe make this better with a minimum read count)
#the first filter I tried only removes genes that have no reads in either the gene body of promoter
#make a more stringent one that looks forat least 25 reads in each

minreadcount = 24
GoodData = pi.output$LacZ_pr_counts > minreadcount & pi.output$LacZ_gb_counts > minreadcount & pi.output$BEAF_pr_counts > minreadcount & pi.output$BEAF_gb_counts > minreadcount & pi.output$PBro_pr_counts > minreadcount & pi.output$PBro_gb_counts > minreadcount

Good_pi.output = pi.output[GoodData,]


#####Test if the pausing indices change between samples.  With this more filtered data set, I now see a change in the wt cells with copper added
wilcox.test(Good_pi.output$PBro_pause_index, Good_pi.output$LacZ_pause_index, paired = TRUE)
#Wilcoxon signed rank test with continuity correction

#data:  Good_pi.output$PBro_pause_index and Good_pi.output$LacZ_pause_index
#V = 411731, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(Good_pi.output$BEAF_pause_index, Good_pi.output$LacZ_pause_index, paired = TRUE)
#Wilcoxon signed rank test with continuity correction

#data:  Good_pi.output$BEAF_pause_index and Good_pi.output$LacZ_pause_index
#V = 715624, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0


plot(density(log10(Good_pi.output$PBro_pause_index)), lwd=4, col = PBro_color, main='log10(Pausing Index) for All Genes\nwith > 25 Counts', cex.main=1.25, cex.lab=1.5, cex.axis=1.5, xlab='log10(Pausing Index)')
lines(density(log10(Good_pi.output$LacZ_pause_index)), col=LacZ_color, lwd=4)

lines(density(log10(Good_pi.output$BEAF_pause_index)), col=BEAF_color, lwd=4)
legend("topleft", c("LacZ", "BEAF", "PBro"), cex=1, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(BEAF_color, LacZ_color, PBro_color), y.intersp = 1.5)
text(1.8, .4, paste('N =', dim(Good_pi.output)[1]), cex = 1.25)



plot(ecdf(log10(Good_pi.output$BEAF_pause_index)), xlim=c(-2.5,2.25), col=BEAF_color, lwd = 4, xlab='log10(Pausing Index)', ylab='Cumulative Density', cex.axis=1.25, cex.lab=1.25, main='log10(Pausing Index) for All Genes \nWith > 25 Counts', cex.main=1.25)
lines(ecdf(log10(Good_pi.output$LacZ_pause_index)), col=LacZ_color, lwd = 4)
lines(ecdf(log10(Good_pi.output$PBro_pause_index)), col=PBro_color, lwd = 4, lty = 3)
legend("topleft", c("LacZ", "BEAF", "PBro"), cex=1.125, bty= "n", ncol=1, inset =0.01, lwd=3, col = c(LacZ_color, BEAF_color, PBro_color), y.intersp = 1.5)
grid()
abline(v=mean(log10(Good_pi.output$PBro_pause_index)), col=PBro_color, lty=3, lwd=3)
abline(v=mean(log10(Good_pi.output$LacZ_pause_index)), col=LacZ_color, lty=3, lwd=3)
text(1.5, .1, paste('N =', length(Good_pi.output$PBro_pause_index)), cex = 1.25)



#now show changes in pausing index with vioplot

#Plot the normalized reads... I think I need to filter out the lowly expressed guys

meanlog10PI = mean(c(log10(Good_pi.output$PBro_pause_index),log10(Good_pi.output$LacZ_pause_index),log10(Good_pi.output$BEAF_pause_index)))

plot(1,1,xlim=c(0.5,3.5),ylim=range(c(log10(Good_pi.output$PBro_pause_index),log10(Good_pi.output$LacZ_pause_index),log10(Good_pi.output$BEAF_pause_index))),type="n", xlab="",ylab="log10(Pausing Index)",axes=FALSE, cex.lab=1.25)
## bottom axis, with user-specified labels
title(main='log10(Pausing Index) for Genes \nwith > 25 counts', cex=1.25)
axis(side=1,at=1:3,labels=c( "LacZ", "BEAF", "PBro"), cex.axis=1.25)
axis(side=2, cex.axis=1.25)
vioplot(log10(Good_pi.output$LacZ_pause_index),at=1,col=LacZ_color,add=TRUE)
vioplot(log10(Good_pi.output$BEAF_pause_index),at=2,col=BEAF_color,add=TRUE)
vioplot(log10(Good_pi.output$PBro_pause_index),at=3,col=PBro_color,add=TRUE)
abline(h=meanlog10PI, lty=3)
text(4, 3.2, paste('N =', length(Good_pi.output$PBro_pause_index)), cex = 1.25)
#grid()


plot(log10(Good_pi.output$PBro_pause_index), log10(Good_pi.output$LacZ_pause_index), xlab = 'log10(Pausing Index) in PBro RNAi', ylab = 'log10(PausIn) in LacZ RNAi cells', cex.lab=1.25, cex.axis=1.25, xlim=c(-3, 3), ylim=c(-3, 3))
abline(a=0, b=1, col='red', lwd=2, lty=3)
title(main='Pausing Index in LacZ and PBro \nRNAi-treated S2 cells', cex=1.25)
text(2, -2.8, paste('N =', length(Good_pi.output$PBro_pause_index)), cex = 1.25)


plot(log10(Good_pi.output$PBro_pause_index), log10(Good_pi.output$BEAF_pause_index), xlab = 'log10(Pausing Index) PBro RNAi', ylab = 'log10(PausIn) BEAF RNAi', cex.lab=1.25, cex.axis=1.25, xlim=c(-3, 3), ylim=c(-3, 3))
abline(a=0, b=1, col='red', lwd=2, lty=3)
title(main='Pausing Index in BEAF and PBro \nRNAi-treated S2 Cells', cex=1.25)
text(2, -2.8, paste('N =', length(Good_pi.output$PBro_pause_index)), cex = 1.25)


plot(log10(Good_pi.output$BEAF_pause_index), log10(Good_pi.output$LacZ_pause_index), xlab = 'log10(Pausing Index) BEAF RNAi', ylab = 'log10(PausIn) LacZ RNAi', cex.lab=1.25, cex.axis=1.25, xlim=c(-3, 3), ylim=c(-3, 3))
abline(a=0, b=1, col='red', lwd=2, lty=3)
title(main='Pausing Index in LacZ and BEAF \nRNAi treated S2 cells', cex=1.25)
text(2, -2.8, paste('N =', length(Good_pi.output$PBro_pause_index)), cex = 1.25)

