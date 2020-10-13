


get.DFS <- function(base_counts, site_counts=NULL, Ns=NULL){
  
  #if the number of haplotypes per population is not specified, assume it is the maximum value observed - NOT RECOMMENDED
  if (is.null(Ns) == TRUE) Ns <- apply(base_counts, 2, max)
  
  #site counts are used to compress the input data. They give the number of sites observed with those base counts
  #if not specified, assume all patterns are represented once
  if (is.null(site_counts) == TRUE) site_counts <- rep(1, nrow(base_counts))
  
  if (!(((ncol(base_counts) == 4 & length(Ns) ==4) | (ncol(base_counts) == 3 & length(Ns) == 3)) & Ns[1] == Ns[2])){
    print("Incorrect specification")
    return()
  }
  
  #convert base_counts to frequencies
  freqs = base_counts / t(replicate(nrow(base_counts), Ns))
  
  N = Ns[1]
  
  #identify sites where P1 and P2 each have a specified allele frequency
  idx_1 = lapply(1:N, function(i) which(base_counts[,1] == i))
  idx_2 = lapply(1:N, function(i) which(base_counts[,2] == i))
  
  #get total counts of each pattern
  pattern_sums_by_count_1 <- sapply(1:N, function(i) apply(get.patterns(freqs[idx_1[[i]],1],
                                                                        freqs[idx_1[[i]],2],
                                                                        freqs[idx_1[[i]],3],
                                                                        freqs[idx_1[[i]],4]) * site_counts[idx_1[[i]]],2,sum))
  
  pattern_sums_by_count_2 <- sapply(1:N, function(i) apply(get.patterns(freqs[idx_2[[i]],1],
                                                                        freqs[idx_2[[i]],2],
                                                                        freqs[idx_2[[i]],3],
                                                                        freqs[idx_2[[i]],4]) * site_counts[idx_2[[i]]],2,sum))
  
  ABBA_by_count <- pattern_sums_by_count_2["ABBA",]
  BABA_by_count <- pattern_sums_by_count_1["BABA",]
  
  DFS <- (ABBA_by_count - BABA_by_count) /
    (ABBA_by_count + BABA_by_count)
  
  weights <-     (ABBA_by_count + BABA_by_count)/
    sum((ABBA_by_count + BABA_by_count))
  
  data.frame(DFS=DFS, weights=weights, ABBA=ABBA_by_count, BABA=BABA_by_count)
}


#function to compute the overall D statistic
get.D.from.base.counts <- function(base_counts, site_counts=NULL, Ns=NULL, full=FALSE){
  
  #if the number of haplotypes per population is not specified, assume it is the maximum value observed - NOT RECOMMENDED
  if (is.null(Ns) == TRUE) Ns <- apply(base_counts, 2, max)
  
  #site counts are used to compress the input data. They give the number of sites observed with those base counts
  #if not specified, assume all patterns are represented once
  if (is.null(site_counts) == TRUE) site_counts <- rep(1, nrow(base_counts))
  
  if (!((ncol(base_counts) == 4 & length(Ns) ==4) | (ncol(base_counts) == 3 & length(Ns) == 3)) ){
    print("Incorrect specification")
    return()
  }
  
  #convert base_counts to frequencies
  freqs = base_counts / t(replicate(nrow(base_counts), Ns))
  
  idx <- 1:nrow(freqs)
  
  #get total counts of each pattern
  pattern_sums <- apply(get.patterns(freqs[idx,1],freqs[idx,2],freqs[idx,3],freqs[idx,4]) * site_counts,2,sum)
  
  if (full == FALSE){
    ABBA <- pattern_sums["ABBA"]
    BABA <- pattern_sums["BABA"]
    return(as.numeric((ABBA - BABA) / (ABBA + BABA)))
  }
  else{
    ABBA_BAAB <- pattern_sums["ABBA_BAAB"]
    BABA_ABAB <- pattern_sums["BABA_ABAB"]
    return(as.numeric((ABBA_BAAB - BABA_ABAB) / (ABBA_BAAB + BABA_ABAB)))
  }
}


plotDFS <- function(DFS, weights, method="lines", ylim=c(-1,1), show_D=TRUE,
                    col="black", col_D="black", width_scale=100, no_xlab=FALSE, add=FALSE){
  
  if (method == "lines"){
    N = length(DFS)
    if (add == FALSE){
      plot(0, xlim = c(1,N), ylim = ylim, cex=0, xlab = "", ylab = "", xaxt="n", bty="n")
      abline(h=0)
    }
    segments(1:N, 0, 1:N, DFS, lwd = width_scale*weights, lend=1, col=col)
  }
  
  if (method == "bars") barplot(DFS, col= rgb(0,0,0,weights), ylim = ylim, add=add)
  
  if (method == "scaled_bars") barplot(DFS*weights, ylim = ylim, add=add)
  
  if (no_xlab == FALSE & add == FALSE) mtext(1,text="Derived allele frequency", line = 0)
  
  if (add == FALSE) mtext(2,text=expression(italic("D")), line = 2.8, las=2)
  if (show_D == TRUE) abline(h= sum(DFS * weights), lty = 2, col=col_D)
  
}


plot.dcfs <- function(dcfs){
  plot(dcfs, type="b")
}



#Example code to run the function
#want to calculate Dfs with P1=cali, P2=island, P3=woodhouse, P4=FL
#also with P1=woodhouse, P2=texana, P3=sumichrasti, P4=FL
#also with P1=sumichrasti, P2=remota, P3=woodhouse, P4=FL
#also with P1=obscura+californica, P2=immanis, P3=island, P4=FL
#compare the distributions to compare the timing/extent of introgression between woodhouse/cali and woodhouse/sumi

### import the frequency spectrum
FS <- read.table("empirical_data/Arabidopsis/arn_lyr_72.DP5MIN58MAC2.lyrata2_lyrata4_arenosa4.sfs")

### get Dfs

dfs_data <- get.DFS(base_counts=FS[,-4], #base counts are the first three columns (i.e everything minus column 4)
                    site_counts=FS[,4], # site counts are column 4
                    Ns = c(14,14,4)) # Ns provide the haploid sample sizes of each population (1 and 2 must always be equal) 

### plot

png("images/DFS_arabidopsis.lyr2_lyr4_arn4.png", width = 2000, height = 1000, res=300)
par(mar=c(1,4,1,1))
plotDFS(dfs_data$DFS, dfs_data$weights, method="lines", col_D="red", no_xlab=T)
dev.off()
