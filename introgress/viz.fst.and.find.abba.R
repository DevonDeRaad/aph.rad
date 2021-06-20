#search for ABBA loci

library(vcfR)
vcf<-read.vcfR("/Users/devder/Desktop/aph.data/unzipped.filtered.vcf")
#read in locality info for samples
locs<-read.csv("~/Desktop/aph.data/rad.sampling.csv")
#subset locs to include only samples that passed filtering, and have been retained in the vcf
locs<-locs[locs$id %in% colnames(vcf@gt),]
rownames(locs)<-1:95

#remove sumi, non-us wood
locs<-locs[c(1:46,57:73),]
vcf@gt<-vcf@gt[,c(1:47,58:74)]
colnames(vcf@gt)[-1] == locs$id
rownames(locs)<-1:nrow(locs)
dim(vcf)
#remove SNPs with no data
source("~/Desktop/snpfiltR/min_mac.R")
vcf<-min_mac(vcf, min.mac = 1)
dim(vcf)

#convert to genlight
gen<-vcfR2genlight(vcf)

#now have only california from CA & oregon, and wood from Nevada, CO, utah, NM, and AZ
#identify SNPs that are fixed away from the hybrid zone (Oregon vs. NM/CO)
mat<-extract.gt(vcf)
mat[1:5,1:5]
conv.mat<-mat
conv.mat[conv.mat == "0/0"]<-0
conv.mat[conv.mat == "0/1"]<-1
conv.mat[conv.mat == "1/1"]<-2
conv.mat<-as.data.frame(conv.mat)
#convert to numeric
for (i in 1:ncol(conv.mat)){
  conv.mat[,i]<-as.numeric(as.character(conv.mat[,i]))
}

#calc AF for each pop
cal.wood.af<-(rowSums(conv.mat[,c(1:32,47:63)], na.rm=T)/(rowSums(is.na(conv.mat[,c(1:32,47:63)]) == FALSE)))/2
isl.fla.af<-(rowSums(conv.mat[,c(33:46)], na.rm=T)/(rowSums(is.na(conv.mat[,c(33:46)]) == FALSE)))/2

#find fixed SNPs
diff<-abs((cal.wood.af)-(isl.fla.af))
#how many SNPs are fixed
table(is.na(diff) == FALSE & diff == 1)
vcf@fix[,1][(is.na(diff) == FALSE & diff == 1)]
conv.mat[(is.na(diff) == FALSE & diff == 1),]

#write vcf out to use for Fst calculation in vcftools
write.vcf(vcf, "~/Downloads/cal.wood.fst.vcf.gz")

#write pop files for each pop to feed into vcftools
write.table(locs$id[locs$species == "californica"], "~/Downloads/cali.txt", quote = F, row.names = F, col.names = F)
write.table(locs$id[locs$species == "woodhouseii"], "~/Downloads/wood.txt", quote = F, row.names = F, col.names = F)

#read in Fst file
fst<-read.table("~/Desktop/aph.data/introgress/cali.wood.fst/cal.wood.weir.fst", header = T)
fst$diff<-diff

fst$CHROM<-gsub("Pseudochr","",fst$CHROM,)
table(fst$CHROM)
fst$CHROM[fst$CHROM == "13_EQ832958_random"]<-"unplaced"
fst$CHROM[fst$CHROM == "21_EQ833162_random"]<-"unplaced"
fst$CHROM[fst$CHROM == "3_EQ832594_random"]<-"unplaced"
fst$CHROM[fst$CHROM == "4_EQ832640_random"]<-"unplaced"
fst$CHROM[fst$CHROM == "8_EQ832819_random"]<-"unplaced"
fst$CHROM[fst$CHROM == "scaffold_936"]<-"unplaced"
fst$CHROM[fst$CHROM == "Un_EQ835415"]<-"unplaced"
fst$CHROM[fst$CHROM == "Z_EQ833367_random"]<-"unplaced"

fst$CHROM[fst$CHROM == "1A"]<-"1.5"
fst$CHROM[fst$CHROM == "4A"]<-"4.5"
fst$CHROM[fst$CHROM == "M"]<-"30"
fst$CHROM[fst$CHROM == "Z"]<-"31"
table(fst$CHROM)
fst<-fst[fst$CHROM != "unplaced",]
fst$CHROM<-as.numeric(as.character(fst$CHROM))
table(fst$CHROM)

#make BPcum
nCHR <- length(unique(fst$CHROM))
ca$BPcum <- NA
s <- 0
nbp <- c()
for (i in levels(as.factor(fst$CHROM))){
  nbp[i] <- max(fst[fst$CHROM == i,]$POS)
  fst[fst$CHROM == i,"BPcum"] <- fst[fst$CHROM == i,"POS"] + s
  s <- s + nbp[i]
}


#set plotting params
axis.set <- fst %>% 
  dplyr::group_by(CHROM) %>% 
  dplyr::summarize(center = (max(BPcum) + min(BPcum)) / 2)

#
fst$WEIR_AND_COCKERHAM_FST[is.na(fst$WEIR_AND_COCKERHAM_FST)]<-0
fst$WEIR_AND_COCKERHAM_FST[fst$WEIR_AND_COCKERHAM_FST < 0]<-0

#isolate abbas
abbas<-fst[is.na(fst$diff) == FALSE & fst$diff == 1,]

#ggplot ca/wood
ggplot(fst, aes(x = BPcum, y =WEIR_AND_COCKERHAM_FST, 
               color = as.factor(CHROM))) +
  geom_point(alpha = 0.75, cex=.75) +
  #geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) +
  #scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values=c(rep(c("black","darkgrey"), times=14)))+
  scale_size_continuous(range = c(.2,1)) +
  labs(x = NULL, 
       y = "Fst") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))+
  geom_point(data=abbas, aes(x=BPcum, y=WEIR_AND_COCKERHAM_FST), col="red")


#
ggsave("~/Desktop/aph.data/introgress/cali.wood.fst/ca.wood.manhattan.pdf", units="in", width=8.5, height=1.7, dpi=300)


