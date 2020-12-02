library(tidyr)
library(dplyr)

#Fst comparisons
setwd("~/Desktop/aph.data/")
recomb.map<-read.table("ZF.recomb.map.txt", header=T)

cali.wood<-read.table("cali.wood.weir.fst", header=T)
cali.island<-read.table("cali.island.weir.fst", header=T)
wood.sumi<-read.table("wood.sumi.weir.fst", header=T)
wood.tex<-read.table("wood.tex.weir.fst", header=T)

fsts<-as.data.frame(rbind(cali.wood,cali.island,wood.sumi,wood.tex))
fsts$comp<-c(rep("cali.wood",times=nrow(cali.wood)),
             rep("cali.island",times=nrow(cali.wood)),
             rep("wood.sumi",times=nrow(cali.wood)),
             rep("wood.tex",times=nrow(cali.wood)))

fsts$WEIR_AND_COCKERHAM_FST[fsts$WEIR_AND_COCKERHAM_FST < 0]<-0

ggplot(fsts, aes(WEIR_AND_COCKERHAM_FST, colour = comp)) +
  geom_freqpoly()+
  theme_classic()

ggplot(fsts, aes(WEIR_AND_COCKERHAM_FST, fill = comp)) +
  geom_histogram(bins=50)+
  theme_classic()

ggplot(fsts, aes(WEIR_AND_COCKERHAM_FST, after_stat(density), colour = comp)) +
  geom_freqpoly(bins=50)+
  theme_classic()


ca<-data.frame(CHROM=cali.wood$CHROM, POS=cali.wood$POS, ca.wood.fst=cali.wood$WEIR_AND_COCKERHAM_FST,
               ca.island.fst=cali.island$WEIR_AND_COCKERHAM_FST, wood.tex.fst=wood.tex$WEIR_AND_COCKERHAM_FST,
               wood.sumi.fst=wood.sumi$WEIR_AND_COCKERHAM_FST)
ca$ca.wood.fst[ca$ca.wood.fst < 0]<-0
ca$ca.island.fst[ca$ca.island.fst < 0]<-0
ca$wood.tex.fst[ca$wood.tex.fst < 0]<-0
ca$wood.sumi.fst[ca$wood.sumi.fst < 0]<-0


ca$chrom<- NA
ca$chrom[ca$CHROM == "Pseudochr1"]<-"1"
ca$chrom[ca$CHROM == "Pseudochr1A"]<-"1A"
ca$chrom[ca$CHROM == "Pseudochr2"]<-"2"
ca$chrom[ca$CHROM == "Pseudochr3"]<-"3"
ca$chrom[ca$CHROM == "Pseudochr4"]<-"4"
ca$chrom[ca$CHROM == "Pseudochr4A"]<-"4A"
ca$chrom[ca$CHROM == "Pseudochr5"]<-"5"
ca$chrom[ca$CHROM == "Pseudochr6"]<-"6"
ca$chrom[ca$CHROM == "Pseudochr7"]<-"7"
ca$chrom[ca$CHROM == "Pseudochr8"]<-"8"
ca$chrom[ca$CHROM == "Pseudochr9"]<-"9"
ca$chrom[ca$CHROM == "Pseudochr10"]<-"10"
ca$chrom[ca$CHROM == "Pseudochr11"]<-"11"
ca$chrom[ca$CHROM == "Pseudochr12"]<-"12"
ca$chrom[ca$CHROM == "Pseudochr13"]<-"13"
ca$chrom[ca$CHROM == "Pseudochr14"]<-"14"
ca$chrom[ca$CHROM == "Pseudochr15"]<-"15"
ca$chrom[ca$CHROM == "Pseudochr17"]<-"17"
ca$chrom[ca$CHROM == "Pseudochr18"]<-"18"
ca$chrom[ca$CHROM == "Pseudochr19"]<-"19"
ca$chrom[ca$CHROM == "Pseudochr20"]<-"20"
ca$chrom[ca$CHROM == "Pseudochr21"]<-"21"
ca$chrom[ca$CHROM == "Pseudochr23"]<-"23"
ca$chrom[ca$CHROM == "Pseudochr24"]<-"24"
ca$chrom[ca$CHROM == "Pseudochr26"]<-"26"
ca$chrom[ca$CHROM == "Pseudochr27"]<-"27"
ca$chrom[ca$CHROM == "PseudochrZ"]<-"Z"
ca$chrom[ca$CHROM == "PseudochrM"]<-"Mt"

ca<-ca[is.na(ca$chrom) == FALSE,]
ca$chrom<-as.factor(ca$chrom)
levels(ca$chrom)

ca$nuc<-NA
ca$nuc[ca$chrom != "Z"]<-"nuclear"
ca$nuc[ca$chrom == "Z"]<-"Z"
ggplot(ca, aes(ca.wood.fst, after_stat(density), colour = nuc)) +
  geom_freqpoly(bins=50)+
  theme_classic()

ggplot(ca, aes(ca.wood.fst, after_stat(count), colour = nuc)) +
  geom_freqpoly(bins=50)+
  theme_classic()

#make by chromosome mean Fsts
ca.wood.chroms<-aggregate(ca$ca.wood.fst[is.na(ca$ca.wood.fst) == FALSE], list(ca$chrom[is.na(ca$ca.wood.fst) == FALSE]), mean)
ca.island.chroms<-aggregate(ca$ca.island.fst[is.na(ca$ca.island.fst) == FALSE], list(ca$chrom[is.na(ca$ca.island.fst) == FALSE]), mean)
wood.tex.chroms<-aggregate(ca$wood.tex.fst[is.na(ca$wood.tex.fst) == FALSE], list(ca$chrom[is.na(ca$wood.tex.fst) == FALSE]), mean)
wood.sumi.chroms<-aggregate(ca$wood.sumi.fst[is.na(ca$wood.sumi.fst) == FALSE], list(ca$chrom[is.na(ca$wood.sumi.fst) == FALSE]), mean)

#fix up recomb map
recomb.map$left_snp<-as.numeric(as.character(recomb.map$left_snp))
recomb.map$right_snp<-as.numeric(as.character(recomb.map$right_snp))
recomb.map<-recomb.map[is.na(recomb.map$right_snp) == FALSE,]

#map recombination rate onto the Fsts we already calculated
caz<-ca[ca$chrom %in% recomb.map$chr,]
caz$chrom<-droplevels(caz$chrom)
recomb.map<-recomb.map[recomb.map$chr %in% ca$chrom,]
recomb.map$chr<-droplevels(recomb.map$chr)
recomb.map$mean<-as.numeric(as.character(recomb.map$mean))
z<-NULL
for (i in 1:nrow(caz)){
  if(length(recomb.map$mean[recomb.map$chr == caz$chrom[i] & recomb.map$left_snp < caz$POS[i] & recomb.map$right_snp > caz$POS[i] ]) == 0){
    z[i]<-NA
  }
  else{
    z[i]<-recomb.map$mean[recomb.map$chr == caz$chrom[i] & recomb.map$left_snp < caz$POS[i] & recomb.map$right_snp > caz$POS[i] ]
  }
}

caz$recomb<-z
car<-caz[is.na(caz$recomb) == FALSE,]

#test for correlation btwn recombination rate and divergence
#ca.wood
cor.test(log(car$recomb), car$ca.wood.fst)
ggplot(car, aes(x=log(recomb), y=ca.wood.fst, col=nuc))+
  geom_point(alpha=.3, cex=2.5)+
  scale_color_manual(values = c("black","red"))+
  theme_classic()+
  geom_smooth(method=lm,  linetype="dashed",
              color="darkred", fill="grey")

#ca.island
cor.test(log(car$recomb), car$ca.island.fst)
ggplot(car, aes(x=log(recomb), y=ca.island.fst, col=nuc))+
  geom_point(alpha=.3, cex=2.5)+
  scale_color_manual(values = c("black","red"))+
  theme_classic()+
  geom_smooth(method=lm,  linetype="dashed",
              color="darkred", fill="grey")

#wood.sumi
cor.test(log(car$recomb), car$wood.sumi.fst)
ggplot(car, aes(x=log(recomb), y=wood.sumi.fst, col=nuc))+
  geom_point(alpha=.3, cex=2.5)+
  scale_color_manual(values = c("black","red"))+
  theme_classic()+
  geom_smooth(method=lm,  linetype="dashed",
              color="darkred", fill="grey")

#wood.tex
cor.test(log(car$recomb), car$wood.tex.fst)
ggplot(car, aes(x=log(recomb), y=wood.tex.fst, col=nuc))+
  geom_point(alpha=.3, cex=2.5)+
  scale_color_manual(values = c("black","red"))+
  theme_classic()+
  geom_smooth(method=lm,  linetype="dashed",
              color="darkred", fill="grey")



