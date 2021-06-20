

#define a function that will process the raw LDhelmet output and create a minimally sized dataframe per chrom
process.file<-function(x){
chrom<-readLines(x)
chrom<-chrom[-c(1,2)]
chrom[1]<-"left_snp right_snp mean p0.025 p0.975"
chrom<-read.table(textConnection(chrom), header = TRUE, stringsAsFactors = FALSE)
par(mfrow=c(1,2))
plot(chrom$left_snp,chrom$mean)
left_snp<-NULL
right_snp<-NULL
mean<-NULL
left_snp[1]<-chrom$left_snp[1]
f<-1
for (i in 1:(nrow(chrom)-1)){
  if (chrom$mean[i] == chrom$mean[i+1]){}
  else{right_snp[f]<-chrom$right_snp[i]
  mean[f]<-chrom$mean[i]
  left_snp[f+1]<-chrom$left_snp[i+1]
  f<-f+1
  }
}
right_snp[length(left_snp)]<-chrom$right_snp[nrow(chrom)]
mean[length(left_snp)]<-chrom$mean[nrow(chrom)]
chrom.red<-data.frame(left_snp=left_snp,right_snp=right_snp,mean=mean)
plot(chrom.red$left_snp,chrom.red$mean)
return(chrom.red)
}


chr.1<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr1_recombination_bpen100.txt.gz")
chr.1A<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr1A_recombination_bpen100.txt.gz")
chr.2<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr2_recombination_bpen100.txt.gz")
chr.3<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr3_recombination_bpen100.txt.gz")
chr.4<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr4_recombination_bpen100.txt.gz")
chr.4A<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr4A_recombination_bpen100.txt.gz")
chr.5<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr5_recombination_bpen100.txt.gz")
chr.6<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr6_recombination_bpen100.txt.gz")
chr.7<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr7_recombination_bpen100.txt.gz")
chr.8<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr8_recombination_bpen100.txt.gz")
chr.9<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr9_recombination_bpen100.txt.gz")
chr.10<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr10_recombination_bpen100.txt.gz")
chr.11<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr11_recombination_bpen100.txt.gz")
chr.12<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr12_recombination_bpen100.txt.gz")
chr.13<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr13_recombination_bpen100.txt.gz")
chr.14<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr14_recombination_bpen100.txt.gz")
chr.15<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr15_recombination_bpen100.txt.gz")
chr.16<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr16_recombination_bpen100.txt.gz")
chr.17<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr17_recombination_bpen100.txt.gz")
chr.18<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr18_recombination_bpen100.txt.gz")
chr.19<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr19_recombination_bpen100.txt.gz")
chr.20<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr20_recombination_bpen100.txt.gz")
chr.21<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr21_recombination_bpen100.txt.gz")
chr.22<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr22_recombination_bpen100.txt.gz")
chr.23<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr23_recombination_bpen100.txt.gz")
chr.24<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr24_recombination_bpen100.txt.gz")
chr.25<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr25_recombination_bpen100.txt.gz")
chr.26<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr26_recombination_bpen100.txt.gz")
chr.27<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr27_recombination_bpen100.txt.gz")
chr.28<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chr28_recombination_bpen100.txt.gz")
chr.Z<-process.file("/Users/devder/Downloads/doi_10.5061_dryad.fd24j__v1/DataDryad/recombination_maps/ZF/chrZ_recombination_bpen100.txt.gz")



#add chromosome column to each df
chr.1$chr<-rep(1, times=nrow(chr.1))
chr.1A<-rep("1A", times=nrow(chr.1A))
chr.2$chr<-rep(2, times=nrow(chr.2))
chr.3$chr<-rep(3, times=nrow(chr.3))
chr.4$chr<-rep(4, times=nrow(chr.4))
chr.4A$chr<-rep("4A", times=nrow(chr.4A))
chr.5$chr<-rep(5, times=nrow(chr.5))
chr.6$chr<-rep(6, times=nrow(chr.6))
chr.7$chr<-rep(7, times=nrow(chr.7))
chr.8$chr<-rep(8, times=nrow(chr.8))
chr.9$chr<-rep(9, times=nrow(chr.9))
chr.10$chr<-rep(10, times=nrow(chr.10))
chr.11$chr<-rep(11, times=nrow(chr.11))
chr.12$chr<-rep(12, times=nrow(chr.12))
chr.13$chr<-rep(13, times=nrow(chr.13))
chr.14$chr<-rep(14, times=nrow(chr.14))
chr.15$chr<-rep(15, times=nrow(chr.15))
chr.16$chr<-rep(16, times=nrow(chr.16))
chr.17$chr<-rep(17, times=nrow(chr.17))
chr.18$chr<-rep(18, times=nrow(chr.18))
chr.19$chr<-rep(19, times=nrow(chr.19))
chr.20$chr<-rep(20, times=nrow(chr.20))
chr.21$chr<-rep(21, times=nrow(chr.21))
chr.22$chr<-rep(22, times=nrow(chr.22))
chr.23$chr<-rep(23, times=nrow(chr.23))
chr.24$chr<-rep(24, times=nrow(chr.24))
chr.25$chr<-rep(25, times=nrow(chr.25))
chr.26$chr<-rep(26, times=nrow(chr.26))
chr.27$chr<-rep(27, times=nrow(chr.27))
chr.28$chr<-rep(28, times=nrow(chr.28))
chr.Z$chr<-rep("Z", times=nrow(chr.Z))

#combine all dfs
recomb.map<-as.data.frame(rbind(chr.1,chr.1A,chr.2,chr.3,chr.4,chr.4A,chr.5,chr.6,chr.7,
                                chr.8,chr.9,chr.10,chr.11,chr.12,chr.13,chr.14,chr.15,chr.16,
                                chr.17,chr.18,chr.19,chr.20,chr.21,chr.22,chr.23,chr.24,chr.25,
                                chr.26,chr.27,chr.28,chr.Z))

write.table(recomb.map, "~/Desktop/aph.data/ZF.recomb.map.txt", quote=F, row.names = F, col.names = T)







