#make subsampled vcf files for snapp input

#read in linkage and quality filtered vcf
vcfR <- read.vcfR("~/Desktop/aph.data/unlinked.filtered.recode.vcf")
vcfR
#read in locs
locs<-read.csv("~/Desktop/aph.data/rad.sampling.csv")
#subset locs to include only samples that passed filtering, and have been retained in the vcf
locs<-locs[locs$id %in% gen@ind.names,]
colnames(vcfR@gt)

#define the new vcfs
vcfR1<-vcfR
vcfR2<-vcfR
vcfR3<-vcfR
vcfR4<-vcfR
vcfR5<-vcfR

samps<-colnames(vcfR@gt)
#make a vector that randomly samples 2 individuals from within each species group
sample.specs<-c(sample(c(1:20,22:33),size = 2),sample(c(34:39),size = 2),sample(c(40:47),size = 2),
  sample(c(48:57),size = 2),sample(c(58:95,21),size = 2))
 
#subset vcfR1 
vcfR1@gt <- vcfR@gt[,c(1,sample.specs)]

#resample
sample.specs<-c(sample(c(1:20,22:33),size = 2),sample(c(34:39),size = 2),sample(c(40:47),size = 2),
                sample(c(48:57),size = 2),sample(c(58:95,21),size = 2))

#subset vcfR2 
vcfR2@gt <- vcfR@gt[,c(1,sample.specs)]

#resample
sample.specs<-c(sample(c(1:20,22:33),size = 2),sample(c(34:39),size = 2),sample(c(40:47),size = 2),
                sample(c(48:57),size = 2),sample(c(58:95,21),size = 2))

#subset vcfR3
vcfR3@gt <- vcfR@gt[,c(1,sample.specs)]

#resample
sample.specs<-c(sample(c(1:20,22:33),size = 2),sample(c(34:39),size = 2),sample(c(40:47),size = 2),
                sample(c(48:57),size = 2),sample(c(58:95,21),size = 2))

#subset vcfR4 
vcfR4@gt <- vcfR@gt[,c(1,sample.specs)]

#resample
sample.specs<-c(sample(c(1:20,22:33),size = 2),sample(c(34:39),size = 2),sample(c(40:47),size = 2),
                sample(c(48:57),size = 2),sample(c(58:95,21),size = 2))

#subset vcfR5
vcfR5@gt <- vcfR@gt[,c(1,sample.specs)]

#write out each vcf for conversion
write.vcf(vcfR1, "~/Desktop/aph.data/snapp/rep1.vcf.gz")
write.vcf(vcfR2, "~/Desktop/aph.data/snapp/rep2.vcf.gz")
write.vcf(vcfR3, "~/Desktop/aph.data/snapp/rep3.vcf.gz")
write.vcf(vcfR4, "~/Desktop/aph.data/snapp/rep4.vcf.gz")
write.vcf(vcfR5, "~/Desktop/aph.data/snapp/rep5.vcf.gz")









