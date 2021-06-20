#make subsampled vcf files for snapp input

#read in linkage and quality filtered vcf
vcfR <- read.vcfR("~/Desktop/aph.data/unlinked.filtered.recode.vcf")
#vcfR <- read.vcfR("~/Desktop/aph.data/populations.snps.vcf")

vcfR
#read in locs
locs<-read.csv("~/Desktop/aph.data/rad.sampling.csv")
#subset locs to include only samples that passed filtering, and have been retained in the vcf
locs<-locs[locs$id %in% colnames(vcfR@gt),]
colnames(vcfR@gt)
rownames(locs)<-1:95

#loop to generate randomly subsampled vcfs
for (i in 1:5){
#randomly sample 2 samples from each of the K=7 species groups
sample.specs<-c(sample(c(2:20,22:33),size = 2),sample(c(34:39),size = 2),sample(c(40:47),size = 2),
                sample(c(48:57),size = 2),sample(c(58:74,21),size = 2),sample(c(75:81),size = 2),
                sample(c(82:96),size = 2))
#duplicate full vcfR
vcfR.sub<-vcfR
#subset the full vcf and assign it to the new empty vcfR
vcfR.sub@gt <- vcfR@gt[,c(1,sample.specs)]
#write to an output file
#write.vcf(vcfR.sub, paste0("~/Desktop/aph.data/snapp/7spec/rep",i,".vcf.gz"))
}

#bring this into beauti and prepare five separate .xml input files for snapp


#read in nexus file containing all 95 samples
dna<-read.nexus.data("~/Desktop/aph.data/snapp/7spec/binary.nex")
names(dna)
#loop to generate randomly subsampled vcfs
for (i in 1:5){
  #randomly sample 2 samples from each of the K=7 species groups
  sample.specs<-c(sample(c(1:19,21:32),size = 2),sample(c(33:38),size = 2),sample(c(39:46),size = 2),
                  sample(c(47:56),size = 2),sample(c(57:73,20),size = 2),sample(c(74:80),size = 2),
                  sample(c(81:95),size = 2))
  #subset the full nex and assign it to the new empty vcfR
  dna.sub <- dna[sample.specs]
  #write to an output file
#  write.nexus.data(dna.sub, paste0("~/Desktop/aph.data/snapp/7spec/rep",i,".nex"))
  #clear file
  rm(dna.sub)
}

vcfR@gt[1,]
vcfR@fix


vcfR@fix<-vcfR@fix[rowSums(is.na(vcfR@gt)) == 0,]
vcfR@gt<-vcfR@gt[rowSums(is.na(vcfR@gt)) == 0,]
vdna<-vcfR2DNAbin(vcfR, extract.haps = F, consensus = T)
vdna 
write.nexus.data(x = vdna, file="~/Downloads/fuck.miss.unlinked.nex", format = "dna", missing = "n")
vcfR

#read in nex
dna<-read.nexus.data("~/Desktop/aph.data/snapp/7spec/binary.nex")
names(dna)

#subsample nexus to one sample per locality
dna.sub <- dna[c(1,6,10,13,17,22,25,26,29,30,33,40,50,55,57,60,62,65,67,68,71,77,82,83,84,87)]

write.nexus.data(x = dna.sub, file="~/Desktop/aph.data/snapp/nex.by.locality.nex", format = "dna", missing = "?")








#bring in full vcf file
#read in vcf as vcfR
vcfR <- read.vcfR("~/Desktop/aph.data/unzipped.filtered.vcf")
dim(vcfR@gt)

#
head(vcfR@fix)
levels(as.factor(vcfR@fix[,1]))

j=10000
fix<-as.numeric(vcfR@fix[,2])
g<-c()

for (t in 1:length(levels(as.factor(vcfR@fix[,1])))){
  fix.sub<-fix[vcfR@fix[,1] == levels(as.factor(vcfR@fix[,1]))[t]]
  fix.sub<-fix.sub[order(fix.sub)]
  prev<-fix.sub[1]
  k<-c()
  k[1]<-TRUE
    for (i in 2:length(fix.sub)){
      k[i]<- fix[i] > prev+j
        if (fix[i] > prev+j){
          prev<-fix[i]
        }
    }
  keep<-fix.sub[k]
  g<-c(g, fix.sub %in% keep)
}

vcfR@fix<-vcfR@fix[g,]
dim(vcfR@fix)
vcfR@gt<-vcfR@gt[g,]
vcfR

#check
fix.sub<-as.numeric(vcfR@fix[,2][vcfR@fix[,1] == levels(as.factor(vcfR@fix[,1]))[1]])
z<-c()
k=1
for (i in 2:length(fix.sub)){
  z[k]<-fix.sub[i]-fix.sub[i-1]
  k=k+1
}
hist(z, breaks=10000, xlim=c(0,10000))


#write nexus
vdna<-vcfR2DNAbin(vcfR, extract.haps = F, consensus = T)
vdna 
write.nexus.data(x = vdna, file="~/Downloads/fuck.miss.unlinked.nex", format = "dna", missing = "n")



