library(vcfR)
library(adegenet)
library(StAMPP)


#read in vcfR
vcfR <- read.vcfR("~/Desktop/aph.data/unzipped.filtered.vcf")
#convert to gen
gen<-vcfR2genlight(vcfR)
#read in locs
locs<-read.csv("~/Desktop/aph.data/rad.sampling.csv")
#subset locs to include only samples that passed filtering, and have been retained in the vcf
locs<-locs[locs$id %in% gen@ind.names,]
locs[20,5]<-"woodhouseii"

#separate nuclear and Z
nuclear<-gen[,gen$chromosome != "PseudochrZ" & gen$chromosome != "PseudochrZ_EQ833367_random" & gen$chromosome != "PseudochrM"]
z<-gen[,gen$chromosome == "PseudochrZ" | gen$chromosome == "PseudochrZ_EQ833367_random"]

#define populations (requirement for the next line)
pop(nuclear)<-locs$species
pop(z)<-locs$species
#calculate pairwise genetic distance matrix among all samples in the genlight object
sample.div <- stamppNeisD(nuclear, pop = FALSE)
sample.div.z <- stamppNeisD(z, pop = FALSE)
rownames(sample.div)<-gsub(".*_", "", rownames(sample.div))
rownames(sample.div.z)<-gsub(".*_", "", rownames(sample.div.z))

#export for splitstree
stamppPhylip(distance.mat=sample.div, file="~/Downloads/nuclear.scrub.splits.txt")
stamppPhylip(distance.mat=sample.div.z, file="~/Downloads/z.scrub.splits.txt")

#write out these vcf files to pass to stacks to calculate pi for Z and autosomes separately in each of the six putatively delimited species
vcf.z<-vcfR[vcfR@fix[,1] == "PseudochrZ",]
vcf.nuc<-vcfR[vcfR@fix[,1] != "PseudochrZ" & vcfR@fix[,1] != "PseudochrZ_EQ833367_random" & vcfR@fix[,1] != "PseudochrM",]
write.vcf(vcf.nuc, "~/Downloads/nuclear.vcf.gz")
write.vcf(vcf.z, "~/Downloads/z.vcf.gz")

#make popmap for calculating pi for each of the 6 identified species using stacks
samps<-colnames(vcfR@gt)[2:96]
pops<-gsub("A_","",samps)
popmap<-data.frame(id=samps,spec=gsub("_.*","",pops), stringsAsFactors = FALSE)
#manually add texas as unique species
popmap[c(74:80),2]<-"texana"
write.table(popmap, "~/Downloads/six.species.popmap.txt", quote = F, row.names = F, col.names = F, sep = '\t')

