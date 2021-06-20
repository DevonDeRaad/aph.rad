#18 September 2020
library(vcfR)
library(ggplot2)
library(gridExtra)
library(ggridges)
library(adegenet)
library(dplyr)
library(StAMPP)
library(gplots)
library(scatterpie)
library(geosphere)
library(ggtree)

#read in locality info for samples
locs<-read.csv("~/Desktop/aph.data/rad.sampling.csv")
#fix mislabeled sample
locs$id<-as.character(locs$id)
locs$id[locs$id == "A_californica_334171"]<-"A_woodhouseii_334171"

#read in vcf
vcfR <- read.vcfR("~/Desktop/aph.data/filtered.vcf.gz")
#fix mislabeled sample
colnames(vcfR@gt)[colnames(vcfR@gt) == "A_californica_334171"]<-"A_woodhouseii_334171"

#convert to genlight
gen<- vcfR2genlight(vcfR)
#reorder gen
gen<-gen[c(1:19,21:56,20,57:95),]

#subset locs to include only samples that passed filtering, and have been retained in the vcf
locs<-locs[locs$id %in% gen@ind.names,]
locs<-locs[c(1:19,21:56,20,57:95),]
locs[56,5]<-"woodhouseii"

table(locs$species)
table(locs$subspecies)
table(gen$chromosome)

#make a nuclear subset and a Z chromosome subset
gen.z<-gen[,gen$chromosome == "PseudochrZ"]
gen.nuc<-gen[,gen$chromosome != "PseudochrZ"]
gen.nuc<-gen[,gen$chromosome != "PseudochrZ_EQ833367_random"]
gen.nuc<-gen[,gen$chromosome != "PseudochrM"]

#compare pairwise divergence between subspecies in nuclear and Z SNPs
pop(gen.z)<-locs$subspecies
z.fst<-stamppFst(gen.z)
pop(gen.nuc)<-locs$subspecies
nuc.fst<-stamppFst(gen.nuc)
hist(nuc.fst$Fsts, breaks=10, xlim=c(0,1), col=rgb(0,0,1,1/4), main="hist pairwise Fst between subspecies nuclear & Z")
hist(z.fst$Fsts, breaks=20, xlim=c(0,1), col=rgb(1,0,0,1/4), add=TRUE)

#compare between species
pop(gen.z)<-locs$species
z.fst<-stamppFst(gen.z)
pop(gen.nuc)<-locs$species
nuc.fst<-stamppFst(gen.nuc)
hist(nuc.fst$Fsts, breaks=4, xlim=c(0,1), col=rgb(0,0,1,1/4), main="hist pairwise Fst between species nuclear & Z")
hist(z.fst$Fsts, breaks=4, xlim=c(0,1), col=rgb(1,0,0,1/4), add=TRUE)

#make eems of nuclear and Z data separately
#PCA of nuclear vs Z data for california only
cali.gen.nuc<-new("genlight",(as.matrix(gen.nuc)[locs$species == "californica", ]))
cali.gen.z<-new("genlight",(as.matrix(gen.z)[locs$species == "californica", ]))
#subset loc info
cali.locs<-locs[locs$species == "californica", ]
cali.locs<-droplevels(cali.locs)

nuc.pca<-glPca(cali.gen.nuc, nf=6)
#pull pca scores out of df
nuc.pca.scores<-as.data.frame(nuc.pca$scores)
nuc.pca.scores$subspecies<-cali.locs$subspecies
#ggplot color by subspecies
ggplot(nuc.pca.scores, aes(x=PC1, y=PC2, color=subspecies)) +
  geom_point(cex = 2)

#zpca
z.pca<-glPca(cali.gen.z, nf=6)
#pull pca scores out of df
z.pca.scores<-as.data.frame(z.pca$scores)
z.pca.scores$subspecies<-cali.locs$subspecies
#ggplot color by subspecies
ggplot(z.pca.scores, aes(x=PC1, y=PC2, color=subspecies)) +
  geom_point(cex = 2)

pop(cali.gen.z)<-cali.locs$subspecies
z.fst<-stamppFst(cali.gen.z)
z.fst$Fsts
pop(cali.gen.nuc)<-cali.locs$subspecies
nuc.fst<-stamppFst(cali.gen.nuc)
nuc.fst$Fsts

#PCA of nuclear vs Z data for woodhouse only
wood.gen.nuc<-new("genlight",(as.matrix(gen.nuc)[locs$species == "woodhouseii" | locs$species == "sumichrasti", ]))
wood.gen.z<-new("genlight",(as.matrix(gen.z)[locs$species == "woodhouseii" | locs$species == "sumichrasti", ]))
#subset loc info
wood.locs<-locs[locs$species == "woodhouseii" | locs$species == "sumichrasti", ]
wood.locs<-droplevels(wood.locs)

nuc.pca<-glPca(wood.gen.nuc, nf=6)
#pull pca scores out of df
nuc.pca.scores<-as.data.frame(nuc.pca$scores)
nuc.pca.scores$subspecies<-wood.locs$subspecies
#ggplot color by subspecies
ggplot(nuc.pca.scores, aes(x=PC1, y=PC2, color=subspecies)) +
  geom_point(cex = 2)

#zpca
z.pca<-glPca(wood.gen.z, nf=6)
#pull pca scores out of df
z.pca.scores<-as.data.frame(z.pca$scores)
z.pca.scores$subspecies<-wood.locs$subspecies
#ggplot color by subspecies
ggplot(z.pca.scores, aes(x=PC1, y=PC2, color=subspecies)) +
  geom_point(cex = 2)

pop(wood.gen.z)<-wood.locs$subspecies
z.fst<-stamppFst(wood.gen.z)
z.fst$Fsts
pop(wood.gen.nuc)<-wood.locs$subspecies
nuc.fst<-stamppFst(wood.gen.nuc)
nuc.fst$Fsts


write.table(gen@ind.names[locs$species == "insularis"], "island.txt", row.names = F, col.names = F, quote=F)
write.table(gen@ind.names[locs$subspecies == "texana"], "texas.txt", row.names = F, col.names = F, quote=F)
write.table(gen@ind.names[locs$species == "sumichrasti"], "sumi.txt", row.names = F, col.names = F, quote=F)
write.table(gen@ind.names[locs$species == "californica"], "cali.txt", row.names = F, col.names = F, quote=F)
write.table(gen@ind.names[locs$species == "woodhouseii"], "wood.txt", row.names = F, col.names = F, quote=F)

