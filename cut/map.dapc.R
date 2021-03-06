#geographic investigation
#18 September 2020
library(vcfR)
library(ggplot2)
library(gridExtra)
library(ggridges)
library(adegenet)
library(dplyr)
library(StAMPP)
library(gplots)

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
gen

#subset locs to include only samples that passed filtering, and have been retained in the vcf
locs<-locs[locs$id %in% gen@ind.names,]

table(locs$species)
table(locs$subspecies)
table(locs$decimallatitude)


#split df by species
spec.dfs<-split(locs, locs$species)

#init sampling.df which will be a df of samples grouped by unique lat/long
sampling.df<-data.frame(NULL)
for (i in names(spec.dfs)){
  samps<-spec.dfs[[i]] %>% group_by(decimallatitude, decimallongitude) %>% summarize(count=n())
  df<-cbind(rep(i, times=nrow(samps)), samps)
  sampling.df<-as.data.frame(rbind(sampling.df, df))
}
#fix colnames
colnames(sampling.df)<-c("species","lat","long","count")

#make map#
pac<-map_data("world")
#"palegreen","palegreen4","powderblue","skyblue1","royalblue3","navyblue","lightsalmon","sienna3")
#make full map
pac<-map_data("world")
#
ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(-123, -80), ylim = c(14, 47)) + 
  geom_point(data = sampling.df, aes(x = long, y = lat, col=species, size=count), alpha =.9, show.legend=TRUE) +
  theme_classic()+
  scale_color_manual(values=funky(5))+
  scale_size_continuous(range = c(2,8))+
  guides(colour = guide_legend(override.aes = list(size = 4), order=1, label.theme = element_text(face = "italic")),
         size = guide_legend(nrow = 1, order = 2))+
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0.01, 0.01))

#Do dapc and choose groups
#assign samples to the number of groups present in popmap, retain all PCAs
#grp<-find.clusters(gen)
#set manually the values I chose based on looking at the scree plots
grp<-find.clusters(gen, n.pca=100, n.clust=6)

#check how well that assignment matched up to the provided popmap
samps<-merge(popmap, data.frame(group=grp$grp, id=labels(grp$grp)), by='id')
table(samps$pop, samps$group)

#run dapc, retain all discriminant axes, and enough PC axes to explain 75% of variance
#dapc1<-dapc(gen, grp$grp)
#set manually the values I chose from scree plots
dapc1<-dapc(gen, grp$grp, n.da = 5, n.pca =10)

#plot in two dimensions
scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=funky(6), solid=.4,
        cex=5,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:6))

#plot compoplot
compoplot(dapc1, legend=FALSE, col=funky(6), show.lab =TRUE, cex.names=.4)

#plot map colored by species
ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(-123, -80), ylim = c(14, 47)) + 
  geom_point(data = locs, aes(x = decimallongitude, y = decimallatitude, col=species), alpha =.9, size=3, show.legend=TRUE) +
  theme_classic()+
  scale_color_manual(values=funky(5))+
  guides(colour = guide_legend(override.aes = list(size = 4), order=1, label.theme = element_text(face = "italic")),
         size = guide_legend(nrow = 1, order = 2))+
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0.01, 0.01))

#plot map colored by group assignment
ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(-123, -80), ylim = c(14, 47)) + 
  geom_point(data = locs, aes(x = decimallongitude, y = decimallatitude, col=grp$grp), alpha =.9, size=3, show.legend=TRUE) +
  theme_classic()+
  scale_color_manual(values=funky(6))+
  guides(colour = guide_legend(override.aes = list(size = 4), order=1, label.theme = element_text(face = "italic")),
         size = guide_legend(nrow = 1, order = 2))+
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0.01, 0.01))

#calculate divergence btwn inds and species
gen@pop<-grp$grp
inds<-stamppNeisD(gen, pop = FALSE)
pops<-stamppNeisD(gen, pop = TRUE)

#plot heatmap by group assignment
heatmap.2(pops, trace="none", cexRow=0.8, cexCol=0.4)
#plot nj tree unrooted
plot(nj(pops), type = "unrooted", cex = 1.5)

#plot heatmap by individual
heatmap.2(inds, trace="none", cexRow=0.4, cexCol=0.4)
plot(nj(inds), type = "unrooted", cex = .4)

#plot map colored by group assignment
ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(-123, -80), ylim = c(14, 47)) + 
  geom_point(data = locs, aes(x = decimallongitude, y = decimallatitude, col=grp$grp), alpha =.9, size=3, show.legend=TRUE) +
  theme_classic()+
  scale_color_manual(values=funky(6))+
  guides(colour = guide_legend(override.aes = list(size = 4), order=1, label.theme = element_text(face = "italic")),
         size = guide_legend(nrow = 1, order = 2))+
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0.01, 0.01))

#calculate pairwise Fst between groups identified in dapc
fst<-stamppFst(gen, nboots = 1, percent =95)
fst

#subset genlight to only california scrub jays
cali.gen<-new("genlight",(as.matrix(gen)[locs$species == "californica", ]))
#subset loc info
cali.locs<-locs[locs$species == "californica", ]

#subset genlight to only woodhouse's scrub jays
wood.gen<-new("genlight",(as.matrix(gen)[locs$species == "woodhouseii" | locs$species == "sumichrasti", ]))
#subset loc info
wood.locs<-locs[locs$species == "woodhouseii" | locs$species == "sumichrasti", ]

#PCA of cali
cali.pca<-glPca(cali.gen, nf=6)

#pull pca scores out of df
cali.pca.scores<-as.data.frame(cali.pca$scores)

#plot map of california colored by subspecies
ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(-125, -105), ylim = c(20, 50)) + 
  geom_point(data = cali.locs, aes(x = decimallongitude, y = decimallatitude, col=subspecies), alpha =.9, size=3, show.legend=TRUE) +
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size = 4), order=1, label.theme = element_text(face = "italic")),
         size = guide_legend(nrow = 1, order = 2))+
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0.01, 0.01))

#ggplot color by subspecies
ggplot(cali.pca.scores, aes(x=PC1, y=PC2, color=cali.locs$subspecies)) +
  geom_point(cex = 2)

#ggplot compare PC1 to latitude
ggplot(cali.pca.scores, aes(x=PC1, y=cali.locs$decimallatitude, color=cali.locs$subspecies)) +
  geom_point(cex = 2)

#ggplot compare PC2 to latitude
ggplot(cali.pca.scores, aes(x=PC2, y=cali.locs$decimallatitude, color=cali.locs$subspecies)) +
  geom_point(cex = 2)


#PCA of wood
wood.pca<-glPca(wood.gen, nf=6)

#pull pca scores out of df
wood.pca.scores<-as.data.frame(wood.pca$scores)

#plot map of woodfornia colored by subspecies
ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(-120, -95), ylim = c(15, 45)) + 
  geom_point(data = wood.locs, aes(x = decimallongitude, y = decimallatitude, col=subspecies), alpha =.9, size=3, show.legend=TRUE) +
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size = 4), order=1, label.theme = element_text(face = "italic")),
         size = guide_legend(nrow = 1, order = 2))+
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0.01, 0.01))

#ggplot color by subspecies
ggplot(wood.pca.scores, aes(x=PC1, y=PC2, color=wood.locs$subspecies)) +
  geom_point(cex = 2)

#ggplot compare PC3 to latitude
ggplot(wood.pca.scores, aes(x=PC3, y=wood.locs$decimallatitude, color=wood.locs$subspecies)) +
  geom_point(cex = 2)

#two intermediate/admixed samples: "A_sumichrasti_343513", "A_sumichrasti_393636"
#calculate missing data in all samples to see if this lack of cohesion is due to true admixture or a lack of power to categorize these samples due to missing data
dp<- extract.gt(vcfR, element='DP', as.numeric=TRUE)
#calculate missingness by individual
miss<-colSums(is.na(dp))/nrow(dp)
miss[locs$species == "sumichrasti"]

#pca of woodhouse with only SNPs with no missing data
wood.gen<-new("genlight", (as.matrix(wood.gen))[,(colSums(is.na(as.matrix(wood.gen))) < 1)])
wood.gen

#PCA
wood.pca<-glPca(wood.gen, nf=6)
#pull pca scores out of df
wood.pca.scores<-as.data.frame(wood.pca$scores)

#ggplot color by subspecies
ggplot(wood.pca.scores, aes(x=PC1, y=PC2, color=wood.locs$subspecies)) +
  geom_point(cex = 2)+
  ggtitle("No missing data")


