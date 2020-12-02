#geographic investigation
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
library(ggplotify)

#read in vcf
vcfR <- read.vcfR("~/Desktop/aph.data/unzipped.filtered.vcf")
vcfR
#convert to genlight
gen<- vcfR2genlight(vcfR)
#reorder gen to match species assignments
gen<-gen[c(33:46,1:19,21:32,57:73,81:95,20,74:80,47:56),]

#read in locality info for samples
locs<-read.csv("~/Desktop/aph.data/rad.sampling.csv")
#subset locs to include only samples that passed filtering, and have been retained in the vcf
locs<-locs[locs$id %in% gen@ind.names,]
#order by species assingments
locs<-locs[c(33:46,1:19,21:32,57:73,81:95,20,74:80,47:56),]
rownames(locs)<-1:95
#combine lat/longs for nearby sampling localities 
locs[20,c(2,3)]<-locs[21,c(2,3)]
locs[34,c(2,3)]<-locs[37,c(2,3)]
locs[43,c(2,3)]<-locs[44,c(2,3)]
locs[79,c(2,3)]<-locs[80,c(2,3)]
locs[70,c(2,3)]<-locs[68,c(2,3)]
locs[71,c(2,3)]<-locs[68,c(2,3)]
locs[74,c(2,3)]<-locs[64,c(2,3)]

locs$species<-as.character(locs$species)
locs[c(79:85),4]<-c(rep("texana", times=7))
locs$species<-as.factor(locs$species)

table(locs$species)
table(locs$subspecies)
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

#make full map
pac<-map_data("world")
#
ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(-123, -80), ylim = c(14, 47)) + 
  geom_point(data = sampling.df, aes(x = long, y = lat, col=species, size=count), alpha =.9, show.legend=TRUE) +
  theme_classic()+
  scale_color_manual(values=funky(5))+
  scale_size_continuous(range = c(4,8))+
  guides(colour = guide_legend(override.aes = list(size = 4), order=1, label.theme = element_text(face = "italic")),
         size = guide_legend(nrow = 1, order = 2))+
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0.01, 0.01))

#Do dapc and choose groups
#assign samples to the number of groups present in popmap, retain all PCAs
#grp<-find.clusters(gen, max.n.clust=15)
#set manually the values I chose based on looking at the scree plots
grp<-find.clusters(gen, n.pca=100, n.clust=6)
#run dapc, retain all discriminant axes, and enough PC axes to explain 75% of variance
#dapc1<-dapc(gen, grp$grp)
#set manually the values I chose from scree plots
dapc1<-dapc(gen, grp$grp, n.pca = 6, n.da = 100)

#create tidy format df for plotting the compoplot as a ggplot object
#fix names to order plot alphabetically
names<-rownames(dapc1$posterior)
names[79:85]<-c("A_texana_334230","A_texana_334241","A_texana_334242","A_texana_334243",
                "A_texana_334244","A_texana_334247","A_texana_334250")
#create df
dapc.assignments<-data.frame(sample=rep(names, times=6),
                             cluster=c(rep(1, times=95),rep(2, times=95),rep(3, times=95),
                                       rep(4, times=95),rep(5, times=95),rep(6, times=95)),
                             proportion=c(dapc1$posterior[,1],dapc1$posterior[,2],dapc1$posterior[,3],
                                          dapc1$posterior[,4],dapc1$posterior[,5],dapc1$posterior[,6]))
dapc.assignments$cluster<-as.factor(dapc.assignments$cluster)

ggplot(dapc.assignments, aes(fill=cluster, y=proportion, x=sample))+ 
  geom_bar(position="fill", stat="identity")+
  theme(axis.text.x = element_text(angle = 90))
  
#make bar plot
com<-ggplot(dapc.assignments, aes(fill=cluster, y=proportion, x=sample))+ 
  geom_bar(position="fill", stat="identity")+
  theme_classic()+
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank())+
  labs(x="")+
  #scale_x_discrete(breaks = c("A_californica_334012","A_coerulescens_396256","A_insularis_334032",
  #                            "A_sumichrasti_343512","A_texana_334243","A_woodhouseii_334211"),
  #                 label = c("California","Florida","Island","Sumichrast","Texas","Woodhouse")) +
  scale_fill_manual(name="Cluster",
                     limits = c(1,6,3,4,2,5),
                     values=c("red","blue","green","orange","purple","pink"),
                     labels = c("Island","California","Woodhouse","Texas","Sumichrast","Florida"))+
  labs(y="Assignment\nProbability")
  
#make scatterpie plot
post<-as.data.frame(dapc1$posterior)
post$id<-rownames(post)
plot.pies<-merge(post, locs, by= "id")
#
map<-ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(-123, -80), ylim = c(14, 47)) + 
  geom_scatterpie(aes(x=decimallongitude, y=decimallatitude,size=decimallatitude), 
                  data = plot.pies, cols = colnames(plot.pies[,c(2:7)]), lwd=.2, pie_scale = 1.1)+
  scale_fill_manual(name="Cluster",
                    limits = c(1,6,3,4,2,5),
                    values=c("red","blue","green","orange","purple","pink"),
                    labels = c("Island","California","Woodhouse","Texas","Sumichrast","Florida"))+
  theme_classic()+
  theme(legend.position="none")+
  labs(x="Longitude",y="Latitude")

#plot DF axis loadings
pca.loadings<-as.data.frame(dapc1$ind.coord)
df1<-ggplot(data=pca.loadings, aes(x=LD1, y=LD2, color=grp$grp))+
  geom_point(alpha=.75, cex=5)+
  theme_classic()+
  labs(x="Discriminant Function 1",y="Discriminant Function 2")+
  scale_color_manual(name="Cluster",
                     limits = c(1,6,3,4,2,5),
                     values=c("red","blue","green","orange","purple","pink"),
                     labels = c("Island","California","Woodhouse","Texas","Sumichrast","Florida"))+
  theme(legend.position = "none")

df2<-ggplot(data=pca.loadings, aes(x=LD3, y=LD4, color=grp$grp))+
  geom_point(alpha=.75, cex=5)+
  theme_classic()+
  labs(x="Discriminant Function 3",y="Discriminant Function 4")+
  scale_color_manual(name="Cluster",
                     limits = c(1,6,3,4,2,5),
                     values=c("red","blue","green","orange","purple","pink"),
                     labels = c("Island","California","Woodhouse","Texas","Sumichrast","Florida"))+
  theme(legend.position = "none")
  

lay <- rbind(c(1,1),
             c(2,3))
grid.arrange(map, df1, df2, layout_matrix = lay, heights=c(2.5,1))


#make same figure and add compoplot
lay <- rbind(c(1,1),
             c(2,2),
             c(3,4))
grid.arrange(map, com, df1, df2, layout_matrix = lay, heights=c(2,.5,1.2))


#logic for calculating shared, private, and fixed SNPs
#if target AF == 1, and rest AF==0, fixed non-ref
#if target AF == 0, and rest AF==1, fixed ref
#if target AF > 0 & < 1, and rest AF==0, private non-ref
#if target AF > 0 & < 1, and rest AF==1, private ref
#if rest AF > 0 & < 1, shared

#generate popmap file. Two column popmap with the same format as stacks, and the columns must be named 'id' and 'pop'
popmap<-data.frame(id=colnames(vcfR@gt)[2:length(colnames(vcfR@gt))],pop=substr(colnames(vcfR@gt)[2:length(colnames(vcfR@gt))], 3,11))
popmap$pop<-as.character(popmap$pop)
popmap[74:80,2]<-rep("texas", times=7)
popmap$pop<-as.factor(popmap$pop)

fix.private.share<-function(vcfR, popmap){

  #make vcfR into gt matrix
  gt.matrix<-extract.gt(vcfR)
  #convert values for easier computation
  gt.matrix[gt.matrix=="0/0"]<-0
  gt.matrix[gt.matrix=="0/1"]<-1
  gt.matrix[gt.matrix=="1/1"]<-2
  gt.matrix<-as.data.frame(gt.matrix)
  for (i in 1:ncol(gt.matrix)){
    gt.matrix[,i]<-as.numeric(as.character(gt.matrix[,i]))
  }
  
  #open df to fill up with information by pop
  df<-data.frame(pop=character(), class=character(), snps=numeric())
  #separate gt.matrix based on the population of interest
  for (i in 1:length(levels(popmap$pop))){
  target<-rowSums(gt.matrix[,as.character(popmap$id[popmap$pop==levels(popmap$pop)[i]])], na.rm = T)/
    (rowSums(!is.na(gt.matrix[,as.character(popmap$id[popmap$pop==levels(popmap$pop)[i]])]))*2)
  rest<-rowSums(gt.matrix[,as.character(popmap$id[popmap$pop!=levels(popmap$pop)[i]])], na.rm = T)/
    (rowSums(!is.na(gt.matrix[,as.character(popmap$id[popmap$pop!=levels(popmap$pop)[i]])]))*2)
  #remove sites where AF can't be calculated due to missing data
  missing.sites<-c(!is.na(target) & !is.na(rest))
  rest<-rest[missing.sites]
  target<-target[missing.sites]
  
  ##calculate # of fixed SNPs
  fix<-sum(abs(target-rest) == 1) #488
  
  #calc # of private non ref 
  y<-0
  for (j in 1:length(target)){
    if(target[j] > 0 & target[j] < 1 & rest[j] == 0){
      y<-y+1
    }
  }
  #calc # of private ref 
  for (k in 1:length(target)){
    if(target[k] > 0 & target[k] < 1 & rest[k] == 1){
      y<-y+1
    }
  }
  
  #calc # of shared SNPs
  z<-0
  for (m in 1:length(target)){
    if(rest[m] > 0 & rest[m] < 1){
      z<-z+1
    }
  }
  
  #combine this pop's info
  x<-cbind(rep(levels(popmap$pop)[i], times=3),
           c("fixed","private","shared"),
           c(fix,y,z))
  
  #add it to communal tidy df
  df<-as.data.frame(rbind(df,x))
  }
  #fix df colnames
  colnames(df)<-c("pop","class","snps")
  return(df)
}

#calc pie df
pie.df<-fix.private.share(vcfR, popmap)

#plot pies
#tex
tex<-ggplot(pie.df[pie.df$pop == "texas",], aes(x="", y=as.numeric(as.character(snps)), fill=class)) +
  geom_bar(stat="identity", width=1, color="black")+
  coord_polar("y", start=0)+
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values=c("black","grey","white"))
#cali
cal<-ggplot(pie.df[pie.df$pop == "californi",], aes(x="", y=as.numeric(as.character(snps)), fill=class)) +
  geom_bar(stat="identity", width=1, color="black")+
  coord_polar("y", start=0)+
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values=c("black","grey","white"))
#wood
wood<-ggplot(pie.df[pie.df$pop == "woodhouse",], aes(x="", y=as.numeric(as.character(snps)), fill=class)) +
  geom_bar(stat="identity", width=1, color="black")+
  coord_polar("y", start=0)+
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values=c("black","grey","white"))
#Fla
fla<-ggplot(pie.df[pie.df$pop == "coerulesc",], aes(x="", y=as.numeric(as.character(snps)), fill=class)) +
  geom_bar(stat="identity", width=1, color="black")+
  coord_polar("y", start=0)+
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values=c("black","grey","white"))
#Island
isl<-ggplot(pie.df[pie.df$pop == "insularis",], aes(x="", y=as.numeric(as.character(snps)), fill=class)) +
  geom_bar(stat="identity", width=1, color="black")+
  coord_polar("y", start=0)+
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values=c("black","grey","white"))
#sumi
sumi<-ggplot(pie.df[pie.df$pop == "sumichras",], aes(x="", y=as.numeric(as.character(snps)), fill=class)) +
  geom_bar(stat="identity", width=1, color="black")+
  coord_polar("y", start=0)+
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values=c("black","grey","white"))

#plot legend
legend <- cowplot::get_legend(ggplot(pie.df[pie.df$pop == "sumichras",], aes(x="", y=as.numeric(as.character(snps)), fill=class)) +
                                geom_bar(stat="identity", width=1, color="black")+
                                coord_polar("y", start=0)+
                                theme(legend.title = element_blank())+
                                scale_fill_manual(values=c("black","grey","white")))

grid.arrange(isl,cal,wood,sumi,tex,fla,legend, ncol=1)

#calc mean heterozygosity per sample and make df for plotting
gen.mat<-as.matrix(gen)
loci<-rowSums(is.na(gen.mat) == FALSE)
het<-rowSums(gen.mat == 1, na.rm = TRUE)/loci
het.df<-data.frame(id=locs$id,subspecies=locs$subspecies,species=locs$species,het=het)

calc.nuc.div<-function(vcfR, popmap){
  
  #make vcfR into gt matrix
  gt.matrix<-extract.gt(vcfR)
  #convert values for easier computation
  gt.matrix[gt.matrix=="0/0"]<-0
  gt.matrix[gt.matrix=="0/1"]<-1
  gt.matrix[gt.matrix=="1/1"]<-2
  gt.matrix<-as.data.frame(gt.matrix)
  for (i in 1:ncol(gt.matrix)){
    gt.matrix[,i]<-as.numeric(as.character(gt.matrix[,i]))
  }
  
  #open up vector to fill up with information by pop
  nucs<-c()
  #separate gt.matrix based on the population of interest
  for (i in 1:length(levels(popmap$pop))){
    target<-gt.matrix[,as.character(popmap$id[popmap$pop==levels(popmap$pop)[i]])]
    #remove sites where AF can't be calculated due to missing data leaving one or fewer called genotypes
    target<-target[rowSums(is.na(target)) < ncol(target)-1,]
    #calculate nucleotide diversity  within a given pop as proportion of pairwise differences between genotypes, per-SNP
    pi<-c() #open pi
      for (j in 1:nrow(target)){
        pi[j]<-sum(dist(target[j,][is.na(target[j,])==FALSE]) != 0)/length(dist(target[j,][is.na(target[j,])==FALSE]))
      }
  nucs[i]<-mean(pi)
  }
  df<-data.frame(pop=levels(popmap$pop),nuc.div=nucs)
  return(df)
}

nucs<-calc.nuc.div(vcfR, popmap)

#plot heterozygosity as violin plots for each subspecies
div.plot<-ggplot(het.df, aes(x=species, y=het)) + 
  #geom_violin(trim = FALSE)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 3, alpha=.75, aes(fill=species, color=species))+
  theme_classic()+
  scale_fill_manual(values=c("blue","pink","red","purple","orange","green"))+
  scale_color_manual(values=c("blue","pink","red","purple","orange","green"))+
  scale_x_discrete(labels=c("California","Florida","Island","Sumichrast","Texas","Woodhouse"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = .5, size = 14),
        legend.position = "none")+
  geom_point(nucs, mapping=aes(x=1:6, y=nuc.div/1.45), pch=8, cex=3)+
  labs(x="",y="heterozygosity")+
  scale_y_continuous(sec.axis = sec_axis(trans = (~.*1.45), name="nucleotide diversity", breaks = c(0,.05,.1)))+
  coord_flip()

lay <- rbind(c(NA,7),
             c(8,7),
             c(8,1),
             c(8,2),
             c(8,3),
             c(8,4),
             c(8,5),
             c(8,6),
             c(8,NA))
g2<-arrangeGrob(wood,tex,sumi,isl,fla,cal,legend,div.plot, layout_matrix = lay,
             widths =c(1,.35), heights=c(.1,.48,1,1,1,1,1,1,.48))


#make same figure and add compoplot
lay <- rbind(c(1,1),
             c(2,2),
             c(3,4))
g1<-arrangeGrob(map, com, df1, df2, layout_matrix = lay, heights=c(2,.5,1))


g<-grid.arrange(g1, g2, ncol=2, widths=c(1,.6))
#reorder compoplot and dotchart as Island, Cali, Woodhouse, Tex, Sumi, Florida
ggsave(filename = "~/Desktop/aph.data/Fig.1.pdf", plot = g, width = 11, height=8.5, units="in")
ggsave(filename = "~/Desktop/aph.data/Fig.1.png", plot = g, width = 11, height=8.5, units="in")

