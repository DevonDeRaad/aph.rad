###############################################
###############################################
## script adapted from Derkarabetian S., Castillo S., Peter K.K., Ovchinnikov S., Hedin M. "An Empirical Demonstration of Unsupervised Machine Learning in Species Delimitation"
###############################################
###############################################
#required packages
library("adegenet")
library("randomForest")
library("PCDimension")
library("mclust")
library("cluster")
library("MASS")
library("factoextra")
library("tsne")

# import str file. Adjust input file name, n.ind, and n.loc for specific file/dataset.
# example dataset used in this study
# data <- read.structure("Metano_UCE_SNPs_70percent_random_structure-adegenet.str", n.ind=30, n.loc=316, onerowperind=FALSE, col.lab=1, col.pop=3, col.others=NULL, row.marknames=NULL, NA.char="-9", pop=NULL, ask=FALSE, quiet=FALSE)
# data <- read.structure("input.str", n.ind=XX, n.loc=XX, onerowperind=FALSE, col.lab=1, col.pop=0, col.others=NULL, row.marknames=NULL, NA.char="-9", pop=NULL, ask=FALSE, quiet=FALSE)

#read in vcf as vcfR
#vcfR <- read.vcfR("~/Desktop/aph.data/unlinked.filtered.recode.vcf")
vcfR <- read.vcfR("~/Desktop/aph.data/unzipped.filtered.vcf")
dim(vcfR@gt) #1779 out of 16307 SNPs contain no missing data
#filter out SNPs with missing data
#filter to only Z chrom SNPs
e<-vcfR2genlight(vcfR)
e<-as.matrix(e)
dim(e)
e<-e[,vcfR@fix[,1] == "PseudochrZ"]
#calculate Z chrom heterozygosity for all samples
f<-c()
for (i in 1:nrow(e)){
  f[i]<-sum(e[i,] == 1, na.rm=T)
}
hist(f, nclass = 20)
table(f)
#vcfR@fix<-vcfR@fix[rowSums(is.na(vcfR@gt)) == 0,]
#vcfR@gt<-vcfR@gt[rowSums(is.na(vcfR@gt)) == 0,]

#convert vcfR into a 'genind' object
data<-vcfR2genind(vcfR)

#scale the genind
data_scaled <- scaleGen(data, center=FALSE, scale=FALSE, NA.method=c("mean"), nf)
data_scaled <- scaleGen(data, center=FALSE, scale=FALSE)

###############################################
###############################################
# PCA and DAPC
###############################################
###############################################

# DAPC (interactive, requires input)
#clusters <- find.clusters(data, max.n.clust=10, n.iter=1e6, n.start=10)
#with the appropriate settings
clusters <- find.clusters(data, max.n.clust=10, n.iter=1e6, n.start=10, n.pca = 50, n.clust = 6)
#results <- dapc(data, clusters$grp, perc.pca = NULL)
#with appropriate settings
results <- dapc(data, clusters$grp, perc.pca = NULL, n.pca = 6, n.da = 4)
compoplot(results)
scatter.dapc(results, xax = 1, yax=2)
scatter.dapc(results, xax = 3, yax=4)
dap<-results$tab
dap$clusters<-clusters$grp
ggplot(data=dap, aes(x=`PCA-pc.3`, y=`PCA-pc.4`, col=clusters))+
  geom_point(cex=3)+
  theme_classic()
#prefers 6 groups with texas as a unique cluster

grp_k <- nlevels(clusters$grp)

# PCA, can adjust nf to include more components
pca1 <- dudi.pca(data_scaled, center=TRUE, scale=TRUE, scannf=FALSE, nf=5)
# PCA with DAPC groups
pc<-pca1$li
ggplot(data=pc, aes(x=Axis3, y=Axis4, col=dap$clusters))+
  geom_point(cex=3)+
  theme_classic()

# pam clustering on pca output
for (i in 2:10){
  print(paste(i, mean(silhouette(pam(pc, i))[, "sil_width"])))
}
#prefers 6 groups with non-perfect split between US and Mexico

pam(pc, 6)
# determine optimal k of PCA via hierarchical clustering with BIC
# adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
pca_clust <- Mclust(pc, G=1:10)
pca_clust$classification
#prefers 6 groups, with non-perfect split between US and Mexico


###############################################
###############################################
# into the Random Forest, unsupervised
###############################################
###############################################

# convert genind scaled data to factors for randomForest
data_conv <- as.data.frame(data_scaled)
data_conv[is.na(data_conv)] <- ""
data_conv[sapply(data_conv, is.integer)] <- lapply(data_conv[sapply(data_conv, is.integer)], as.factor)
data_conv[sapply(data_conv, is.character)] <- lapply(data_conv[sapply(data_conv, is.character)], as.factor)
nsamp <- nrow(data_conv)

# unsupervised random forest
rftest <- randomForest(data_conv, ntree=5000)
#rftest <- randomForest(pca1$tab, ntree=500)
#rftest <- randomForest(data_scaled, ntree=500)

###############
# classic MDS
###############
# cMDS with optimal number of components to retain using broken-stick
cmdsplot1 <- MDSplot(rf=rftest, fac=results$grp, k=10) # may need to adjust number of dimensions if given error
cmdsplot_bstick <- PCDimension::bsDimension(cmdsplot1$eig)
cmdsplot2 <- MDSplot(rftest, results$grp, cmdsplot_bstick)

#cMDS plot with dapc groups
cmds<-as.data.frame(cmdsplot2$points)
ggplot(data=cmds, aes(x=`Dim 1`, y=`Dim 2`, col=dap$clusters))+
  geom_point(cex=3)+
  theme_classic()

# pam clustering on cMDS output
for (i in 2:10){
  print(paste(i, mean(silhouette(pam(cmdsplot1$points, i))[, "sil_width"])))
}
#prefers 6 groups, matching dapc
DAPC_pam_clust_prox <- pam(cmdsplot1$points, 6)
DAPC_pam_clust_prox$clustering

# cMDS with optimal k and clusters via PAM
cmds$clusters<-as.factor(DAPC_pam_clust_prox$clustering)
ggplot(data=cmds, aes(x=`Dim 1`, y=`Dim 2`, col=clusters))+
  geom_point(cex=3)+
  theme_classic()

# determine optimal k from cMDS via hierarchical clustering with BIC
# adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
cmdsplot_clust <- Mclust(cmdsplot2$points)
cmdsplot_clust$classification
#hierarchical clustering of random forest identifies 7 groups 

# cMDS with optimal k and clusters of RF via hierarchical clustering
cmds$clusters<-as.factor(cmdsplot_clust$classification)
ggplot(data=cmds, aes(x=`Dim 1`, y=`Dim 2`, col=clusters))+
  geom_point(cex=3)+
  theme_classic()

###############
# isotonic MDS
###############

# isoMDS
isomdsplot <- isoMDS(1-rftest$proximity)
# "The output of cmdscale on 1 - rf$proximity is returned invisibly" (MDSplot documentation)
#plot isomds with dapc groups
df<-as.data.frame(isomdsplot$points)
ggplot(data=df, aes(x=V1, y=V2, col=results$grp))+
  geom_point(cex=3)+
  theme_classic()

# pam clustering on isomds with optimal k from DAPC
for (i in 2:10){
  print(paste(i, mean(silhouette(pam(isomdsplot$points, i))[, "sil_width"])))
}
#pam prefers only 2 groups, Florida and everything else

# determine optimal k of RF via hierarchical clustering with BIC
# adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
isomdsplot_clust <- Mclust(isomdsplot$points, G =1:10)
isomdsplot_clust$classification
#prefers only 2 groups

# isoMDS with optimal k and clusters of RF via hierarchical clustering
ggplot(data=df, aes(x=V1, y=V2, col=as.factor(isomdsplot_clust$classification)))+
  geom_point(cex=3)+
  theme_classic()


###############################################
###############################################
# t-SNE
###############################################
###############################################

# prepare plot labels and such
# this makes it so it is grouped by DAPC clusters
colors = rainbow(length(unique(results$grp)))
names(colors) = unique(results$grp)
ecb = function(x,y){plot(x,t='n'); text(x, labels=results$grp, col=colors[results$grp])}

# t-SNE on principal components of scaled data
# adjust perplexity, initial_dims
# can do k=3 for 3D plot
# should do only <50 variables
# can do it on pca$li (if you reduce the number of components), or on cmdsplot2$points
tsne_p5 = tsne(pca1$tab, epoch_callback=ecb, max_iter=5000, perplexity=5, initial_dims=5)

# tSNE plot with DAPC groups
plot(tsne_p5, main="t-SNE perplexity=5 with DAPC optimal k and clusters", col=results$grp, pch=16)

# pam clustering with optimal k from DAPC
for (i in 2:10){
  print(paste(i, mean(silhouette(pam(tsne_p5, i))[, "sil_width"])))
}
#pam prefers same 6 groups as DAPC
pam(tsne_p5, 6)

# determine optimal k of tSNE via hierarchical clustering with BIC
# adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
tsne_p5_clust <- Mclust(tsne_p5)
mclust_grps_tsne_p5 <- as.numeric(tsne_p5_clust$classification)
max(mclust_grps_tsne_p5)
# t-SNE p5 with optimal k and clusters of RF via hierarchical clustering
plot(tsne_p5, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="t-SNE p5 RF optimal K and clusters (hierarchical clustering)", col=mclust_grps_tsne_p5, pch=16)
mclust_grps_tsne_p5
f<-as.data.frame(tsne_p5)
# tSNE with optimal k and clusters via hierarchical clustering
ggplot(data=f, aes(x=V1, y=V2, col=as.factor(mclust_grps_tsne_p5)))+
  geom_point(cex=3)+
  theme_classic()

cbind(rownames(pca1$tab), mclust_grps_tsne_p5)
#prefers 6 groups where woodhouseii is split into US and Mexico groups, with Texas lumped with the US


#bring in results from VAE and visualize here:


