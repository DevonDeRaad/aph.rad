
library(ggplot2)
#plot admixture results

#setwd to admixture directory run on the cluster
setwd("~/Downloads/admixture/")

#read in log error values to determine optimal K
log<-read.table("log.errors.txt")[,c(3:4)]
#use double backslash to interpret the opening parentheses literally in the regular expression
log$V3<-gsub("\\(K=", "", log$V3)
log$V3<-gsub("):", "", log$V3)
#interpret K values as numerical
log$V3<-as.numeric(log$V3)
#rename columns
colnames(log)<-c("Kvalue","cross.validation.error")

#make plot showing the cross validation error across K values 1:10
ggplot(data=log, aes(x=Kvalue, y=cross.validation.error, group=1)) +
  geom_line(linetype = "dashed")+
  geom_point()+
  ylab("cross-validation error")+
  xlab("K")+
  scale_x_continuous(breaks = c(1:10))+
  theme_classic()

#save file
#ggsave(filename = "~/Desktop/aph.data/admixture/cross.validation.error.pdf", height = 3, width = 4)

#read in input file in order to get list of input samples in order
samps<-read.table("binary_fileset.fam")[,1]

#read in all ten runs and save each dataframe in a list
runs<-list()
#read in log files
for (i in 1:10){
  runs[[i]]<-read.table(paste0("binary_fileset.", i, ".Q"))
}

dev.off()
par(mfrow=c(5,1))
#plot each run
for (i in 1:5){
barplot(t(as.matrix(runs[[i]])), col=rainbow(i), ylab="Ancestry", border="black")
}

#plot each run
for (i in 6:10){
  barplot(t(as.matrix(runs[[i]])), col=rainbow(i), ylab="Ancestry", border="black")
}

#define function to get ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#reorder samples to reflect clustering at K=6
runs[[5]]<-runs[[5]][c(33:46,1:19,21,23:32,22,20,57:73,81:95,74:80,47:56),]
runs[[6]]<-runs[[6]][c(33:46,1:19,21,23:32,22,20,57:73,81:95,74:80,47:56),]

#plot barplots for the two most relevant runs
par(mfrow=c(2,1))
for (i in 5:6){
  barplot(t(as.matrix(runs[[i]])), col=gg_color_hue(i), ylab="Ancestry", border="black")
}

#FMNH 342048 is the only sample with > .1 mis-assignment at K=5
runs[[5]]



#Read in vcf and subset to only woodhouse to do zoomed in admixture run
vcfR <- read.vcfR("~/Desktop/aph.data/unlinked.filtered.recode.vcf")
colnames(vcfR@gt)
vcfR@fix[1:5,1:5]
#subset to only woodhouse samples
vcfR@gt <- vcfR@gt[,c(1,21,58:96)]

#remove samples with no variable sites left
source("~/Desktop/snpfiltR/min_mac.R")
vcf<-min_mac(vcfR = vcfR, min.mac = 1)
vcf
#write out for admixture analysis on the cluster
#write.vcf(vcf, file="~/Downloads/woodhouse.only.vcf.gz")



#####
#investigate results from only woodhouse's scrub-jay (sumichrasti aslo removed)
#setwd to admixture directory from run on the cluster
setwd("~/Downloads/woodhouse.admixture/")

#read in log error values to determine optimal K
log<-read.table("log.errors.txt")[,c(3:4)]
#use double backslash to interpret the opening parentheses literally in the regular expression
log$V3<-gsub("\\(K=", "", log$V3)
log$V3<-gsub("):", "", log$V3)
#interpret K values as numerical
log$V3<-as.numeric(log$V3)
#rename columns
colnames(log)<-c("Kvalue","cross.validation.error")

#make plot showing the cross validation error across K values 1:10
ggplot(data=log, aes(x=Kvalue, y=cross.validation.error, group=1)) +
  geom_line(linetype = "dashed")+
  geom_point()+
  ylab("cross-validation error")+
  xlab("K")+
  scale_x_continuous(breaks = c(1:10))+
  theme_classic()

#read in all ten runs and save each dataframe in a list
wood.runs<-list()
#read in log files
for (i in 1:10){
  wood.runs[[i]]<-read.table(paste0("binary_fileset.", i, ".Q"))
}

#plot
par(mfrow=c(5,1))
#plot each run
for (i in 1:5){
  barplot(t(as.matrix(wood.runs[[i]])), col=rainbow(i), ylab="Ancestry", border="black")
}

#plot each run
for (i in 6:10){
  barplot(t(as.matrix(wood.runs[[i]])), col=rainbow(i), ylab="Ancestry", border="black")
}

samps<-read.table("binary_fileset.fam")[,1]

#reorder samples to reflect clustering at K=6
wood.runs[[2]]<-wood.runs[[2]][c(1:18,26:40,19:25),]
wood.runs[[3]]<-wood.runs[[3]][c(1:18,26:40,19:25),]

#save all plots together
pdf("~/Desktop/aph.data/admixture/all.admix.plots.pdf", width = 8, height=6)
#set number of rows/columns
par(mfrow=c(4,1))
#set margins
par(mar = c(3, 3, 0, 0), oma = c(1, 1, 1, 1))
#set line width
#opar <- par(lwd = 1)
for (i in 5:6){
  barplot(t(as.matrix(runs[[i]])), col=gg_color_hue(i), ylab="Ancestry", border="black")
}
for (i in 2:3){
  barplot(t(as.matrix(wood.runs[[i]])), col=gg_color_hue(i), ylab="Ancestry", border="black")
}
dev.off()
