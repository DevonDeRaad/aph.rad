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

ggplot(fsts, aes(WEIR_AND_COCKERHAM_FST, colour = comp)) +
  geom_freqpoly()+
  theme_classic()+
  xlab("Fst")+
  scale_color_discrete(name = "comparison", labels = c("californica/insularis",
                                                    "californica/woodhouseii",
                                                  "woodhouseii/sumichrasti",
                                                "woodhouseii/texana"))

ggplot(fsts, aes(WEIR_AND_COCKERHAM_FST, fill = comp)) +
  geom_histogram(bins=50)+
  theme_classic()

ggplot(fsts, aes(WEIR_AND_COCKERHAM_FST, after_stat(density), colour = comp)) +
  geom_freqpoly(bins=50)+
  theme_classic()


ca<-fsts[fsts$comp == "cali.wood",]
ca<- ca %>% drop_na()
table(ca$CHROM[ca$WEIR_AND_COCKERHAM_FST > .8])
length(ca$CHROM[ca$WEIR_AND_COCKERHAM_FST > .8])
table(ca$CHROM[ca$WEIR_AND_COCKERHAM_FST > .5])
length(ca$CHROM[ca$WEIR_AND_COCKERHAM_FST > .5])

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

#make BPcum
nCHR <- length(unique(ca$chrom))
ca$BPcum <- NA
s <- 0
nbp <- c()
for (i in c("Mt","1","1A","2","3","4","4A","5","6","7","8","9","10",
            "11","12","13","14","15","17","18","19","20","21",
            "23","24","26","27","Z")){
  nbp[i] <- max(ca[ca$chrom == i,]$POS)
  ca[ca$chrom == i,"BPcum"] <- ca[ca$chrom == i,"POS"] + s
  s <- s + nbp[i]
}

#check df
ca[1:10,]

#set plotting params
axis.set <- ca %>% 
  group_by(chrom) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
sig <- .01

#ggplot ca/wood
ggplot(ca, aes(x = BPcum, y =WEIR_AND_COCKERHAM_FST, 
                         color = chrom, size = WEIR_AND_COCKERHAM_FST)) +
  geom_point(alpha = 0.75) +
  #geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chrom, breaks = axis.set$center) +
  #scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  #scale_color_manual(values=c("grey","black","grey","black","grey","black","grey","black","grey","black","black","grey","grey",
  #                     "black","grey","black","grey","black","black","grey","black","grey","black",
  #                     "grey","black","grey","black","grey")) +
  scale_color_manual(values=c("black","darkgrey","black","darkgrey","black","darkgrey","black","darkgrey","black","darkgrey","darkgrey","black","black",
                              "darkgrey","black","darkgrey","black","darkgrey","darkgrey","black","darkgrey","black","darkgrey",
                              "black","darkgrey","black","darkgrey","black")) +
  
  scale_size_continuous(range = c(.2,1)) +
  labs(x = NULL, 
       y = "Fst") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))
#
#ggsave("ca.wood.manhattan.png", units="in", width=6, height=1.5, dpi=300)

ca$nuc<-NA
ca$nuc[ca$chrom != "Z"]<-"nuclear"
ca$nuc[ca$chrom == "Z"]<-"Z"
ggplot(ca, aes(WEIR_AND_COCKERHAM_FST, after_stat(density), colour = nuc)) +
  geom_freqpoly(bins=50)+
  theme_classic()

ggplot(ca, aes(WEIR_AND_COCKERHAM_FST, after_stat(count), colour = nuc)) +
  geom_freqpoly(bins=50)+
  theme_classic()

h = hist(ca$WEIR_AND_COCKERHAM_FST[ca$chrom != "Z"], plot=FALSE, breaks=25) # or hist(x,plot=FALSE) to avoid the plot of the histogram
h$density = h$counts/sum(h$counts)
g = hist(ca$WEIR_AND_COCKERHAM_FST[ca$chrom == "Z"], plot=FALSE, breaks=25) # or hist(x,plot=FALSE) to avoid the plot of the histogram
g$density = g$counts/sum(g$counts)
plot(g, col=rgb(1,0,0,1/4), freq=FALSE, ylim=c(0,.7))
plot(h, col=rgb(0,0,1,1/4), freq=FALSE, add=T)  # first histogram

png("Fst.hist.png", units="in", width=4, height=4, res=300)
plot(g, col=rgb(1,0,0,1/4), freq=FALSE, border=F, ylim=c(0,.65), ylab="proportion", xlab="Fst", main=NULL)
plot(h, col=rgb(0,0,1,1/4), freq=FALSE, border=F, add=T)  # first histogram
legend("topright", 
       legend = c("nuclear genome", "Z chromosome"), 
       col = c(rgb(0,0,1,1/4), 
               rgb(1,0,0,1/4)), 
       pch = 15, 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))
dev.off()

ca<-ca[is.na(ca$WEIR_AND_COCKERHAM_FST) == FALSE,]
fst.chrom<-aggregate(ca$WEIR_AND_COCKERHAM_FST, list(ca$chrom), mean)
table(ca$chrom)
#subset to only chroms > 10 SNPs
fst.chrom<-fst.chrom[c(1:14,19:26,28),]
fst.chrom$length<-c(118548696,20806668,21403021,21576510,16962381,16419078,14428146,
                    11648728,11201131,11587733,73657157,156412533,15652063,5979137,
                    112617285,69780378,20704505,62374962,36305782,39844632,27993427,
                    27241186,72861351)
plot(log(fst.chrom$length),fst.chrom$x, col=c(rep("black", times=nrow(fst.chrom)-1),"red"))
fst.chrom$gen<-c(rep("nuclear", times=nrow(fst.chrom)-1),"Z")
ggplot(fst.chrom, aes(x=log(length), y=x, color=gen))+
  geom_point(alpha=.3, cex=4)+
  scale_color_manual(values=c("black","red"))+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, lty=2)+
  theme_classic()+
  labs(x="log chromosome length", y="mean Fst")+
  theme(legend.title = element_blank(), text = element_text(size=14))+
  guides(color = guide_legend(override.aes = list(linetype = 0)))

table(ca$chrom[ca$WEIR_AND_COCKERHAM_FST > .8])

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
plot(car$recomb, car$WEIR_AND_COCKERHAM_FST, xlim=c(0,1))
points(car$recomb[car$chrom == "Z"],car$WEIR_AND_COCKERHAM_FST[car$chrom == "Z"], col="red")
cor.test(car$recomb, car$WEIR_AND_COCKERHAM_FST)
ggplot(car, aes(x=recomb, y=WEIR_AND_COCKERHAM_FST, col=nuc))+
  geom_point(alpha=.3, cex=2.5)+
  scale_color_manual(values = c("black","red"))+
  theme_classic()+
  geom_smooth(method=lm,  linetype="dashed",
              color="darkred", fill="grey")






#read in locality info for samples
locs<-read.csv("~/Desktop/aph.data/rad.sampling.csv")
#fix mislabeled sample
locs$id<-as.character(locs$id)
locs$id[locs$id == "A_californica_334171"]<-"A_woodhouseii_334171"

#subset locs to include only samples that passed filtering, and have been retained in the vcf
locs<-locs[locs$id %in% gen@ind.names,]
locs<-locs[c(1:19,21:56,20,57:95),]
locs[56,5]<-"woodhouseii"

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

caw<-sampling.df[c(1:13,18:23,26:32),]
#make full map
pac<-map_data("world")
#
ggplot()+
  geom_polygon(data = pac, aes(x=long, y = lat, group = group), fill="grey", col="black", cex=.1)+
  coord_sf(xlim = c(-123, -97), ylim = c(19, 45)) + 
  geom_point(data = caw, aes(x = long, y = lat, col=species, size=count), alpha =.5, show.legend=TRUE) +
  theme_void()+
  scale_color_manual(values=c("navy","forestgreen"))+
  scale_size_continuous(range = c(2,8))+
  guides(colour = guide_legend(override.aes = list(size = 4), order=1, label.theme = element_text(face = "italic")),
         size = guide_legend(nrow = 1, order = 2))+
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0.01, 0.01), legend.background=element_blank(), text = element_text(size=14))
#plot map colored by species

