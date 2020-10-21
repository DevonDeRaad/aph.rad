library(tidyr)
#Fst comparisons
setwd("~/Desktop/aph.data/")
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
ggsave("ca.wood.manhattan.png", units="in", width=6, height=1.5, dpi=300)

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
  theme(legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(linetype = 0)))


