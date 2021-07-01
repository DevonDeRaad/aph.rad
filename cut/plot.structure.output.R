#plot structure output
library(pophelper)
library(gridExtra)

#step 1: move into the structure output directory
setwd("~/Desktop/aph.data/structure/")

sfiles1 <- list.files(path=system.file("files/structure-ci",package="pophelper"),full.names=TRUE)
sfiles1 <- list.files(path=getwd(),full.names=TRUE)
slist1 <- readQ(files=sfiles1, readci=TRUE)
head(tabulateQ(slist))


sfiles1 <- list.files(path=getwd(),full.names=TRUE)[c(seq(from=2, to=72, by=2),76,78)]
slist <- readQ(sfiles1)
tabulateQ(slist)
sr1 <- summariseQ(tabulateQ(slist))
p <- evannoMethodStructure(data=sr1,exportplot=F,returnplot=T,returndata=F,basesize=12,linesize=0.7)
grid.arrange(p)
evannoMethodStructure(data=sr1,exportplot=F,returnplot=T,returndata=T,basesize=12,linesize=0.7)

plotQ(slist,imgoutput="join",returnplot=T,exportplot=F,quiet=T,basesize=11)
plotQ(slist[1],returnplot=T,exportplot=F,quiet=T,basesize=11, showindlab = T)


#create whitelist for raxml
vcfR <- read.vcfR("~/Desktop/aph.data/unlinked.filtered.recode.vcf")
locus.list<-vcfR@fix[,3]
locus.list<-gsub(":.*","",locus.list)
length(locus.list)
write.table(locus.list, "~/Downloads/whitelist.txt", quote = F, row.names = F, col.names = F)


