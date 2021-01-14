
#make map.txt for converting vcf into treemix file
vcfR <- read.vcfR("~/Desktop/aph.data/unzipped.filtered.vcf")
samps<-colnames(vcfR@gt)[2:96]
pops<-gsub("A_","",samps)
popmap<-data.frame(id=samps,
                   raxml.spec=gsub("_.*","",pops))
write.table(popmap, "~/Downloads/map.txt", quote = F, row.names = F, col.names = F, sep = '\t')


#script to plot treemix output
#step 1: copy the entire treemix outdirectory from the KU cluster to my local machine
#scp -r d669d153@hpc.crc.ku.edu:/home/d669d153/work/aph.rad/treemix /Users/devder/Desktop/aph.data/

#step 2: source plotting functions that are distributed with treemix
source("~/Downloads/plotting_funcs.R")

#step 3: move into the treemix output directory and plot trees
setwd("~/Desktop/aph.data/treemix/")

#0 edge
plot_tree("treem0")
#1 edge
plot_tree("treem1", plus = 0.02, arrow=.1, ybar = 0.3, scale=F, lwd=1.2)
#2 edges
plot_tree("treem2")
#3 edges
plot_tree("treem3")


#plot to see how much variance is explained by each edge
#.994->m0
m=NULL
for(i in 0:3){
  m[i+1] <- get_f(paste0("treem",i))
}

m
plot(seq(0,3),m,pch="*",cex=2,col="blue", type="b",xlab="migration edge number", ylab="% explained variance")

pdf(file="tree.mix.pdf", width=4.5, height = 4)
plot_tree("treem1", plus = 0.02, arrow=.1, ybar = 0.3, scale=F, lwd=1.2)
dev.off()

pdf(file="tree.mix.suppmat.pdf", width=4, height = 4)
plot_tree("treem2", cex = .7, ybar = F)
dev.off()

pdf(file="variance.explained.pdf", width=4, height = 4)
plot(seq(0,3),m,pch="*",cex=2,col="blue", type="b",xlab="migration edge number", ylab="% explained variance")
dev.off()






pop.uniq <- c("coerulescens","sumichrasti","texana","woodhouseii","insularis","californica")
write.table(pop.uniq, "poporder", sep = " ", quote = F, row.names = F, col.names = F)

plot_resid("tree.1.2", "poporder")
plot_resid("treemix.1", "poporder")
plot_resid("treemix.2", "poporder")

install.packages("SiZer")
install.packages("OptM")
library(OptM)
test.optM = optM("~/Desktop/aph.data/treemix")
#Finally, plot the results:
plot_optM(test.optM, method = "Evanno")
folder <- system.file("~/Desktop/aph.data/treemix/", package = "OptM")
test.linear = optM(folder, method = "linear")
plot_optM(test.linear, method = "linear")


