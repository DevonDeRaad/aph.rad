
#workflow describing the arduous process of generating individual gene trees when you did your genotype filtering in R
#meaning you have to go all the way back to gstacks to get your whole locus sequences

#read in vcfR with 2725 filtered, unlinked loci that we want to use for gene tree estimation
vcfR<-read.vcfR("~/Desktop/aph.data/unlinked.filtered.recode.vcf")
head(vcfR@fix)
vcfR@fix[,3]

#get list of just locus names to use as whitelist for stacks
whitelist<-sub(":.*", "", vcfR@fix[,3])
whitelist

#write out whitelist for stacks
write.table(whitelist, file = "~/Downloads/2725.whitelist.txt", quote = F, row.names = F, col.names = F)

#run this in terminal in the directory where you ran gstacks (whitelist includes the loci you want and popmap includes the samples you want based on filtering)
#/home/d669d153/work/stacks-2.41/populations -P /home/d669d153/work/aph.rad -O . -M phylip.popmap.txt --whitelist whitelist.txt --phylip-var-all

#####
#read in the phylip with all loci concatenated
dna<-read.dna(file = "~/Desktop/aph.data/astral/populations.all.phylip")
as.list(dna)
dim(dna)
#write.FASTA(x = dna, file="~/Desktop/aph.data/svdquartets/whole.locus.2725.fasta")

#read in the partition file
parts<-read.table("~/Desktop/aph.data/astral/populations.all.partitions.phylip")
parts$parts<-sub("=.*","",parts$V2)
parts$begin<-sub(".*=","",parts$V2)
parts$begin<-sub("-.*","",parts$begin)
parts$end<-sub(".*=","",parts$V2)
parts$end<-sub(".*-","",parts$end)

#write a loop that generates a separate fasta for each partition
for (i in 1:nrow(parts)){
  locus<- dna[,c(parts$begin[i]:parts$end[i])]
  write.FASTA(x = locus, file=paste0("~/Desktop/aph.data/astral/fastas/locus.",i,".fasta"))
  print(i)
}

#following steps come from: https://github.com/mmatschiner/tutorials/blob/master/ml_species_tree_inference/README.md#astral
#To remove sequences that contain only missing information from all alignments,
#and at the same time translate all alignments into Nexus format,
#we can use the Python script convert.py with the following commands: (run in terminal)
#mkdir nex
#for i in fastas/*.fasta
#do
#gene_id=`basename ${i%.fasta}`
#python3 convert.py ${i} nex/${gene_id}.nex -f nexus -m 0.9
#done

#run this in the directory to make sure that missing data is correctly recognized as being encoded by 'N' rather than '?'
#for i in *
#  do
#sed -i '' 's/ missing=?/ missing=N/' $i
#done

#We can now use IQ-TREE to generate maximum-likelihood gene trees for all alignments. To be able to later use bootstrapping with ASTRAL,
#we will also generate sets of trees from bootstrapped gene alignments in the same IQ-TREE analysis, by using the option -B 
#we do not specify a substition model with the -m option and therefore allow IQ-TREE to automatically select the best-fitting model.
#we will now ensure that the bootstrap trees are actually written to a file, and not just summarized in the output, by specifying the option --wbt.
#To start IQ-TREE in a loop so that it analyzes one gene alignment after the other, use the following commands: (run on the cluster, 2725 gene trees take ~8 hours to generate on my laptop)
#SWITCH TO CLUSTER AT THIS POINT

#!/bin/sh
#
#SBATCH --job-name=aph.gene.trees               # Job Name
#SBATCH --nodes=1             # 40 nodes
#SBATCH --ntasks-per-node=1               # 1 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/work/aph.rad/astral     # Set working d$
#SBATCH --mem-per-cpu=1gb            # memory requested
#SBATCH --time=10000

#for i in /home/d669d153/work/aph.rad/astral/nex/*.nex
#do
#/Users/devder/Downloads/iqtree-1.6.12-mac/bin/iqtree -s ${i} -bb 1000 -wbt
#done

#The IQ-TREE analyses will have generated a number of files for each gene, containing run info, a distance matrix, starting trees, and so on.
#The only output files required for our further analyses are those ending in .treefile (the maximum-likelihood gene trees with branch lengths)
#and .ufboot (the set of bootstrap trees without branch lengths). To clean up the directory and keep only the important files,
#use the following commands: (run in terminal)
  
rm nex/*.bionj
rm nex/*.ckp.gz
rm nex/*.contree
rm nex/*.iqtree
rm nex/*.log
rm nex/*.mldist
rm nex/*.model.gz
rm nex/*.splits.nex
rm nex/*.uniqueseq.phy

#Next, have a look at one of the files containing the maximum-likelihood trees, e.g. with the following command:
less nex/locus.1.nex.treefile

#Since ASTRAL will require as input a single file containing all gene trees, combine all files with maximum-likelihood trees into a single file named ml_best.trees, using the following command:
cat nex/*.treefile > ml_best.trees
#To further clean up the directory, you could then also remove all files that contain the maximum-likelihood trees for single genes, using
#rm nex/*.treefile

#We'll first use the set of bootstrapped trees to estimate node support on the species tree.
#To do so, ASTRAL requires as input a single file with the names of all the files containing bootstrapped trees. We can generate such a file with the following command:
ls nex/*.ufboot > ml_boot.txt


#we need a namemap in order to assign tips to species in astral
#easiest to just do by hand unfortunately as the required format is very un-table-like
#save as namemap.txt

#We can then run ASTRAL with two input files: The file containing the maximum-likelihood trees for each gene (ml_best.trees), and the file containing the names of all files with bootstrapped trees (ml_boot.txt).
#The first of these is to be specified with ASTRAL's option -i, and the second should be given with option -b.
#In addition, we'll use option -o to set the name of the output file to species_boot.trees:
module load java
java -jar /home/d669d153/work/Astral/astral.5.7.7.jar -i ml_best.trees -a namemap.txt -b ml_boot.txt -o species_boot.trees

#Have a look at the output file species_boot.trees using a text editor (or again the command less). You'll see that it contains 102 lines. 
#The first 100 of these lines represent species trees in Newick format estimated for each the first 100 bootstrapped trees of each gene.
#On line 101 is a consensus tree for these 100 trees. Finally, the last line contains the species tree estimated from the maximum-likelihood gene trees,
#annotated with node support based on the bootstrap trees. Thus, the last line contains the species tree that we'll use for interpretation.
#However, before visualizing the species tree in FigTree, first conduct the second ASTRAL analysis based on the maximum-likelihood trees alone.
#Do so using the following command:
java -jar /home/d669d153/work/Astral/astral.5.7.7.jar -i ml_best.trees -a namemap.txt -o species_pp.tre

java -jar /home/d669d153/work/Astral/astral.5.7.7.jar -q species_pp.tre -i ml_best.trees -t 16 -a namemap.txt -o species_pp_quartets.tre

#The output file named species_pp.tre should contain a single species tree, annotated with posterior probabilities as node support.
  
#Since we are interested only in the last of the trees from file species_boot.trees as well as the tree from file species_pp.tre,
#we'll generate a new file named species.trees that contains both of these two trees using the following commands:
tail -n 1 species_boot.trees > species.trees
cat species_pp.tre >> species.trees


#now try contracting branches with extremely low support using newick utilities
/panfs/pfs.local/work/bi/bin/newick_utils/bin/nw_ed  ml_best.trees 'i & b<=10' o > ml_best_BS10.trees

#re-run astral
java -jar /home/d669d153/work/Astral/astral.5.7.7.jar -i ml_best_BS10.trees -a namemap.txt -o species_pp_BS10.tre

java -jar /home/d669d153/work/Astral/astral.5.7.7.jar -q species_pp_BS10.tre -i ml_best_BS10.trees -t 16 -a namemap.txt -o species_pp_BS10_quartets.tre






