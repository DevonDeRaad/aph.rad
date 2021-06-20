#!/bin/sh
#
#SBATCH --job-name=aph.gene.trees               # Job Name
#SBATCH --nodes=1             # 40 nodes
#SBATCH --ntasks-per-node=1               # 1 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/work/aph.rad/astral/bs10.7species     # Set working d$
#SBATCH --mem-per-cpu=10gb            # memory requested
#SBATCH --time=10000

#for i in /home/d669d153/work/aph.rad/astral/nex/*.nex
#do
#/Users/devder/Downloads/iqtree-1.6.12-mac/bin/iqtree -s ${i} -bb 1000 -wbt
#done

#The IQ-TREE analyses will have generated a number of files for each gene, containing run info, a distance matrix, starting trees, and so on.
#The only output files required for our further analyses are those ending in .treefile (the maximum-likelihood gene trees with branch lengths)
#and .ufboot (the set of bootstrap trees without branch lengths). To clean up the directory and keep only the important files,
#use the following commands: (run in terminal)
  
#rm nex/*.bionj
#rm nex/*.ckp.gz
#rm nex/*.contree
#rm nex/*.iqtree
#rm nex/*.log
#rm nex/*.mldist
#rm nex/*.model.gz
#rm nex/*.splits.nex
#rm nex/*.uniqueseq.phy

#Next, have a look at one of the files containing the maximum-likelihood trees, e.g. with the following command:
#less nex/locus.1.nex.treefile

#Since ASTRAL will require as input a single file containing all gene trees, combine all files with maximum-likelihood trees into a single file named ml_best.trees, using the following command:
#cat nex/*.treefile > ml_best.trees
#To further clean up the directory, you could then also remove all files that contain the maximum-likelihood trees for single genes, using
#rm nex/*.treefile

#We'll first use the set of bootstrapped trees to estimate node support on the species tree.
#To do so, ASTRAL requires as input a single file with the names of all the files containing bootstrapped trees. We can generate such a file with the following command:
#ls nex/*.ufboot > ml_boot.txt


#we need a namemap in order to assign tips to species in astral
#make namemap by hand

#load java
module load java

#We can then run ASTRAL with two input files: The file containing the maximum-likelihood trees for each gene (ml_best.trees), and the file containing the names of all files with bootstrapped trees (ml_boot.txt).
#The first of these is to be specified with ASTRAL's option -i, and the second should be given with option -b.
#In addition, we'll use option -o to set the name of the output file to species_boot.trees:
#java -jar /panfs/pfs.local/work/bi/bin/Astral_5.6.1/astral.5.6.1.jar -i ml_best.trees -a namemap.txt -b ml_boot.txt -o species_boot.trees

#Have a look at the output file species_boot.trees using a text editor (or again the command less). You'll see that it contains 102 lines. 
#The first 100 of these lines represent species trees in Newick format estimated for each the first 100 bootstrapped trees of each gene.
#On line 101 is a consensus tree for these 100 trees. Finally, the last line contains the species tree estimated from the maximum-likelihood gene trees,
#annotated with node support based on the bootstrap trees. Thus, the last line contains the species tree that we'll use for interpretation.
#However, before visualizing the species tree in FigTree, first conduct the second ASTRAL analysis based on the maximum-likelihood trees alone.
#Do so using the following command:
#java -jar /home/d669d153/work/Astral/astral.5.7.7.jar -i ml_best_BS10.trees -a namemap.7spec.txt -o species_7spec_BS10_pp.tre

#export quartet support for each internal branch using the '-t 16' command
#java -jar /home/d669d153/work/Astral/astral.5.7.7.jar -q species_7spec_BS10_pp.tre -i ml_best_BS10.trees -t 16 -o species_7spec_BS10_pp_quartets.tre -a namemap.7spec.txt

#run polytomy test using the option '-t 10'
java -jar /home/d669d153/work/Astral/astral.5.7.7.jar -q species_7spec_BS10_pp.tre -i ml_best_BS10.trees -t 10 -o species_7spec_BS10_pp_polytomy.tre -a namemap.7spec.txt

#The output file named species_pp.tre should contain a single species tree, annotated with posterior probabilities as node support.
  
#Since we are interested only in the last of the trees from file species_boot.trees as well as the tree from file species_pp.tre,
#we'll generate a new file named species.trees that contains both of these two trees using the following commands:
#tail -n 1 species_boot.trees > species.trees
#cat species_pp.tre >> species.trees


#now try contracting branches with extremely low support using newick utilities
#nw_ed  ml_best.trees 'i & b<=10' o > ml_best_BS10.trees

#re-run astral with low support branches concatenated
#java -jar /panfs/pfs.local/work/bi/bin/Astral_5.6.1/astral.5.6.1.jar -i ml_best_BS10.trees -a namemap.txt -o species_pp_BS10.tre

#java -jar /panfs/pfs.local/work/bi/bin/Astral_5.6.1/astral.5.6.1.jar -q species_pp_BS10.tre -i ml_best_BS10.trees -t 16 -a namemap.txt -o species_pp_BS10_quartets.tre

