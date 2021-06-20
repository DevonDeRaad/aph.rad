#!/bin/sh
#
#SBATCH --job-name=aph.treemix               # Job Name
#SBATCH --nodes=1             # 40 nodes
#SBATCH --ntasks-per-node=1               # 40 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/work/aph.rad/treemix/7spec/boots     # Set working d$
#SBATCH --mem-per-cpu=1000            # memory requested
#SBATCH --time=2000

#bootstraps over 200s snps with migration edges
for i in {1..100}; do
    /panfs/pfs.local/work/bi/bin/treemix-1.13/src/treemix -i /home/d669d153/work/aph.rad/treemix/7spec/unlinked.filtered.recode.p.treemix.gz -m 1 -g /home/d669d153/work/aph.rad/treemix/7spec/treem1.vertices.gz /home/d669d153/work/aph.rad/treemix/7spec/treem1.edges.gz -bootstrap -k 100 -o $i.treemix
done;

# unzip the tree files
for i in { ls *treeout.gz }; do
    gzip -d $i
done;

# in R:
#summarize bootstraps
module load R
R --no-save
setwd("/home/d669d153/work/aph.rad/treemix/7spec/boots/")

x <- list.files(pattern="*treeout")
for(a in 1:length(x)) {
	if (a==1) {
		output <- scan(x[a], what="character")[1]
	} else {
		output <- c(output, scan(x[a], what="character")[1])
	}
}
write(output, file="aph.100.bootstraps.trees", ncolumns=1)
quit()

# in bash 
# summarize bootstraps
/home/d669d153/work/DendroPy/applications/sumtrees/sumtrees.py --output=aph.boots.summed.tre --min-clade-freq=0.05 aph.100.bootstraps.trees
