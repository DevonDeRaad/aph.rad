#!/bin/sh
#
#SBATCH --job-name=treemix              # Job Name
#SBATCH --nodes=1             # 40 nodes
#SBATCH --ntasks-per-node=1               # 40 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/work/aph.rad/treemix/7spec     # Set working d$
#SBATCH --mem-per-cpu=1000           # memory requested
#SBATCH --time=1000

#convert vcf into treemix file
/home/d669d153/work/stacks-2.41/populations --in_vcf unlinked.filtered.recode.vcf -O . --treemix -M 7spec.map.txt
#remove stacks header
echo "$(tail -n +2 unlinked.filtered.recode.p.treemix)" > unlinked.filtered.recode.p.treemix
#gzip file for input to treemix
gzip unlinked.filtered.recode.p.treemix

#run treemix with m0
/panfs/pfs.local/work/bi/bin/treemix-1.13/src/treemix -i unlinked.filtered.recode.p.treemix.gz -root coerulescens -o treem0

#add 1 migration edge
/panfs/pfs.local/work/bi/bin/treemix-1.13/src/treemix -i unlinked.filtered.recode.p.treemix.gz -m 1 -g treem0.vertices.gz treem0.edges.gz -o treem1

#add 2 migration edges
/panfs/pfs.local/work/bi/bin/treemix-1.13/src/treemix -i unlinked.filtered.recode.p.treemix.gz -m 1 -g treem1.vertices.gz treem1.edges.gz -o treem2

#add 3 migration edges
/panfs/pfs.local/work/bi/bin/treemix-1.13/src/treemix -i unlinked.filtered.recode.p.treemix.gz -m 1 -g treem2.vertices.gz treem2.edges.gz -o treem3

 
