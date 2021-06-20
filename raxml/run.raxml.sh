#!/bin/sh
#
#SBATCH --job-name=raxml.aph               # Job Name
#SBATCH --nodes=1             # 40 nodes
#SBATCH --ntasks-per-node=15               # 40 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/work/aph.rad/raxml     # Set working d$
#SBATCH --mem-per-cpu=1500            # memory requested
#SBATCH --time=10000

## this is an example call to run raxml tree inference w/ bootstrapping
#'-f a' Tell RAxML to conduct a rapid Bootstrap analysis and search for the best-scoring ML tree in one single program run
#'-T 5' run 5 separate threads (need at least 5 CPUs)
#'-m ASC_GTRGAMMA' GTR + Optimization of substitution rates + GAMMA model of rate heterogeneity (alpha parameter will be estimated). 
# The ASC prefix will correct the likelihood for ascertainment bias. You will also need to specify the correction type via ­­asc­corr!
#'-n all.samps' name output files
#'-o sac' option to set outgroup, not necessary
#'-s populations.var.phylip' specify phylip file for input
#/panfs/pfs.local/work/bi/bin/standard-RAxML-8.2.11/raxmlHPC-PTHREADS-SSE3 -f a -T 15 -m ASC_GTRGAMMA --asc-corr=lewis -# 1000 -x 12345 -p 54321 -n all.samps -s unlinked.filtered.phylip.txt

#use stacks to output the whitelisted loci from samples that passed filtering in phylip format
#/home/d669d153/work/stacks-2.41/populations -P /home/d669d153/work/aph.rad -O . -M phylip.popmap.txt --whitelist whitelist.txt --phylip-var-all

#run raxml on the concatenated phylip file
/panfs/pfs.local/work/bi/bin/standard-RAxML-8.2.11/raxmlHPC-PTHREADS-SSE3 -f a -T 15 -m GTRGAMMA -# 1000 -x 12345 -p 54321 -n full.locus -s populations.all.phylip






