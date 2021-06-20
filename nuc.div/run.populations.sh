#!/bin/sh
#
#SBATCH --job-name=combine.g.vcf               # Job Name
#SBATCH --nodes=1              # 40 nodes
#SBATCH --ntasks-per-node=1               # 40 CPU allocation per Task
#SBATCH --partition=sixhour            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/work/aph.rad/nuc.div/indiv.samps/        # Set working d$
#SBATCH --mem-per-cpu=350000            # memory requested
#SBATCH --time=360


/home/d669d153/work/stacks-2.41/populations -P /home/d669d153/work/aph.rad/ -M indiv.popmap.txt -O .
