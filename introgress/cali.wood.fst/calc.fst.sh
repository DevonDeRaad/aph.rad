#!/bin/sh
#
#SBATCH --job-name=aph.fst            # Job Name
#SBATCH --nodes=1             # 40 nodes
#SBATCH --ntasks-per-node=1               # 40 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/scratch/aph.rad/fst  # Set working d$
#SBATCH --mem-per-cpu=5gb            # memory requested
#SBATCH --time=10000



#
vcftools --vcfgz cal.wood.fst.vcf.gz --weir-fst-pop cali.txt --weir-fst-pop wood.txt --out cali.wood
