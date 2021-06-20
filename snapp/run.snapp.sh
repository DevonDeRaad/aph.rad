#!/bin/sh
#
#SBATCH --job-name=snapp               # Job Name
#SBATCH --nodes=1             # 40 nodes
#SBATCH --ntasks-per-node=15               # 40 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/work/aph.rad/snapp/7spec    # Set working d$
#SBATCH --mem-per-cpu=800            # memory requested
#SBATCH --array=1-5
#SBATCH --time=20000

/home/d669d153/work/beast/bin/beast -threads 15 rep$SLURM_ARRAY_TASK_ID.xml
