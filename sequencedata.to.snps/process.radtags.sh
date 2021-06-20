#!/bin/sh
#
#SBATCH --job-name=process.radtags            # Job Name
#SBATCH --nodes=1              # 40 nodes
#SBATCH --ntasks-per-node=1             # 40 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/scratch/plate1/      # Set working d$
#SBATCH --mem-per-cpu=5gb            # memory requested
#SBATCH --time=5000

/home/d669d153/work/stacks-2.3b/process_radtags -p .  -o . -b plate.1.barcodes.txt -e ndeI -r -c -q
