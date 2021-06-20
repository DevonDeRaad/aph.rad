#!/bin/sh
#
#SBATCH --job-name=radseq.aphelocoma               # Job Name
#SBATCH --nodes=1             # 40 nodes
#SBATCH --ntasks-per-node=1               # 40 CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/scratch/aph.rad  # Set working d$
#SBATCH --mem-per-cpu=20gb            # memory requested
#SBATCH --time=10000

files="A_californica_333849
A_californica_333854
A_californica_333855
A_californica_333857
A_californica_333860
A_californica_333862
A_californica_333868
A_californica_333871
A_californica_333873
A_californica_333874
A_californica_333914
A_californica_333917
A_californica_333925
A_californica_333931
A_californica_333932
A_californica_333934
A_californica_334002
A_californica_334012
A_californica_334015
A_californica_334017
A_californica_334019
A_californica_334171
A_californica_342037
A_californica_342048
A_californica_342050
A_californica_342051
A_californica_342066
A_californica_342072
A_californica_343428
A_californica_343430
A_californica_343432
A_californica_343438
A_californica_343442
A_californica_343451
A_californica_393615
A_californica_393616
A_californica_393721
A_coerulescens_396251
A_coerulescens_396254
A_coerulescens_396256
A_coerulescens_396259
A_coerulescens_396262
A_coerulescens_396263
A_coerulescens_396264
A_coerulescens_396265
A_insularis_334025
A_insularis_334029
A_insularis_334031
A_insularis_334032
A_insularis_334033
A_insularis_334034
A_insularis_334037
A_insularis_334038
A_sumichrasti_343502
A_sumichrasti_343503
A_sumichrasti_343510
A_sumichrasti_343511
A_sumichrasti_343512
A_sumichrasti_343513
A_sumichrasti_343514
A_sumichrasti_343515
A_sumichrasti_393633
A_sumichrasti_393635
A_sumichrasti_393636
A_sumichrasti_393637
A_sumichrasti_393638
A_sumichrasti_393639
A_sumichrasti_393640
A_woodhouseii_334059
A_woodhouseii_334062
A_woodhouseii_334063
A_woodhouseii_334086
A_woodhouseii_334088
A_woodhouseii_334096
A_woodhouseii_334132
A_woodhouseii_334133
A_woodhouseii_334134
A_woodhouseii_334142
A_woodhouseii_334148
A_woodhouseii_334153
A_woodhouseii_334156
A_woodhouseii_334161
A_woodhouseii_334170
A_woodhouseii_334172
A_woodhouseii_334188
A_woodhouseii_334190
A_woodhouseii_334196
A_woodhouseii_334210
A_woodhouseii_334211
A_woodhouseii_334217
A_woodhouseii_334230
A_woodhouseii_334240
A_woodhouseii_334241
A_woodhouseii_334242
A_woodhouseii_334243
A_woodhouseii_334244
A_woodhouseii_334247
A_woodhouseii_334250
A_woodhouseii_343453
A_woodhouseii_343458
A_woodhouseii_343461
A_woodhouseii_343476
A_woodhouseii_343480
A_woodhouseii_343481
A_woodhouseii_343483
A_woodhouseii_343497
A_woodhouseii_393605
A_woodhouseii_393606
A_woodhouseii_393697
A_woodhouseii_393698
A_woodhouseii_393699
A_woodhouseii_393702
A_woodhouseii_393712
A_woodhouseii_393713
A_woodhouseii_395768"

#
src=/home/d669d153/scratch/aph.rad #directory that contains a folder with sample data in fasta format

#index ref
/panfs/pfs.local/work/bi/bin/bwa/bwa index pseudochromosomes.fasta

#Align paired-end data with BWA, convert to BAM and SORT.
for sample in $files
do 
    /panfs/pfs.local/work/bi/bin/bwa/bwa mem pseudochromosomes.fasta ${sample}.fq.gz |
      /panfs/pfs.local/work/bi/bin/samtools-1.3.1/bin/samtools view -b |
      /panfs/pfs.local/work/bi/bin/samtools-1.3.1/bin/samtools sort > ${sample}.bam
done


#Run gstacks to build loci from the aligned paired-end data.
#We have instructed gstacks to remove any PCR duplicates that it finds.
/home/d669d153/work/stacks-2.41/gstacks -I $src/ -M $src/aph.popmap.txt -O $src/ --details

# Run populations and export a vcf. Do filtering steps on the output vcf.
/home/d669d153/work/stacks-2.41/populations -P $src/ -M $src/aph.popmap.txt -O $src/ --vcf

