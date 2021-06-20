Range-wide Scrub-Jay RADseq investigation
==================================================================================

This repository contains a comprehensive compilation of code and expalantory visualization associated with the paper: "Speciation proceeds in forward and reverse in real time in the Scrub-Jay species complex"

### Data generation
*   Protocols used for DNA extraction from fresh tissue samples, and RAD library prep can be found here:
    > <https://github.com/DevonDeRaad/aph.rad/tree/master/lab.protocols>

### Sequence data to SNPs
*   Scripts used to demultiplex samples, map reads to the FSJ reference genome, run stacks, and filter the output vcf file are here:
    > <https://github.com/DevonDeRaad/aph.rad/tree/master/sequencedata.to.snps>
*   The detailed step by step SNP filtering process with accompanying visualizations can be viewed here:
    > <https://devonderaad.github.io/aph.rad/sequencedata.to.snps/filter.ref.aligned.radstackshelpr.html>

### Access to SNP data
*   The set of 16,307 filtered SNPs for all 95 samples can be accessed here:
    > <https://github.com/DevonDeRaad/aph.rad/blob/master/unzipped.filtered.vcf.gz>
*   The set of 2,725 unlinked, filtered SNPs for all 95 samples can be accessed here:
    > <https://github.com/DevonDeRaad/aph.rad/blob/master/unlinked.filtered.recode.vcf.gz>


Visualizations of processes described in the paper
------------

### Machine learning species delimitation pipeline
*   The machine learning based species delimitation approach we followed can be viewed step by step with visualizations here:
    > <https://devonderaad.github.io/aph.rad/ml.species.delim/ml.species.delimitation.html>

### D frequency spectrum
*   The calculation and visualization of the D frequency spectrum can be viewed step by step with visualizations here:
    > <https://devonderaad.github.io/aph.rad/dfs/dfs.scrub.jays.html>

### Searching for introgressed haplotype blocks
*   Calling diagnostic SNPs, calculating interspecific heterozygosity, hybrid index, and visualizing genotypes can be viewed step by step with visualizations here:
    > <https://devonderaad.github.io/aph.rad/introgress/introgress.scrub.ca.wood.html>

### Quantifying diversity and divergence
*   The steps of calculating pariwise Fst, and visualizing heterozygostiy and nucleotide diversity can be viewed step by step with visualizations here:
    > <https://devonderaad.github.io/aph.rad/nuc.div/nuc.div.calc.fst.html>
