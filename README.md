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

### Access to raw data
*   The set of 16,307 filtered SNPs for all 95 samples can be accessed here:
    > <https://github.com/DevonDeRaad/aph.rad/blob/master/unzipped.filtered.vcf.gz>
*   The set of 2,725 unlinked, filtered SNPs for all 95 samples can be accessed here:
    > <https://github.com/DevonDeRaad/aph.rad/blob/master/unlinked.filtered.recode.vcf.gz>
*   The phylip file containing all sites (including invariant) for 2,725 unlinked loci (used as input for raxml, svdquartets, and to generate gene sequences for building ASTRAL tree) can be accessed here:
    > <https://github.com/DevonDeRaad/aph.rad/blob/master/raxml/populations.all.phylip.gz>


Tree-building approaches
------------

### Raxml
*   Input phylip with all sites, code for running raxml with 1,000 bootstrap replicates, and output tree can be accessed here:
    > <https://github.com/DevonDeRaad/aph.rad/blob/master/unzipped.filtered.vcf.gz>

### SVDquartets
*   Input nexus with all sites, and output tree from SVDquartets can be accessed here:
    > <https://github.com/DevonDeRaad/aph.rad/blob/master/unzipped.filtered.vcf.gz>
*   Tree-building was done using the PAUP* GUI following the directions from this tutorial: <https://github.com/mmatschiner/tutorials/blob/master/species_tree_inference_with_snp_data/README.md>

### ASTRAL-III
*   Phylip with all sites used to generate individual gene sequences in nexus format, code for running IQ-TREE to generate 2,725 input gene trees, code for running ASTRAL with input gene trees and performing polytomy tests, the output trees with posterior probabilities and polytomy test scores, and a .csv file detailing quartet support for each internal branch, can all be accessed here:
    > <https://github.com/DevonDeRaad/aph.rad/tree/master/astral>

### SNAPP
*   Code to generate five unique subsampled nexus files used to create .xml files for SNAPP input via the beauti GUI, code for running SNAPP, output trees and .log files to verify effective sample sizes and MCMC convergence can all be accessed here:
    > <https://github.com/DevonDeRaad/aph.rad/tree/master/snapp>

### Treemix
*   Code to run Treemix and successively add migration edges, visualize and assess the ideal number of migration edges, and perform bootstrap replicates, along with output files, can be accessed here:
    > <https://github.com/DevonDeRaad/aph.rad/tree/master/treemix>


Visualizations of pop-gen approaches to species delimitation and identifying introgression
------------

### Machine learning species delimitation pipeline
*   The machine learning based species delimitation approach we followed can be viewed step by step with visualizations here:
    > <https://devonderaad.github.io/aph.rad/ml.species.delim/ml.species.delimitation.html>

### Dsuite
*   Code to perform ABBA/BABA tests, calculate f4 statistics, and visualize f-branch statistics on guidetree (following tutorial outlined here: <https://github.com/millanek/Dsuite>), as well as all output files can be found here:
    > <https://github.com/DevonDeRaad/aph.rad/tree/master/dsuite>

### D frequency spectrum
*   The calculation and visualization of the D frequency spectrum can be viewed step by step with visualizations here:
    > <https://devonderaad.github.io/aph.rad/dfs/dfs.scrub.jays.html>

### Searching for introgressed haplotype blocks
*   Calling diagnostic SNPs, calculating interspecific heterozygosity, hybrid index, and visualizing genotypes can be viewed step by step with visualizations here:
    > <https://devonderaad.github.io/aph.rad/introgress/introgress.scrub.ca.wood.html>

### Quantifying diversity and divergence
*   The steps of calculating pariwise Fst, and visualizing heterozygostiy and nucleotide diversity can be viewed step by step with visualizations here:
    > <https://devonderaad.github.io/aph.rad/nuc.div/nuc.div.calc.fst.html>


Current iteration of the manuscript
------------
*   Resides here: <https://github.com/DevonDeRaad/aph.rad/blob/master/ScrubJay.manuscript.May.2021.docx>

