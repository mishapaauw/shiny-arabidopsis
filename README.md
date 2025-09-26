# Shiny Arabidopsis

Building a shiny app to visualize Arabidopsis chromosomes and annotations. This repo acts as practice to develop a GWAS results visualization dashboard. The idea is to visualize a manhattan plot of significant kmers, which you can use to zoom in to a region of interest in a genome browser.

## igvShiny

The app uses the Shiny implementation of igv. We need some features of the development version (at the time of writing) of igvShiny, which we can install using:

```R
> install_github('gladkia/igvShiny', ref = '4e39a83')
> library(igvShiny)
> packageVersion("igvShiny")
[1] ‘1.5.0’
```

## Data 

- The Arabidopsis genome annotation was downloaded from `https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/arabidopsis_thaliana/`
- The Arabidopsis genome was also downloaded from ENSEMBL, and indexed using `samtools`:

```bash
cd data/
curl -o Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
```

## TODOs

The following things need to be implemented

*Data preparation*:

- Move all the data preparation steps to a seperate R script, and save everything needed by Shiny in one RData file object.
- Count number of significant kmers and add this to the valueboxes (we need the 5per and 10per file from the kmerGWAS output for this)
- Also keep track of how many kmers were mapped/were not mapped.

*BAM file visualization*:

- Add read name to the BAM file, somehow. Probably we should do this in the bowtie mapping step
- Make BAM file load upon opening. Currently there's the need to a button to trigger the load.

*Phenotype visualization*:

- Make phenotype distribution histogram: done
- Make sure that the cmd line scripts copy both files ('pheno.phenotypes' and 'pheno.used_phenotypes') to the results folder 
    
- Make 'candlestick chart' of selected kmers, with option to highlight genome of interest
    - For this, we need a mapping file between nordberg IDs and accession names.
    - how to select kmers? probably not possible in the IGV window.
    - Need to get presence absence pattern for all kmers instead of top100. (and turn of automated plotting as currently implemented)

When everything is implemented, we move to the design phase where we optimize graph appearance and overal dashboard appearance. 
    - for this, use card layout from `bslib` package.

