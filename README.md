# Shiny Arabidopsis

Building a shiny app to visualize Arabidopsis chromosomes and annotations. This repo acts as practice to develop a GWAS results visualization dashboard. The idea is to visualize a manhattan plot of significant kmers, which you can use to zoom in to a region of interest in a genome browser.

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

- Count number of significant kmers and add this to the valueboxes
- Make manhattan plot from actual kmer mappings + p values
- Add .bam file of mapped kmers to the genome browser
- Make phenotype distribution histogram
- Make 'candlestick chart' of selected kmers, with option to highlight genome of interest
    - For this, we need a mapping file between nordberg IDs and accession names.

When everything is implemented, we move to the design phase where we optimize graph appearance and overal dashboard appearance.

