# Shiny Arabidopsis

Building a shiny app to visualize Arabidopsis chromosomes and annotations. Ultimately, I want to develop a GWAS results visualization dashboard. The idea is to visualize a manhattan plot of significant kmers, which you can use to zoom in to a region of interest in a genome browser.

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

## User instructions

Click on a chromosome panel to zoom into that chromosome. Use 'brush selection', i.e. dragging across the x-axis, to zoom the IGV window into the selected area.

**Known bugs**

- If the user starts brushing in the non-zoomed facetted window, it triggers the zoom event already. How to fix this?

## TODOs

**Data preparation**:

- Move all the data preparation steps to a seperate R script, and save everything needed by Shiny in one RData file object.
- Count number of significant kmers and add this to the valueboxes (we need the 5per and 10per file from the kmerGWAS output for this)
- Also keep track of how many kmers were mapped/were not mapped.

**BAM file visualization**:

- Add read name to the BAM file, somehow. Probably we should do this in the bowtie mapping step?
- Make BAM file load upon opening. Currently there's the need to a button to trigger the load.

**Phenotype visualization**:

- Make sure that the cmd line scripts copy both files ('pheno.phenotypes' and 'pheno.used_phenotypes') to the 
results folder
- Highlight Col-0 or other accession in the phenotype distribution graph
- Make 'candlestick chart' of selected kmers, with option to highlight genome of interest
    - For this, we need a mapping file between nordberg IDs and accession names.
    - how to select kmers? probably not possible in the IGV window.
    - Need to get presence absence pattern for all kmers instead of top100. (and turn of automated plotting as currently implemented)

**Download results**:

- Button to download all graphs as currently shown on screen. 

**Finishing touch**:

- When everything is implemented, we move to the design phase where we optimize graph appearance and overal dashboard appearance (using `{bslib}` package).
    - Increase font sizes of the graphs
    - Choose uniform color scheme for all graphs and dashboard
    - Make hexagon shaped logo
    - Come up with a cool name integrating: Arabidopsis, GWAS, kmers, Shiny


