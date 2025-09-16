library(igvShiny)
library(ggplot2)
library(tidyverse)

source("helpers.R")

# Initialization of genome brwoser
options <- parseAndValidateGenomeSpec(
  genomeName   = "Arabidopsis",
  initialLocus = "Chr2:100000-200000",
  stockGenome = FALSE,
  dataMode = 'localFiles',
  fasta = 'data/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa',
  fastaIndex = 'data/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai',
  genomeAnnotation = 'data/Arabidopsis_thaliana.TAIR10.62.gff3'
)

# Theme options for the manhattan plot
theme_manhattan <- theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),  
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
  ) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Define Arabidopsis chromosome length table to initialize the manhattan plot
chromosome_lengths <- read.csv('data/TAIR10_chromosomes.csv')
chromosome_lengths_l <- chromosome_lengths %>% pivot_longer(cols = c(start, end), values_to="bp", names_to = "position_type")

# Read data and manipulate dataframes to make them compatible with the Manhattan plot canvas
results_folder <- 'casper_gr_tryp_results' # this should be replaced by a newly generated results folder

kmers_locations <- read.table(paste0("results/", results_folder, "/kmer_mapping.tsv"))
colnames(kmers_locations) <- c("kmer","sam_classification", "chromosome", "position")
kmers_pvalues <- read.table(paste0("results/", results_folder, "/kmer_pvalues.tsv"))
colnames(kmers_pvalues) <- c("kmer", "p_val")

# Get the thresholds for significance
threshold_10 <- scan(paste0("results/", results_folder, "/threshold_10per"))
threshold_5 <- scan(paste0("results/", results_folder, "/threshold_5per"))

# Remove kmers without mapping ("*"), and transform the NCBI chromosome headers to Chr1, etc
# This code should be removed if we pick the right reference genome for the mapping
kmers_locations <- kmers_locations %>% filter(chromosome != "*") %>%
  mutate(chr = case_when(
    chromosome == "NC_003070.9" ~ "Chr1",
    chromosome == "NC_003071.7" ~ "Chr2",
    chromosome == "NC_003074.8" ~ "Chr3",
    chromosome == "NC_003075.7" ~ "Chr4",
    chromosome == "NC_003076.8" ~ "Chr5",
    TRUE ~ chromosome 
  ))

# Merge the positions table with p-value table
kmers_merged <- left_join(kmers_locations, kmers_pvalues, by = "kmer")


ui <- fluidPage(
  
  ######### Header section with value boxes
  
  h2("kmer GWAS results"),
  valueBox("42", "Significant kmers (5 per)", color = "#E69F00", width = "150px"),
  valueBox("123", "Significant kmers (10 per)", color = "#56B4E9", width = "150px"),
  
  ######### Manhattan plot section (brushable x axis to determine plotting window of IGV)
  h2("Manhattan plot"),
  plotOutput("manhattan", brush = brushOpts(id = "coord_brush", direction = "x")),
  
  ######## IGV window section
  h2("Genome browser"),
  igvShinyOutput('igvShiny')
  )

server <-
  function(input, output, session) {
    
    # The manhattan plot is made here
    output$manhattan <- renderPlot({
      ggplot() + 
        geom_blank(data = chromosome_lengths_l, aes(x = bp, y = threshold_10 - 0.5)) + 
        facet_grid(~chr, scales = "free_x", space = "free_x") +
        geom_hline(yintercept = threshold_10, col = "#56B3E9", linetype = "dashed") +
        geom_hline(yintercept = threshold_5, col = "#E69F00", linetype = "dashed") +
        
        # I use height jitter to avoid overlapping points with nearly the same location and the same p value
        # Possible edit: make this interactive to enable people to make a 'true' plot and jittered plot.
        geom_jitter(data = kmers_merged, aes(x = position, y = -log10(p_val)), alpha = 0.3, width = 0, height = 0.05) + 
        facet_grid(~chr, scales = "free_x", space = "free_x") + 
        theme_manhattan
      
    })
    
    # The genome browser is initiated here
    output$igvShiny <- renderIgvShiny(
      igvShiny(options)
    )
    
    # Debugging observe function
    observe({
      print(input$coord_brush$xmin)
      print(input$coord_brush$xmax)
      print(input$coord_brush$panelvar1)
    })
    
    
     # Update IGV viewer based on brush selection
  observeEvent(input$coord_brush, {
    
    new_xmin <- input$coord_brush$xmin
    new_xmax <-input$coord_brush$xmax
    new_chromosome <- input$coord_brush$panelvar1
    
    new_locus <- paste0(new_chromosome, ":", new_xmin, "-", new_xmax)
      
    # Update IGV viewer to show the selected region
    showGenomicRegion(session, id = "igvShiny", new_locus)
  })
}


shinyApp(ui = ui, server = server)