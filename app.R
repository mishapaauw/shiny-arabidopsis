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

# Fake a GWAS peak, replace with actual data later.
fake_GWAS <- data.frame(
  chr = "Chr2",
  position = runif(40, 2230000,3440000),     
  pval = runif(40, 0.7, 1)
)


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
        geom_blank(data = chromosome_lengths_l, aes(x = bp, y = 0)) + 
        facet_grid(~chr, scales = "free_x", space = "free_x") +
        geom_point(data = fake_GWAS, aes(x = position, y = pval)) + facet_grid(~chr, scales = "free_x", space = "free_x") + theme_manhattan
      
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