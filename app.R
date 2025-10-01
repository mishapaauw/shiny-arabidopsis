library(igvShiny)
library(bslib)
library(ggplot2)
library(tidyverse)
library(rtracklayer)
library(GenomicAlignments)
library(plotly)

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

# Initialize the BAM file
bam_options <- readGAlignments('results/FT10_results/kmer_mapping.bam')
user_config <- list(visibilityWindow = 10000000)

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
results_folder <- 'FT10_results' # this should be replaced by a newly generated results folder

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
    chromosome == "1" ~ "Chr1",
    chromosome == "2" ~ "Chr2",
    chromosome == "3" ~ "Chr3",
    chromosome == "4" ~ "Chr4",
    chromosome == "5" ~ "Chr5",
    TRUE ~ chromosome 
  ))

# Merge the positions table with p-value table
kmers_merged <- left_join(kmers_locations, kmers_pvalues, by = "kmer")

# Load phenotypes (original and used)
original_pheno <- read.table(paste0("results/", results_folder, "/pheno.original_phenotypes"), header = TRUE)
used_pheno <- read.table(paste0("results/", results_folder, "/pheno.phenotypes"), header = TRUE)

original_pheno$dataset = "original"
used_pheno$dataset = "used"

original_pheno$dataset_n = paste0("original (n = ", length(original_pheno$phenotype_value), ")")
used_pheno$dataset_n = paste0("used (n = ", length(used_pheno$phenotype_value), ")")

phenotype <- rbind(original_pheno, used_pheno)
phenotype$dataset <- as.factor(phenotype$dataset)



##### Shiny dashboard starts here

ui <- page_sidebar(
  # Set the theme for the dashboard
  theme = bs_theme(bootswatch = "cosmo"),
  fillable = FALSE,
  
  # Initialize javascript
  useShinyjs(),
  
  # HTML code to be able to define which plots will be able to fade/in out
  tags$head(
    tags$style(HTML("
      .fade-container { transition: opacity 100ms ease-in-out; opacity: 1; }
      .fade-container.fading { opacity: 0 !important; }
    "))
  ),
  
  title = "kmer GWAS dashboard",
  
  # Sidebar contents 
  sidebar = sidebar(
    h4("click to add kmer mapping bam file:"),
    actionButton("addBamLocalFileButton", "BAM local data"),
    h4('Phenotype histogram'),
    checkboxGroupInput("datasets",
                       "Select phenotype dataset(s) to plot:",
                       choices = c("Used" = "used", "Original" = "original"),
                       selected = c("used", "original"))

  ),
  
  # Main dashboard content goes here
  layout_columns(
    
    # Value boxes
    card(card_header("Significant kmers"),
         height = "400px",
         valueBox("42", "Significant kmers (5 per)", color = "#E69F00", width = "150px"),
         valueBox("123", "Significant kmers (10 per)", color = "#56B4E9", width = "150px")),
    
    # Phenotype histogram
    card(card_header("Phenotype histogram"),
         height = "400px",
         plotOutput('histogram')),
    
    # Width ratio for the valueboxes and histogram
    col_widths = c(3,9)),

  # Manhattan plot section
  card(card_header("Manhattan plot"),
       id = "fade-container", class = "fade-container",
       height = "400px",
       plotOutput("manhattan", 
                  brush = brushOpts(id = "coord_brush", direction = "x", delayType = "debounce", delay = 700),
                  click = clickOpts(id = "click_event"),
                  dblclick = dblclickOpts(id = "double_click_event"))
       ),

  # IGV browser
  card(card_header("IGV browser"),
       height = "400px",
       igvShinyOutput('igvShiny')
  )
)

server <-
  function(input, output, session) {
    
    # This code enables the selection of facet via clicking
    selectedFacet <- reactiveVal(NULL)
    fade_ms <- 100
    
    # Definition of fading functions via javascript
    fade_out <- function() {
      runjs("$('#fade-container').addClass('fading');")
    }
    fade_in <- function() {
      runjs("$('#fade-container').removeClass('fading');")
    }
    
    # Upon click, fade out current plot and fade in zoomed in plot
    observeEvent(input$click_event, {
      
      if (!is.null(selectedFacet())) {
        return()     # Ignore single click events while already in zoomed state
      }
    
      facet_clicked <- input$click_event$panelvar1
      
      if (is.null(facet_clicked)) return()
      
      fade_out()
      later(function() {
        selectedFacet(facet_clicked)
        session$onFlush(function() fade_in(), once = TRUE)
      }, fade_ms/1000)
    })
    
    # Double click = zoom back out
    observeEvent(input$double_click_event, {
      fade_out()
      later(function() {
        selectedFacet(NULL)
        session$onFlush(function() fade_in(), once = TRUE)
      }, fade_ms/1000)
    })
    
    # Here we subset the datasets upon clicking
    chromosomes <- reactive({
      if (is.null(selectedFacet())) chromosome_lengths_l else filter(chromosome_lengths_l, chr == selectedFacet())
    })
    sig_kmers <- reactive({
      if (is.null(selectedFacet())) kmers_merged else filter(kmers_merged, chr == selectedFacet())
    })
    
      
    # Phenotype histograms
    output$histogram <- renderPlot({
      selected <- input$datasets
        
        ggplot(phenotype %>% filter(dataset %in% selected), aes(x = phenotype_value, fill = dataset_n)) +
          geom_histogram(position = "identity", alpha = 0.5, bins = 100, colour = 'black') +
          labs(x = "Phenotype value", y = "Count") +
          theme_bw() +
          theme(
            panel.grid.major.x = element_blank(),  
            panel.grid.minor.x = element_blank(),  
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank(),
          )
      
    })
    
    
    
    
    # The manhattan plot is made here
    output$manhattan <- renderPlot({
      ggplot() + 
        geom_blank(data = chromosomes(), aes(x = bp, y = threshold_10 - 0.5)) + 
        facet_grid(~chr, scales = "free_x", space = "free_x") +
        geom_hline(yintercept = threshold_10, col = "#56B3E9", linetype = "dashed") +
        geom_hline(yintercept = threshold_5, col = "#E69F00", linetype = "dashed") +
        
        # I use height jitter to avoid overlapping points with nearly the same location and the same p value
        # Possible edit: make this interactive to enable people to make a 'true' plot and jittered plot.
        geom_jitter(data = sig_kmers(), aes(x = position, y = -log10(p_val)), alpha = 0.3, width = 0, height = 0.05) + 
        facet_grid(~chr, scales = "free_x", space = "free_x") + 
        theme_manhattan
     
    })
    
    
    observeEvent(input$addBamLocalFileButton, {
    loadBamTrackFromLocalData(session,
                              id = "igvShiny",
                              trackName = 'kmer_mapping.bam',
                              data = bam_options,
                              displayMode = "EXPANDED",
                              trackConfig = user_config)
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
      print(input$selected)
    })  
    
    # Update IGV viewer based on brush selection (dragging x coordinates)
    observeEvent(input$coord_brush, {
    
      new_xmin <- input$coord_brush$xmin
      new_xmax <-input$coord_brush$xmax
      new_chromosome <- input$coord_brush$panelvar1
      
      new_locus <- paste0(new_chromosome, ":", new_xmin, "-", new_xmax)
      showGenomicRegion(session, id = "igvShiny", new_locus)
    })
  
}


shinyApp(ui = ui, server = server)
