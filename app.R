library(igvShiny)
library(ggplot2)

options <- parseAndValidateGenomeSpec(genomeName="tair10")


df <- data.frame(
  position = 1:2000,              # 20 Mbp
  pval = runif(2000, 0, 1)        # fake p-values
)



ui <- shinyUI(fluidPage(plotOutput("manhattan", brush = "coord_brush", height = "300px"),
                        igvShinyOutput('igvShiny'), width = 10))

server <-
  function(input, output, session) {
    
    output$manhattan <- renderPlot({
      ggplot(df, aes(x = position, y = -log10(pval))) +
        geom_point(alpha = 0.3, size = 0.5) +
        labs(x = "Genomic position (Chr1)", y = "-log10(p-value)") +
        theme_minimal()
    })
    
    
    output$igvShiny <- renderIgvShiny({
      igvShiny(options)
    })
    
    
    observeEvent(input$coord_brush, {
      sel <- input$coord_brush
      if (!is.null(sel)) {
        start <- round(sel$xmin)
        end   <- round(sel$xmax)
        chr   <- "Chr1"
        
        session$sendCustomMessage(
          "igvShiny_zoomToRegion",
          list(
            id = "igv",
            region = list(chr = chr, start = start, end = end)
          )
        )
      }
    
    })
  }


shinyApp(ui = ui, server = server)