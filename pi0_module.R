# --- pi0 Estimation Module ---
# This file contains all UI and Server logic for
# estimating the proportion of true null hypotheses (pi0).

library(shiny)
# cp4p is loaded conditionally

# --- 1. Module UI Function ---
pi0_ui <- function(id) {
  ns <- NS(id) # Namespace
  
  sidebarLayout(
    # --- Sidebar Panel for Inputs ---
    sidebarPanel(
      h4("pi0 Estimation"),
      p("This module estimates pi0, the proportion of proteins 
        that are *not* differentially expressed (true nulls)."),
      
      selectInput(ns("pi0_method"), "Select Estimation Method:",
                  choices = c(
                    "ABH (abh)" = "abh",
                    "Storey-BH (st-bh)" = "st-bh",
                    "Meinshausen (mein)" = "mein",
                    "Pounds-Cheng (pounds)" = "pounds"
                  ),
                  selected = "abh"),
      
      actionButton(ns("run_pi0"), "Calculate pi0"),
      
      tags$hr(),
      p(strong("Note:")),
      p("This analysis uses the p-values from the Volcano Plot tab. 
        You must run the Volcano Plot analysis first.")
    ),
    
    # --- Main Panel for Outputs ---
    mainPanel(
      h3("Analysis Results"),
      plotOutput(ns("pvalue_histogram")),
      tags$hr(),
      uiOutput(ns("pi0_interpretation"))
    )
  )
}

# --- 2. Module Server Function ---
# This server takes the reactive 'volcano_results' from the main app
pi0_server <- function(id, volcano_results) {
  
  moduleServer(id, function(input, output, session) {
    
    # --- Reactive: Run pi0 Estimation ---
    pi0_data <- eventReactive(input$run_pi0, {
      
      # Check that volcano plot has been run
      req(volcano_results())
      
      # Check that cp4p is installed
      if (!requireNamespace("cp4p", quietly = TRUE)) {
        showNotification("Please install the 'cp4p' package to use this feature.", 
                         type = "error")
        return(NULL)
      }
      
      # Get p-values and remove NAs
      p_values <- volcano_results()$pvalue
      p_values <- p_values[!is.na(p_values)]
      
      if (length(p_values) == 0) {
        showNotification("No valid p-values found from the volcano analysis.", 
                         type = "warning")
        return(NULL)
      }
      
      # Estimate pi0
      tryCatch({
        pi0_est_result <- cp4p::estim.pi0(p_values, pi0.method = input$pi0_method)
        
        return(list(
          pi0 = pi0_est_result$pi0,
          p_values = p_values
        ))
        
      }, error = function(e) {
        showNotification(paste("Error during pi0 estimation:", e$message), 
                         type = "error")
        return(NULL)
      })
    })
    
    # --- Output: P-Value Histogram ---
    output$pvalue_histogram <- renderPlot({
      req(pi0_data())
      
      p_values <- pi0_data()$p_values
      pi0 <- pi0_data()$pi0
      
      hist(p_values, 
           breaks = 20, 
           main = "Histogram of P-Values", 
           xlab = "P-Value", 
           col = "gray",
           freq = FALSE)
      
      # Add the pi0 line
      abline(h = pi0, col = "blue", lwd = 3, lty = 2)
      legend("topright", 
             legend = paste("pi0 Estimate:", round(pi0, 3)), 
             col = "blue", lty = 2, lwd = 3, bty = "n")
    })
    
    # --- Output: Interpretation Text ---
    output$pi0_interpretation <- renderUI({
      req(pi0_data())
      
      pi0 <- pi0_data()$pi0
      
      tagList(
        h4(paste0("Estimated pi0: ", round(pi0, 4))),
        p(paste0("This analysis estimates that approximately ", 
                 strong(round(pi0 * 100, 1), "%"), 
                 " of the proteins are *not* differentially expressed 
                 (i.e., they are 'true nulls').")),
        p(paste0("Conversely, this suggests that up to ", 
                 strong(round((1 - pi0) * 100, 1), "%"), 
                 " of the proteins *are* differentially expressed 
                 ('true positives')."))
      )
    })
    
  }) # end moduleServer
}
