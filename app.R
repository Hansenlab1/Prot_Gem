library(shiny)
library(ggplot2)
library(dplyr)
library(ggrepel) # For volcano plot labels

# --- Load Modules ---
source("loader_module.R")
source("normalization_module.R")
source("volcano_module.R")
source("correlation_module.R")
source("pi0_module.R")
source("marker_module.R") # --- ADDED THIS MODULE ---

# --- 1. UI: Main Application UI ---
ui <- fluidPage(
  
  titlePanel("Proteomics Data Analysis Pipeline"),
  
  navbarPage(
    "Analysis Steps:",
    
    # --- Tab 1: Load Data ---
    tabPanel("1. Load Data",
             loader_ui("loader")
    ),
    
    # --- Tab 2: Normalize Data ---
    tabPanel("2. Normalization",
             normalization_ui("norm")
    ),
    
    # --- Tab 3: Volcano Plot ---
    tabPanel("3. Volcano Plot",
             volcano_ui("volcano")
    ),
    
    # --- Tab 4: Correlation Plot ---
    tabPanel("4. Correlation Plot",
             correlation_ui("correlation")
    ),
    
    # --- Tab 5: pi0 Estimation ---
    tabPanel("5. pi0 Estimation",
             pi0_ui("pi0")
    ),
    
    # --- Tab 6: Marker Signatures --- (NEW)
    tabPanel("6. Marker Signatures",
             marker_ui("marker")
    )
    
  ) # end navbarPage
) # end fluidPage


# --- 2. Server ---
server <- function(input, output, session) {
  
  # --- Module 1: Loader ---
  loaded_data <- loader_server("loader")
  
  
  # --- Module 2: Normalization ---
  normalized_data <- normalization_server("norm", loaded_data = loaded_data)
  
  
  # --- Module 3: Volcano Plot ---
  volcano_analysis_results <- volcano_server("volcano", normalized_data)
  
  
  # --- Module 4: Correlation Plot ---
  correlation_server("correlation", 
                     loaded_data = loaded_data, 
                     normalized_data = normalized_data)
  
  
  # --- Module 5: pi0 Estimation ---
  pi0_server("pi0", volcano_analysis_results)
  
  
  # --- Module 6: Marker Signatures --- (NEW)
  marker_server("marker",
                loaded_data = loaded_data,
                normalized_data = normalized_data)
  
} # end server


# --- 3. Run Application ---
shinyApp(ui = ui, server = server)