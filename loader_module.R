# --- Data Loader Module ---
# This file contains all UI and Server logic for
# loading and parsing the mixed CSV data.

library(shiny)
library(DT) # For interactive tables

# --- 1. Module UI Function ---
loader_ui <- function(id) {
  ns <- NS(id) # Namespace
  
  # Use a sidebar layout
  sidebarLayout(
    
    # --- Sidebar Panel for Inputs ---
    sidebarPanel(
      h3("Data Input"),
      
      # Input: File upload
      fileInput(ns("file_upload"), "Upload CSV File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      # Horizontal line
      tags$hr(),
      
      # Input: Checkbox if file has header
      checkboxInput(ns("header"), "File contains header", TRUE),
      
      # Input: Select separator
      radioButtons(ns("sep"), "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ","),
      
      tags$hr(),
      
      actionButton(ns("load_button"), "Load & Validate Data"),
      
      width = 3
    ),
    
    # --- Main Panel for Outputs ---
    mainPanel(
      # This main panel will show validation results
      # and the sample information table.
      
      h3("Validation & Sample Information"),
      
      # --- Validation Messages ---
      # We'll use a uiOutput to dynamically show
      # success or error messages.
      uiOutput(ns("validation_status")),
      
      tags$hr(),
      
      # --- Sample Info Table ---
      h4("Sample Information Table"),
      p("This table summarizes the metadata (e.g., Group, meta_*) 
        and protein counts for each sample."),
      DT::dataTableOutput(ns("sample_info_table"))
    )
  )
}


# --- 2. Module Server Function ---
loader_server <- function(id) {
  
  moduleServer(id, function(input, output, session) {
    
    # --- Reactive Data Store ---
    # This reactiveVal will hold the final, processed data
    # to be used by other modules in the app.
    # It will be a list: list(sample_info = ..., protein_data = ...)
    processed_data <- reactiveVal(NULL)
    
    
    # --- Event: Load & Validate Button ---
    # This code runs only when the button is clicked
    observeEvent(input$load_button, {
      
      # 1. --- Reset ---
      # Clear previous data and validation messages
      processed_data(NULL)
      output$validation_status <- renderUI(NULL)
      
      # 2. --- Check File Input ---
      req(input$file_upload) # Require a file to be uploaded
      
      # 3. --- Read the Raw File ---
      tryCatch({
        # Force all columns to be read as character to prevent
        # automatic type inference (e.g., "1" becoming 1.0)
        raw_df <- read.csv(input$file_upload$datapath,
                           header = input$header,
                           sep = input$sep,
                           row.names = NULL, # Read first col as data
                           check.names = FALSE, # Allow sample names like '102', '103'
                           colClasses = "character") # <--- FIX APPLIED HERE
        
        
        # 4. --- Validation Step 1: Check format ---
        if (nrow(raw_df) < 2) {
          stop("File must have at least 2 rows (for Group and data).")
        }
        
        # 5. --- Separate Metadata and Protein Data ---
        
        # Metadata rows are "Group" or end in "_meta"
        # The first column name (e.g., "Gene") is ignored
        meta_row_indices <- which(raw_df[, 1] == "Group" | 
                                    grepl("_meta$", raw_df[, 1]))
        
        if (length(meta_row_indices) == 0 || !"Group" %in% raw_df[meta_row_indices, 1]) {
          stop("A 'Group' row is required.")
        }
        
        # Protein data is everything *else*
        protein_row_indices <- setdiff(1:nrow(raw_df), meta_row_indices)
        
        # --- Process Metadata ---
        # We need to transpose it so samples are rows
        metadata_t <- as.data.frame(t(raw_df[meta_row_indices, -1]),
                                    check.names = FALSE)
        colnames(metadata_t) <- raw_df[meta_row_indices, 1]
        
        # Add sample names as a column
        metadata_t$Sample <- rownames(metadata_t)
        rownames(metadata_t) <- NULL
        
        # Re-order to put "Sample" first
        metadata_t <- metadata_t[, c("Sample", "Group", 
                                     grep("_meta$", colnames(metadata_t), value = TRUE))]
        
        
        # --- Process Protein Data ---
        protein_data <- raw_df[protein_row_indices, ]
        
        # Use first column for row names
        # Use make.unique to handle duplicate row names
        prot_names <- protein_data[, 1]
        unique_prot_names <- make.unique(as.character(prot_names))
        
        rownames(protein_data) <- unique_prot_names
        protein_data <- protein_data[, -1] # Remove the name column
        
        
        # 6. --- Validation Step 2: Check consistency ---
        if (ncol(protein_data) != nrow(metadata_t)) {
          stop("Mismatch between metadata sample count and protein data sample count.")
        }
        
        
        # 7. --- Calculate Protein Counts ---
        # Convert to numeric, coercing errors to NA
        # Since we read as "character", this conversion is now explicit and necessary
        protein_data_numeric <- data.frame(apply(protein_data, 2, as.numeric), 
                                           row.names = rownames(protein_data))
        
        # Count non-blank or non-zero entries
        protein_counts <- apply(protein_data_numeric, 2, function(sample_col) {
          sum(!is.na(sample_col) & sample_col != 0)
        })
        
        # Add to the sample info table
        metadata_t$TotalProteins <- protein_counts
        
        
        # 8. --- Show Success Message ---
        output$validation_status <- renderUI({
          tags$div(class = "alert alert-success",
                   "Validation Successful!",
                   tags$br(),
                   paste("Loaded", nrow(protein_data), "proteins and", 
                         nrow(metadata_t), "samples."),
                   tags$br(),
                   paste("Found metadata:", 
                         paste(grep("_meta$", colnames(metadata_t), value = TRUE), 
                               collapse = ", "))
          )
        })
        
        # 9. --- Render Sample Info Table ---
        output$sample_info_table <- DT::renderDataTable({
          DT::datatable(metadata_t, 
                        rownames = FALSE,
                        options = list(scrollX = TRUE))
        })
        
        # 10. --- Store Final Data ---
        # We store the *raw* protein data (still characters here, but that's okay
        # as downstream modules will likely convert it as needed, or you can store
        # protein_data_numeric if preferred).
        processed_data(list(
          sample_info = metadata_t,
          protein_data = protein_data 
        ))
        
      }, 
      # --- Error Handling ---
      error = function(e) {
        # Show a user-friendly error message
        output$validation_status <- renderUI({
          tags$div(class = "alert alert-danger",
                   "Error loading file:",
                   tags$br(),
                   e$message)
        })
        # Clear the table and data
        output$sample_info_table <- DT::renderDataTable(NULL)
        processed_data(NULL)
      })
      
    }) # end observeEvent
    
    
    # --- Return the reactive data ---
    # Return a reactive expression that calls the reactiveVal
    return(reactive({ processed_data() }))
    
  }) # end moduleServer
}