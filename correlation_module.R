# --- Correlation Module ---
# This file contains all UI and Server logic for
# the correlation analysis.

library(shiny)
library(ggplot2)
library(ggrepel)

# --- 1. Module UI Function ---
correlation_ui <- function(id) {
  ns <- NS(id) # Namespace
  
  fluidPage(
    titlePanel("Correlation Analysis"),
    
    sidebarLayout(
      # --- Sidebar Panel for Inputs ---
      sidebarPanel(
        h4("Correlation Settings"),
        
        # --- *** NEW: Data Source Selection *** ---
        radioButtons(ns("data_source"), "1. Select Data Source:",
                     choices = c("Normalized Data (Recommended)" = "normalized",
                                 "Raw Data" = "raw"),
                     selected = "normalized"),
        
        # --- *** NEW: Zero Handling Selection *** ---
        radioButtons(ns("zero_handling"), "2. Zero Value Handling (for Raw Data):",
                     choices = c("Treat 0 as NA (Recommended)" = "na",
                                 "Keep 0 values" = "keep"),
                     selected = "na"),
        
        # --- *** NEW: Explanatory Note *** ---
        tags$div(class = "well",
                 h5("Notes on Data Selection:"),
                 tags$ul(
                   tags$li(strong("Normalized Data:"), "Uses the log2-transformed data from Tab 2. Zeros are already 'NA'. This is best for most statistical comparisons."),
                   tags$li(strong("Raw Data:"), "Uses the original, non-normalized values. This may be useful if log-transformation is not desired."),
                   tags$li(strong("Treat 0 as NA:"), "Excludes 0-value pairings from correlation. Recommended, as '0' often means 'missing'."),
                   tags$li(strong("Keep 0 values:"), "Treats 0 as a true measurement. Use only if 0 is a biologically meaningful, measured value.")
                 )
        ),
        
        tags$hr(),
        
        h4("Analysis Parameters"),
        
        # We are placing the selectInput directly,
        # and we will update it from the server.
        selectInput(ns("meta_variable"), "3. Select Metadata:",
                    choices = c("Loading..." = "")),
        
        selectInput(ns("cor_method"), "4. Correlation Method:",
                    choices = c("Pearson" = "pearson",
                                "Spearman" = "spearman",
                                "Kendall" = "kendall"),
                    selected = "pearson"),
        
        numericInput(ns("pval_cutoff_1"), "P-Value Cutoff 1 (Orange/Blue):",
                     value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput(ns("pval_cutoff_2"), "P-Value Cutoff 2 (Red/Dark Blue):",
                     value = 0.001, min = 0, max = 1, step = 0.01),
        
        # --- *** NEW: Added numeric input for N labels *** ---
        numericInput(ns("num_labels"), "Number of Top Hits to Label:",
                     value = 15, min = 0, max = 100, step = 1),
        # --- *** END OF NEW UI ELEMENT *** ---
        
        actionButton(ns("run_correlation"), "Run Correlation Analysis"),
        
        width = 3
      ),
      
      # --- Main Panel for Outputs ---
      mainPanel(
        h4("Correlation U-Plot"),
        plotOutput(ns("correlation_plot")),
        
        # --- *** NEW: U-PLOT DOWNLOAD BUTTONS *** ---
        tags$div(style="margin-top: 10px;",
                 downloadButton(ns("download_u_plot_pdf"), "Download Plot (PDF)"),
                 downloadButton(ns("download_u_plot_svg"), "Download Plot (SVG)")
        ),
        # --- *** END OF NEW UI *** ---
        
        tags$hr(),
        
        h4("Correlation Results Table"),
        DT::dataTableOutput(ns("correlation_table")),
        
        tags$hr(),
        
        h4("Download Results"),
        downloadButton(ns("download_all"), "Download All Results"),
        downloadButton(ns("download_pos"), "Download Significant Positive"),
        downloadButton(ns("download_neg"), "Download Significant Negative")
        
        # --- *** NEW: SINGLE PROTEIN SCATTER PLOT UI *** ---
        ,
        tags$hr(),
        h4("Single Protein Scatter Plot"),
        p("Type the exact name of a protein (e.g., from the table above) 
          to see its individual correlation with the selected metadata."),
        
        textInput(ns("protein_name_input"), "Protein Name:",
                  placeholder = "Type exact protein name (e.g., FGA)"),
        
        # --- *** NEW: Color By Dropdown *** ---
        selectInput(ns("color_by_variable"), "Color Points By:",
                    choices = c("None" = "none")),
        # --- *** END OF NEW UI *** ---
        
        # --- *** MODIFIED: Constrained width for square-ish plot *** ---
        div(style = "max-width: 500px; margin-left: 0;",
            plotOutput(ns("scatter_plot"), height = "450px")
        ),
        # --- *** END OF MODIFICATION *** ---
        
        # --- *** NEW: SCATTER PLOT DOWNLOAD BUTTONS *** ---
        tags$div(style="margin-top: 10px; margin-bottom: 30px;",
                 downloadButton(ns("download_scatter_plot_pdf"), "Download Plot (PDF)"),
                 downloadButton(ns("download_scatter_plot_svg"), "Download Plot (SVG)")
        )
        # --- *** END OF NEW UI *** ---
        
      ) # end mainPanel
      
    ) # end sidebarLayout
  ) # end fluidPage
}


# --- 2. Module Server Function ---
correlation_server <- function(id, loaded_data, normalized_data) {
  
  moduleServer(id, function(input, output, session) {
    
    # --- *** NEW: Reactive values to store plot objects *** ---
    u_plot_obj <- reactiveVal(NULL)
    scatter_plot_obj <- reactiveVal(NULL)
    # --- *** END OF NEW REACTIVES *** ---
    
    
    # --- Helper: Check if a column is robustly numeric ---
    is_column_numeric <- function(col) {
      if (is.null(col)) return(FALSE)
      
      # 1. Convert to character (handles factors)
      col_char <- as.character(col)
      
      # 2. *** NEW: Trim whitespace ***
      # This converts " 1 " to "1" and " " to ""
      col_trimmed <- trimws(col_char)
      
      # 3. *** NEW: Convert empty strings "" to NA ***
      # This handles cells that were just spaces
      col_trimmed[col_trimmed == ""] <- NA
      
      # 4. Omit TRUE NAs (now includes blanks, spaces, and empty strings)
      col_no_na <- na.omit(col_trimmed)
      
      # 5. If all values were NA (a blank column), it's fine.
      if (length(col_no_na) == 0) {
        return(TRUE)
      }
      
      # 6. Check if converting the rest produces any new NAs
      # This will now only fail on actual text like "Male" or "Y/N"
      all(!is.na(suppressWarnings(as.numeric(col_no_na))))
    }
    
    
    # --- Reactive: Update Metadata Dropdown ---
    # This observer watches the *loaded_data* reactive.
    observeEvent(loaded_data(), {
      
      load_data <- loaded_data()
      
      if (is.null(load_data) || is.null(load_data$sample_info)) {
        updateSelectInput(session, "meta_variable",
                          choices = c("Load data on Tab 1 first" = ""),
                          selected = "")
        return()
      }
      
      sample_info <- load_data$sample_info
      
      # --- *** NEW: Populate Color By Dropdown *** ---
      # Get ALL metadata columns (including "Group" and non-numeric)
      all_meta_cols <- setdiff(colnames(sample_info), "Sample")
      updateSelectInput(session, "color_by_variable",
                        choices = c("None" = "none", all_meta_cols),
                        selected = "none")
      # --- *** END OF NEW LOGIC *** ---
      
      # Find all metadata columns (exclude "Sample" and "Group")
      meta_cols <- setdiff(colnames(sample_info), c("Sample", "Group"))
      
      numeric_meta_cols <- c() 
      if (length(meta_cols) > 0) {
        # Only keep columns that are numeric
        for (col_name in meta_cols) {
          if (!is.null(sample_info[[col_name]]) && is_column_numeric(sample_info[[col_name]])) {
            numeric_meta_cols <- c(numeric_meta_cols, col_name)
          }
        }
      }
      
      # Now, update the selectInput
      if (is.null(numeric_meta_cols) || length(numeric_meta_cols) == 0) {
        updateSelectInput(session, "meta_variable",
                          choices = c("No numeric metadata found" = ""),
                          selected = "")
      } else {
        updateSelectInput(session, "meta_variable",
                          choices = numeric_meta_cols,
                          selected = numeric_meta_cols[1])
      }
      
    }, ignoreNULL = FALSE, ignoreInit = FALSE) # Run on load
    
    
    # --- *** NEW: Reactive for Data Preparation *** ---
    # This reactive prepares the correct data matrix
    # based on the user's UI selections.
    data_to_correlate <- reactive({
      
      # 1. Require the base data
      req(loaded_data(), normalized_data())
      
      data_source <- input$data_source
      zero_handling <- input$zero_handling
      
      prot_data <- NULL
      
      # 2. Select the Data Source
      if (data_source == "normalized") {
        # Use the log2, normalized data from Tab 2
        # In this data, 0s are already NA
        prot_data <- normalized_data()$protein_data_log2
        
      } else {
        # Use the raw, non-log2, non-normalized data from Tab 1
        prot_data <- loaded_data()$protein_data
        
        # 3. Handle Zeros (if using raw data)
        # Must convert to numeric matrix first
        prot_data <- data.frame(apply(prot_data, 2, as.numeric), 
                                row.names = rownames(prot_data))
        
        if (zero_handling == "na") {
          # Treat 0 values as NA
          prot_data[prot_data == 0] <- NA
        }
        # Else (if "keep"), do nothing and keep the 0s
      }
      
      req(prot_data)
      return(prot_data)
      
    })
    
    
    # --- Reactive: Run Correlation Analysis ---
    # This is an eventReactive, so it only runs when the
    # 'run_correlation' button is pressed.
    correlation_results <- eventReactive(input$run_correlation, {
      
      # 1. --- Validate Inputs ---
      req(loaded_data(), 
          data_to_correlate(),
          input$meta_variable != "")
      
      # 2. --- Get Data ---
      prot_data <- data_to_correlate()
      sample_info <- loaded_data()$sample_info
      meta_name <- input$meta_variable
      
      # Get the metadata vector
      meta_vector <- sample_info[[meta_name]]
      
      # Ensure it's numeric (should be, but as a safeguard)
      meta_vector <- as.numeric(as.character(meta_vector))
      
      # Get protein names (rownames)
      protein_names <- rownames(prot_data)
      
      # 3. --- Run Correlation Loop ---
      # We will loop through each protein (row) and correlate
      # it against the metadata vector.
      
      results_list <- lapply(1:nrow(prot_data), function(i) {
        
        prot_vector <- as.numeric(prot_data[i, ])
        
        # Use tryCatch to handle errors, e.g.,
        # if a protein has 0 variance (all 0s or all NAs)
        # 'cor.test' will fail.
        
        test_result <- tryCatch({
          
          # 'cor.test' is great because it handles 'NA' values
          # automatically by using 'use = "pairwise.complete.obs"'
          cor.test(prot_vector, 
                   meta_vector, 
                   method = input$cor_method,
                   use = "pairwise.complete.obs")
          
        }, error = function(e) {
          # If test fails, return NULL
          return(NULL)
        })
        
        # If the test failed (e.g., 0 variance), skip this protein
        if (is.null(test_result)) {
          return(NULL)
        }
        
        # Return a data.frame (will be row-bound later)
        data.frame(
          Protein = protein_names[i],
          Correlation = test_result$estimate,
          PValue = test_result$p.value
        )
      })
      
      # 4. --- Assemble Final Data Frame ---
      
      # 'rbind' all the individual data.frames from the loop
      results_df <- do.call(rbind, results_list)
      
      # If no results (all proteins failed?), return empty frame
      if (is.null(results_df) || nrow(results_df) == 0) {
        return(data.frame())
      }
      
      # Calculate -log10(p-value) for the plot
      # Handle p-values of 0 (which become Inf)
      results_df$NegLog10PValue <- -log10(results_df$PValue)
      
      # Cap infinite values for plotting
      if (any(is.infinite(results_df$NegLog10PValue))) {
        max_val <- max(results_df$NegLog10PValue[is.finite(results_df$NegLog10PValue)], 100)
        results_df$NegLog10PValue[is.infinite(results_df$NegLog10PValue)] <- max_val * 1.1
      }
      
      # 5. --- Assign Significance ---
      # This adds a column for color-coding the plot
      results_df$Significance <- "Not Significant"
      
      # Cutoff 1 (e.g., p < 0.05)
      results_df$Significance[
        which(results_df$PValue < input$pval_cutoff_1 & results_df$Correlation > 0)
      ] <- "Positive (p < 0.05)"
      
      results_df$Significance[
        which(results_df$PValue < input$pval_cutoff_1 & results_df$Correlation < 0)
      ] <- "Negative (p < 0.05)"
      
      # Cutoff 2 (e.g., p < 0.001)
      results_df$Significance[
        which(results_df$PValue < input$pval_cutoff_2 & results_df$Correlation > 0)
      ] <- "Positive (p < 0.001)"
      
      results_df$Significance[
        which(results_df$PValue < input$pval_cutoff_2 & results_df$Correlation < 0)
      ] <- "Negative (p < 0.001)"
      
      return(results_df)
      
    })
    
    
    # --- Output: Correlation Plot ---
    output$correlation_plot <- renderPlot({
      
      req(correlation_results())
      results_df <- correlation_results()
      
      # Color mapping
      color_map <- c(
        "Not Significant" = "grey",
        "Positive (p < 0.05)" = "orange",
        "Negative (p < 0.05)" = "cyan",
        "Positive (p < 0.001)" = "red",
        "Negative (p < 0.001)" = "blue"
      )
      
      # --- *** NEW: Prepare data for labels *** ---
      # Find the top N most significant proteins
      # Need to use base R 'order' and 'head'
      top_hits <- results_df
      if(nrow(top_hits) > 0) {
        top_hits <- top_hits[order(top_hits$PValue), ]
        top_hits <- head(top_hits, n = input$num_labels)
      }
      # --- *** END OF NEW LABEL DATA *** ---
      
      # Create the U-plot (Correlation vs. -Log10 P-Value)
      p <- ggplot(results_df, aes(x = Correlation, y = NegLog10PValue, color = Significance)) +
        geom_point(alpha = 0.7) +
        
        # --- *** ADDED LINES BACK IN *** ---
        geom_hline(yintercept = -log10(input$pval_cutoff_1),
                   linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = -log10(input$pval_cutoff_2),
                   linetype = "dotted", color = "grey50") +
        # --- *** END OF ADDITION *** ---
        
        scale_color_manual(values = color_map) +
        labs(
          title = "Protein Correlation U-Plot",
          subtitle = paste("Correlating", input$meta_variable, "using", input$cor_method, "method"),
          x = "Correlation (rho or r)",
          y = "-Log10(P-Value)"
        ) +
        theme_minimal() +
        theme(legend.position = "bottom") +
        
        # --- *** MODIFIED: Use top_hits for labels *** ---
        ggrepel::geom_text_repel(
          data = top_hits, # Use the new data frame
          aes(label = Protein),
          size = 3,
          max.overlaps = 15
        )
      
      # --- *** NEW: Store plot object and print *** ---
      u_plot_obj(p)
      print(p)
      # --- *** END OF MODIFICATION *** ---
    })
    
    
    # --- Output: Results Table ---
    output$correlation_table <- DT::renderDataTable({
      req(correlation_results())
      
      # Format for better readability
      display_df <- correlation_results()
      display_df$Correlation <- round(display_df$Correlation, 3)
      display_df$PValue <- format.pval(display_df$PValue, digits = 3)
      
      # Re-order columns
      display_df <- display_df[, c("Protein", "Correlation", "PValue", "NegLog10PValue", "Significance")]
      
      DT::datatable(display_df,
                    rownames = FALSE,
                    options = list(pageLength = 10))
    })
    
    
    # --- Download Handlers ---\
    
    # 1. Download All
    output$download_all <- downloadHandler(
      filename = function() {
        paste0("correlation_all_results_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(correlation_results(), file, row.names = FALSE)
      }
    )
    
    # 2. Download Positive
    output$download_pos <- downloadHandler(
      filename = function() {
        paste0("correlation_positive_sig_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(correlation_results())
        results_df <- correlation_results()
        
        pos_sig <- results_df[
          which(results_df$Correlation > 0 & results_df$PValue < input$pval_cutoff_1), 
        ]
        
        write.csv(pos_sig, file, row.names = FALSE)
      }
    )
    
    # 3. Download Negative
    output$download_neg <- downloadHandler(
      filename = function() {
        paste0("correlation_negative_sig_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(correlation_results())
        results_df <- correlation_results()
        
        neg_sig <- results_df[
          which(results_df$Correlation < 0 & results_df$PValue < input$pval_cutoff_1), 
        ]
        
        write.csv(neg_sig, file, row.names = FALSE)
      }
    )
    
    # --- *** NEW: U-PLOT DOWNLOAD HANDLERS *** ---
    output$download_u_plot_pdf <- downloadHandler(
      filename = function() {
        paste0("correlation_u_plot_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        req(u_plot_obj())
        # Use ggsave to write the stored plot object
        ggsave(file, plot = u_plot_obj(), device = "pdf", width = 10, height = 8)
      }
    )
    
    output$download_u_plot_svg <- downloadHandler(
      filename = function() {
        paste0("correlation_u_plot_", Sys.Date(), ".svg")
      },
      content = function(file) {
        req(u_plot_obj())
        # Use ggsave to write the stored plot object
        ggsave(file, plot = u_plot_obj(), device = "svg", width = 10, height = 8)
      }
    )
    # --- *** END OF NEW HANDLERS *** ---
    
    
    # --- *** NEW: SCATTER PLOT SERVER LOGIC *** ---
    output$scatter_plot <- renderPlot({
      
      # 1. Get user input for protein name
      prot_name <- input$protein_name_input
      
      # 2. Require all data (wait for user to type)
      req(prot_name != "", 
          data_to_correlate(), 
          loaded_data(), 
          # correlation_results(), # Don't req this, so plot can be made first
          input$meta_variable)
      
      # 3. Get data sources
      prot_data <- data_to_correlate()
      meta_data <- loaded_data()$sample_info
      all_results <- correlation_results() # Can be NULL if not run
      meta_name <- input$meta_variable
      color_by_name <- input$color_by_variable # --- *** NEW *** ---
      
      # 4. Check if protein exists
      if (!prot_name %in% rownames(prot_data)) {
        # Show a blank plot with an error message
        p_err <- ggplot() + 
          annotate("text", x = 0, y = 0, 
                   label = paste("Protein not found:", prot_name, 
                                 "\nCheck spelling, case, and data source."), 
                   color = "red", size = 5) +
          theme_void()
        
        scatter_plot_obj(p_err) # Store the error plot
        print(p_err)
        return(NULL)
      }
      
      # 5. Get vectors for plotting
      prot_vec <- as.numeric(prot_data[prot_name, ])
      meta_vec <- as.numeric(meta_data[[meta_name]])
      
      # 6. Get stats from main analysis for the title
      prot_stats <- NULL
      if (!is.null(all_results)) {
        prot_stats <- all_results[all_results$Protein == prot_name, ]
      }
      
      # 7. Create plot data frame
      plot_df <- data.frame(
        Metadata = meta_vec,
        Protein = prot_vec
      )
      
      # --- *** NEW: Add color grouping data *** ---
      if (color_by_name != "none") {
        
        # Get the raw color vector
        raw_color_vector <- meta_data[[color_by_name]]
        
        # --- *** THIS IS THE ROBUST FIX *** ---
        # We will robustly clean the vector first.
        
        # 1. Convert to character, trim whitespace
        cleaned_vector <- trimws(as.character(raw_color_vector))
        
        # 2. Convert ALL forms of blank/NA to a real NA
        cleaned_vector[cleaned_vector == ""] <- NA
        cleaned_vector[cleaned_vector == "NA"] <- NA
        cleaned_vector[cleaned_vector == "na"] <- NA
        cleaned_vector[cleaned_vector == "N/A"] <- NA
        
        # 3. Check if all remaining (non-NA) values are numeric
        is_numeric_like <- all(!is.na(suppressWarnings(as.numeric(na.omit(cleaned_vector)))))
        
        if (is_numeric_like) {
          # --- *** THIS IS THE NEW BINNING LOGIC *** ---
          numeric_vector <- as.numeric(cleaned_vector)
          
          # Check for "Age"
          if (grepl("Age", color_by_name, ignore.case = TRUE)) {
            min_val <- floor(min(numeric_vector, na.rm = TRUE))
            max_val <- ceiling(max(numeric_vector, na.rm = TRUE))
            if (is.finite(min_val) && is.finite(max_val) && min_val < max_val) {
              breaks <- seq(min_val, max_val + 4, by = 5)
              plot_df$ColorGroup <- cut(numeric_vector, breaks = breaks, right = FALSE, include.lowest = TRUE)
            } else {
              plot_df$ColorGroup <- as.factor(numeric_vector) # Fallback
            }
            
            # Check for "BMI"  
          } else if (grepl("BMI", color_by_name, ignore.case = TRUE)) {
            min_val <- floor(min(numeric_vector, na.rm = TRUE))
            max_val <- ceiling(max(numeric_vector, na.rm = TRUE))
            if (is.finite(min_val) && is.finite(max_val) && min_val < max_val) {
              breaks <- seq(min_val, max_val + 2, by = 3)
              plot_df$ColorGroup <- cut(numeric_vector, breaks = breaks, right = FALSE, include.lowest = TRUE)
            } else {
              plot_df$ColorGroup <- as.factor(numeric_vector) # Fallback
            }
            
          } else {
            # It's numeric, but not Age or BMI (e.g., 0/1)
            # Unify 0, 0.0, etc. by converting to factor
            plot_df$ColorGroup <- as.factor(numeric_vector)
          }
          # --- *** END OF BINNING LOGIC *** ---
          
        } else {
          # If it's a true text field (like "Group", "Sex_meta"),
          # just convert the original vector to factor (to preserve "NA" as a string level if it was one)
          plot_df$ColorGroup <- as.factor(raw_color_vector)
        }
        # --- *** END OF FIX *** ---
        
      }
      # --- *** END OF NEW LOGIC *** ---
      
      # 8. Create plot title and subtitle
      plot_title <- paste("Scatter Plot:", prot_name, "vs.", meta_name)
      plot_subtitle <- "Run correlation analysis to see stats."
      
      if(!is.null(prot_stats) && nrow(prot_stats) == 1) {
        plot_subtitle <- paste(
          "Method:", input$cor_method, "|",
          "Correlation:", round(prot_stats$Correlation, 3), "|",
          "P-Value:", format.pval(prot_stats$PValue, digits = 3)
        )
      }
      
      # 9. Generate plot
      # --- *** MODIFIED: Build plot dynamically for color *** ---
      p_scatter <- ggplot(plot_df, aes(x = Metadata, y = Protein))
      
      if (color_by_name != "none") {
        p_scatter <- p_scatter +
          geom_point(aes(color = ColorGroup), alpha = 0.7, na.rm = TRUE) +
          labs(color = color_by_name) # Add legend title
      } else {
        p_scatter <- p_scatter +
          geom_point(alpha = 0.7, na.rm = TRUE)
      }
      
      p_scatter <- p_scatter +
        geom_smooth(method = "lm", formula = y ~ x, na.rm = TRUE, color = "blue") +
        labs(
          title = plot_title,
          subtitle = plot_subtitle,
          x = meta_name,
          y = paste(prot_name, "(Expression)")
        ) +
        theme_minimal()
      # --- *** END OF MODIFICATION *** ---
      
      # --- *** NEW: Store plot object and print *** ---
      scatter_plot_obj(p_scatter)
      print(p_scatter)
      # --- *** END OF MODIFICATION *** ---
      
    }) # end scatter_plot
    
    
    # --- *** NEW: SCATTER PLOT DOWNLOAD HANDLERS *** ---
    output$download_scatter_plot_pdf <- downloadHandler(
      filename = function() {
        prot_name <- gsub("[^a-zA-Z0-9_]", "_", input$protein_name_input)
        meta_name <- gsub("[^a-zA-Z0-9_]", "_", input$meta_variable)
        paste0("scatter_", prot_name, "_vs_", meta_name, "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        req(scatter_plot_obj())
        # Use ggsave to write the stored plot object
        # Using smaller dimensions for a single plot
        ggsave(file, plot = scatter_plot_obj(), device = "pdf", width = 7, height = 6)
      }
    )
    
    output$download_scatter_plot_svg <- downloadHandler(
      filename = function() {
        prot_name <- gsub("[^a-zA-Z0-9_]", "_", input$protein_name_input)
        meta_name <- gsub("[^a-zA-Z0-9_]", "_", input$meta_variable)
        paste0("scatter_", prot_name, "_vs_", meta_name, "_", Sys.Date(), ".svg")
      },
      content = function(file) {
        req(scatter_plot_obj())
        # Use ggsave to write the stored plot object
        ggsave(file, plot = scatter_plot_obj(), device = "svg", width = 7, height = 6)
      }
    )
    # --- *** END OF NEW HANDLERS *** ---
    
  }) # end moduleServer
}





