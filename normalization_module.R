# --- Normalization Module ---
# This file contains all UI and Server logic for
# normalizing the raw protein data.

library(shiny)
library(ggplot2) # Added ggplot2 for the boxplot
# Bioconductor packages are loaded conditionally in the server
# to avoid errors if they aren't installed yet.

# --- 1. Module UI Function ---
normalization_ui <- function(id) {
  ns <- NS(id) # Namespace
  
  fluidPage(
    titlePanel("Data Normalization"),
    
    sidebarLayout(
      # --- Sidebar Panel for Inputs ---
      sidebarPanel(
        h4("Normalization Method"),
        p("Select a method to adjust for systematic technical variation
          between samples before analysis."),
        
        radioButtons(ns("norm_method"), "Select Method:",
                     choices = c(
                       "None (Use Log2 of Raw Data)" = "none",
                       "Total Sum Scaling (TSS)" = "tss",
                       "Quantile Normalization" = "quantile",
                       "Median Normalization" = "median",
                       "Mean Centering" = "mean",
                       "Cyclic Loess (limma)" = "loess",
                       "Variance Stabilizing (VSN)" = "vsn"
                       # Spike-in is omitted for now as it requires
                       # more complex UI to select the spike-in gene
                     ),
                     selected = "none"
        ),
        
        tags$hr(),
        
        p(strong("Note:")),
        p("All methods (except VSN and 'None') will be log2-transformed 
          for the 'After Normalization' plot and downstream analysis."),
        p(strong("Quantile, Loess, and VSN")),
        p("require packages from Bioconductor. If the plots don't appear,
          please check the R console or README for installation instructions."),
        
        # --- *** NEW CODE: DOWNLOAD BUTTON *** ---
        tags$hr(),
        h4("Download Data"),
        p("Download the normalized and log2-transformed data in the
          same format as the input file."),
        downloadButton(ns("download_norm_data"), "Download Normalized Data")
        # --- *** END OF NEW CODE *** ---
        
      ),
      
      # --- Main Panel for Outputs ---
      mainPanel(
        h4("Before Normalization (Log2 Scale)"),
        plotOutput(ns("plot_before"), height = "300px"),
        tags$hr(),
        h4("After Normalization (Log2 Scale)"),
        plotOutput(ns("plot_after"), height = "300px")
      )
    )
  )
}


# --- 2. Module Server Function ---
# This server takes the reactive 'loaded_data' from the main app
normalization_server <- function(id, loaded_data) {
  
  moduleServer(id, function(input, output, session) {
    
    # --- Reactive: Process Raw Data ---
    # (This reactive is unchanged from your file)
    raw_data <- reactive({
      req(loaded_data())
      prot_data <- loaded_data()$protein_data
      raw_matrix <- tryCatch({
        raw_matrix_char <- as.matrix(prot_data)
        rn <- rownames(raw_matrix_char)
        cn <- colnames(raw_matrix_char)
        raw_matrix_num <- apply(raw_matrix_char, 2, as.numeric)
        rownames(raw_matrix_num) <- rn
        colnames(raw_matrix_num) <- cn
        raw_matrix_num
      }, error = function(e) {
        showNotification(paste("Failed to convert data to numeric matrix:", e$message), type = "error")
        return(NULL)
      })
      req(raw_matrix)
      raw_matrix[raw_matrix == 0] <- NA
      return(raw_matrix)
    })
    
    
    # --- Reactive: Perform Normalization ---
    # (This reactive is unchanged from your file)
    normalized_data_reactive <- reactive({
      
      req(loaded_data(), raw_data())
      
      method <- input$norm_method
      raw_matrix <- raw_data()
      correct_dimnames <- dimnames(raw_matrix)
      
      final_data <- NULL 
      log2_final_data <- NULL
      
      tryCatch({
        if (method == "none") {
          log2_final_data <- log2(raw_matrix)
        } else if (method == "tss") {
          raw_matrix_for_sum <- raw_matrix
          raw_matrix_for_sum[is.na(raw_matrix_for_sum)] <- 0
          col_sums <- colSums(raw_matrix_for_sum, na.rm = TRUE)
          mean_sum <- mean(col_sums)
          norm_factors <- col_sums / mean_sum
          final_data <- sweep(raw_matrix, 2, norm_factors, "/")
          dimnames(final_data) <- correct_dimnames
        } else if (method == "quantile") {
          if (!requireNamespace("preprocessCore", quietly = TRUE)) {
            showNotification("Please install 'preprocessCore' from Bioconductor.", type = "warning")
            return(NULL)
          }
          final_data <- preprocessCore::normalize.quantiles(raw_matrix)
          dimnames(final_data) <- correct_dimnames
        } else if (method == "median") {
          col_medians <- apply(raw_matrix, 2, median, na.rm = TRUE)
          mean_median <- mean(col_medians)
          norm_factors <- col_medians / mean_median
          final_data <- sweep(raw_matrix, 2, norm_factors, "/")
          dimnames(final_data) <- correct_dimnames
        } else if (method == "mean") {
          log2_raw <- log2(raw_matrix)
          col_means <- colMeans(log2_raw, na.rm = TRUE)
          log2_final_data <- sweep(log2_raw, 2, col_means, "-")
          dimnames(log2_final_data) <- correct_dimnames
        } else if (method == "loess") {
          if (!requireNamespace("limma", quietly = TRUE)) {
            showNotification("Please install 'limma' from Bioconductor.", type = "warning")
            return(NULL)
          }
          log2_final_data <- limma::normalizeCyclicLoess(log2(raw_matrix), method = "fast")
          dimnames(log2_final_data) <- correct_dimnames
        } else if (method == "vsn") {
          if (!requireNamespace("vsn", quietly = TRUE)) {
            showNotification("Please install 'vsn' from Bioconductor.", type = "warning")
            return(NULL)
          }
          vsn_fit <- vsn::vsnMatrix(raw_matrix)
          log2_final_data <- vsn::exprs(vsn_fit)
          dimnames(log2_final_data) <- correct_dimnames
        }
        
        if (is.null(log2_final_data)) {
          log2_final_data <- log2(final_data)
        }
        
      }, error = function(e) {
        showNotification(paste("Error during normalization:", e$message), type = "error")
        return(NULL)
      })
      
      if (is.null(log2_final_data)) return(NULL)
      
      list(
        sample_info = loaded_data()$sample_info,
        protein_data_raw = raw_data(),
        protein_data_log2 = as.data.frame(log2_final_data, check.names = FALSE)
      )
    })
    
    
    # --- Boxplot Helper Function ---
    # (This function is unchanged from your file)
    create_boxplot <- function(data_matrix, title, log2_transform = FALSE) {
      req(data_matrix)
      plot_matrix <- as.matrix(data_matrix)
      if (log2_transform) {
        plot_matrix <- log2(plot_matrix)
      }
      data_long <- data.frame(
        Sample = rep(colnames(plot_matrix), each = nrow(plot_matrix)),
        Intensity = as.vector(plot_matrix)
      )
      ggplot(data_long, aes(x = Sample, y = Intensity, fill = Sample)) +
        geom_boxplot(na.rm = TRUE) +
        theme_minimal() +
        labs(title = title, x = "Sample", y = "Log2(Intensity)") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.position = "none")
    }
    
    # --- Output Plots ---
    # (These are unchanged from your file)
    output$plot_before <- renderPlot({
      create_boxplot(raw_data(), "Before Normalization (Log2 Scale)", log2_transform = TRUE)
    })
    
    output$plot_after <- renderPlot({
      req(normalized_data_reactive())
      log2_data_matrix <- normalized_data_reactive()$protein_data_log2
      create_boxplot(log2_data_matrix, paste("After", input$norm_method, "Normalization (Log2 Scale)"), log2_transform = FALSE)
    })
    
    
    # --- *** NEW CODE: DOWNLOAD HANDLER *** ---
    output$download_norm_data <- downloadHandler(
      
      filename = function() {
        paste0("normalized_data_", input$norm_method, "_", Sys.Date(), ".csv")
      },
      
      content = function(file) {
        # We need both the original metadata and the new normalized data
        req(loaded_data(), normalized_data_reactive())
        
        # 1. Get the normalized protein data
        norm_prot_data <- normalized_data_reactive()$protein_data_log2
        
        # 2. Get the original sample_info
        sample_info <- loaded_data()$sample_info
        
        # 3. Reconstruct the metadata rows (like "Group", "Age_meta")
        #    - Remove the "Sample" and "TotalProteins" columns
        #    - Set rownames to be the Sample names
        #    - Transpose the dataframe
        meta_df <- sample_info[, !colnames(sample_info) %in% c("Sample", "TotalProteins")]
        rownames(meta_df) <- sample_info$Sample
        meta_rows <- as.data.frame(t(meta_df), check.names = FALSE)
        
        # 4. Format both dataframes to have the "Identifier" column
        meta_rows$Identifier <- rownames(meta_rows)
        
        # Coerce normalized data to data.frame just in case
        norm_prot_data_df <- as.data.frame(norm_prot_data, check.names = FALSE) 
        norm_prot_data_df$Identifier <- rownames(norm_prot_data_df)
        
        # 5. Get sample columns in the correct order
        sample_cols <- colnames(norm_prot_data) 
        
        # 6. Re-order columns to put "Identifier" first
        meta_rows_final <- meta_rows[, c("Identifier", sample_cols)]
        norm_prot_data_final <- norm_prot_data_df[, c("Identifier", sample_cols)]
        
        # 7. Combine them back together
        output_df <- rbind(meta_rows_final, norm_prot_data_final)
        
        # 8. Write to CSV
        #    We use row.names = FALSE because the row names ("Group", "ProteinA")
        #    are now properly in the "Identifier" column.
        write.csv(output_df, file, row.names = FALSE)
      }
    )
    # --- *** END OF NEW CODE *** ---
    
    
    # --- Return the reactive data for other modules ---
    return(normalized_data_reactive)
    
  }) # end moduleServer
}

