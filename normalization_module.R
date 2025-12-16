# --- Normalization Module ---
# This file contains all UI and Server logic for
# normalizing the raw protein data.

library(shiny)
library(ggplot2)
library(dplyr)
library(ggrepel) # Required for non-overlapping labels

# Bioconductor packages (preprocessCore, limma, vsn) are loaded conditionally

# --- HELPER FUNCTIONS (Defined outside ModuleServer) ---

# 1. Helper to calculate plot height based on sample count
get_plot_height <- function(data) {
  if (is.null(data)) return(400)
  n_samples <- ncol(data)
  return(max(400, n_samples * 20))
}

# 2. Helper to create the boxplot
create_boxplot <- function(data_matrix, title, log2_transform = FALSE) {
  if (is.null(data_matrix)) return(NULL)
  
  plot_matrix <- as.matrix(data_matrix)
  
  if (log2_transform) {
    # Handle zeros before log2 to avoid -Inf
    plot_matrix[plot_matrix == 0] <- NA 
    plot_matrix <- log2(plot_matrix)
  }
  
  # Reshape for ggplot
  data_long <- data.frame(
    Sample = rep(colnames(plot_matrix), each = nrow(plot_matrix)),
    Intensity = as.vector(plot_matrix)
  )
  
  # Remove NAs/Infs for plotting stability
  data_long <- data_long[is.finite(data_long$Intensity), ]
  
  if (nrow(data_long) == 0) return(NULL)
  
  ggplot(data_long, aes(x = Intensity, y = Sample, fill = Sample)) +
    geom_boxplot(outlier.size = 0.5, lwd = 0.3) + 
    theme_bw() +
    # Dynamic Label depending on if we transformed it or not
    labs(title = title, 
         x = if(log2_transform) "Log2(Intensity)" else "Intensity", 
         y = NULL) +
    theme(legend.position = "none")
}


# --- 1. Module UI Function ---
normalization_ui <- function(id) {
  ns <- NS(id) # Namespace
  
  tagList(
    titlePanel("Data Normalization"),
    
    sidebarLayout(
      # --- Sidebar Panel for Inputs ---
      sidebarPanel(
        h4("Normalization Method"),
        p("Select a method to adjust for systematic technical variation."),
        
        radioButtons(ns("norm_method"), "Select Method:",
                     choices = c(
                       "None (Data is already Log2)" = "already_log2",
                       "None (Use Log2 of Raw Data)" = "none",
                       "None (Keep Raw Linear - No Log2)" = "linear_pass", # <--- NEW OPTION
                       "Total Sum Scaling (TSS)" = "tss",
                       "Quantile Normalization" = "quantile",
                       "Median Normalization" = "median",
                       "Mean Centering" = "mean",
                       "Cyclic Loess (limma)" = "loess",
                       "Variance Stabilizing (VSN)" = "vsn"
                     ),
                     selected = "none"
        ),
        
        tags$hr(),
        
        # --- PCA Settings ---
        h4("PCA Settings"),
        p("Filter proteins used in the PCA plot based on missing data."),
        sliderInput(ns("pca_missing_threshold"), 
                    "Max Missing Data (%) Allowed:", 
                    min = 0, max = 100, value = 50, step = 5),
        
        # --- Checkbox for Ellipses and Labels ---
        checkboxInput(ns("show_ellipse"), "Show 95% Confidence Ellipses", value = FALSE),
        checkboxInput(ns("show_labels"), "Show Sample Names", value = FALSE),
        p(em("Note: Ellipses require at least 4 replicates per group.")),
        
        tags$hr(),
        p(strong("Note:")),
        p("Quantile, Loess, and VSN require Bioconductor packages."),
        
        tags$hr(),
        h4("Download Data"),
        downloadButton(ns("download_norm_data"), "Download Normalized Data")
      ),
      
      # --- Main Panel for Outputs ---
      mainPanel(
        
        # 1. Quality Metrics Table
        h4("Quality Metrics (Raw vs. Normalized)"),
        p("Lower CV and higher Correlation generally indicate better normalization."),
        tableOutput(ns("metrics_table")),
        tags$hr(),
        
        # 2. PCA Plot
        h4("PCA Analysis (Normalized Data)"),
        p("This plot shows how samples cluster after normalization. Replicates should group together."),
        p(em("Note: Remaining missing values in filtered proteins are imputed with the minimum observed value.")),
        plotOutput(ns("pca_plot"), height = "500px"),
        tags$hr(),
        
        # 3. Boxplots (Scrollable)
        h4("Intensity Distributions"),
        p("Compare the data distribution before and after normalization."),
        
        fluidRow(
          column(6, 
                 h5("Before Normalization"),
                 uiOutput(ns("ui_plot_before"))
          ),
          column(6, 
                 h5("After Normalization"),
                 uiOutput(ns("ui_plot_after"))
          )
        )
      )
    )
  )
}


# --- 2. Module Server Function ---
normalization_server <- function(id, loaded_data) {
  
  moduleServer(id, function(input, output, session) {
    
    # --- Reactive: Process Raw Data ---
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
        return(NULL)
      })
      req(raw_matrix)
      # Treat 0 as NA for log transformation safety later
      raw_matrix[raw_matrix == 0] <- NA
      return(raw_matrix)
    })
    
    
    # --- Reactive: Perform Normalization ---
    normalized_data_reactive <- reactive({
      req(loaded_data(), raw_data())
      
      method <- input$norm_method
      raw_matrix <- raw_data()
      correct_dimnames <- dimnames(raw_matrix)
      
      final_data <- NULL 
      log2_final_data <- NULL
      
      tryCatch({
        if (method == "already_log2") {
          # Log2 Input -> Pass Through
          log2_final_data <- raw_matrix
          dimnames(log2_final_data) <- correct_dimnames
          
        } else if (method == "linear_pass") {
          # Linear Input -> Pass Through (No Log2)
          # Note: We store this in 'log2_final_data' variable to keep structure consistent
          # even though the data is technically linear.
          log2_final_data <- raw_matrix
          dimnames(log2_final_data) <- correct_dimnames
          
        } else if (method == "none") {
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
          if (requireNamespace("preprocessCore", quietly = TRUE)) {
            final_data <- preprocessCore::normalize.quantiles(raw_matrix)
            dimnames(final_data) <- correct_dimnames
          }
          
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
          if (requireNamespace("limma", quietly = TRUE)) {
            log2_final_data <- limma::normalizeCyclicLoess(log2(raw_matrix), method = "fast")
            dimnames(log2_final_data) <- correct_dimnames
          }
          
        } else if (method == "vsn") {
          if (requireNamespace("vsn", quietly = TRUE)) {
            vsn_fit <- vsn::vsnMatrix(raw_matrix)
            log2_final_data <- vsn::exprs(vsn_fit)
            dimnames(log2_final_data) <- correct_dimnames
          }
        }
        
        # Fill missing log2 if the method returned linear data
        if (is.null(log2_final_data) && !is.null(final_data)) {
          log2_final_data <- log2(final_data)
        }
        
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        return(NULL)
      })
      
      if (is.null(log2_final_data)) return(NULL)
      
      list(
        sample_info = loaded_data()$sample_info,
        protein_data_raw = raw_data(),
        protein_data_log2 = as.data.frame(log2_final_data, check.names = FALSE)
      )
    })
    
    
    # --- Helper: Metrics Calculation ---
    calculate_metrics <- function(data_matrix, sample_info, is_log2) {
      groups <- setNames(sample_info$Group, sample_info$Sample)
      groups <- groups[colnames(data_matrix)]
      
      # CV must be calculated on Linear Scale
      linear_data <- if(is_log2) 2^data_matrix else data_matrix
      
      cvs <- c()
      for(g in unique(groups)) {
        samps <- names(groups)[groups == g]
        if(length(samps) > 1) {
          sub_dat <- linear_data[, samps, drop=FALSE]
          row_means <- rowMeans(sub_dat, na.rm=TRUE)
          row_sds <- apply(sub_dat, 1, sd, na.rm=TRUE)
          valid <- row_means > 0 & !is.na(row_means) & !is.na(row_sds)
          if(sum(valid) > 0) {
            g_cvs <- row_sds[valid] / row_means[valid]
            cvs <- c(cvs, g_cvs)
          }
        }
      }
      median_cv <- median(cvs, na.rm=TRUE) * 100
      
      # Correlation is usually calculated on Log Scale (or just Input scale)
      # We just ensure infinite values are handled
      cor_data <- data_matrix
      cor_data[!is.finite(as.matrix(cor_data))] <- NA
      
      cor_vals <- c()
      suppressWarnings({
        cor_mat <- cor(cor_data, use="pairwise.complete.obs")
      })
      
      for(g in unique(groups)) {
        samps <- names(groups)[groups == g]
        if(length(samps) > 1) {
          sub_cor <- cor_mat[samps, samps]
          vals <- sub_cor[upper.tri(sub_cor)]
          cor_vals <- c(cor_vals, vals)
        }
      }
      mean_cor <- mean(cor_vals, na.rm=TRUE)
      
      return(c(CV = median_cv, COR = mean_cor))
    }
    
    
    # --- Output: Quality Metrics Table ---
    output$metrics_table <- renderTable({
      req(raw_data(), normalized_data_reactive(), loaded_data())
      
      # INTELLIGENT METRICS LOGIC:
      # If user says "Data is already Log2", we assume Input is Log2 (Un-log for CV)
      # If user says anything else (Linear Pass, None, TSS, etc), we assume Input is Linear.
      input_is_log2 <- (input$norm_method == "already_log2")
      
      # For Output data: 
      # "linear_pass" output is still Linear.
      # "already_log2" output is Log2.
      # All others output Log2.
      output_is_log2 <- (input$norm_method != "linear_pass")
      
      raw_m <- calculate_metrics(raw_data(), loaded_data()$sample_info, is_log2 = input_is_log2)
      
      norm_m <- calculate_metrics(normalized_data_reactive()$protein_data_log2, 
                                  loaded_data()$sample_info, is_log2 = output_is_log2)
      
      df <- data.frame(
        Metric = c("Median Intragroup CV (%)", "Mean Intragroup Correlation"),
        `Raw Data` = c(sprintf("%.2f%%", raw_m["CV"]), sprintf("%.3f", raw_m["COR"])),
        `Normalized Data` = c(sprintf("%.2f%%", norm_m["CV"]), sprintf("%.3f", norm_m["COR"]))
      )
      df
    }, striped = TRUE, bordered = TRUE, width = "100%")
    
    
    # --- Output: PCA Plot ---
    output$pca_plot <- renderPlot({
      req(normalized_data_reactive(), input$pca_missing_threshold)
      
      data_log2 <- normalized_data_reactive()$protein_data_log2
      sample_info <- normalized_data_reactive()$sample_info
      
      # 1. Filter proteins
      missing_pct <- rowMeans(is.na(data_log2)) * 100
      keep_rows <- missing_pct <= input$pca_missing_threshold
      data_filtered <- data_log2[keep_rows, ]
      
      validate(
        need(nrow(data_filtered) > 2, 
             paste0("No proteins meet the criteria of having <= ", 
                    input$pca_missing_threshold, "% missing data."))
      )
      
      # 2. Impute for PCA
      data_for_pca <- as.matrix(data_filtered)
      min_val <- min(data_for_pca, na.rm = TRUE)
      if(!is.finite(min_val)) min_val <- 0
      data_for_pca[is.na(data_for_pca)] <- min_val
      
      row_vars <- apply(data_for_pca, 1, var)
      data_for_pca <- data_for_pca[row_vars > 0, ]
      
      validate(
        need(nrow(data_for_pca) > 2, "Not enough variable proteins to run PCA.")
      )
      
      # 3. Run PCA
      pca_res <- prcomp(t(data_for_pca), center = TRUE, scale. = TRUE)
      pca_df <- as.data.frame(pca_res$x)
      pca_df$Sample <- rownames(pca_df)
      pca_df <- merge(pca_df, sample_info[, c("Sample", "Group")], by="Sample")
      
      # Ensure Group is a factor for plotting
      pca_df$Group <- as.factor(pca_df$Group)
      
      var_explained <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
      
      # 4. Build Plot
      gg <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
        geom_point(size = 4, alpha = 0.8) +
        theme_bw(base_size = 14) +
        labs(
          title = paste0("PCA (Method: ", input$norm_method, ")"),
          subtitle = paste0("Using proteins with <= ", input$pca_missing_threshold, "% missing data"),
          x = paste0("PC1 (", var_explained[1], "%)"),
          y = paste0("PC2 (", var_explained[2], "%)")
        )
      
      # --- Ellipses ---
      if (input$show_ellipse) {
        group_counts <- table(pca_df$Group)
        valid_groups <- names(group_counts[group_counts >= 4])
        ellipse_data <- pca_df %>% dplyr::filter(Group %in% valid_groups)
        
        if (nrow(ellipse_data) > 0) {
          tryCatch({
            gg <- gg + stat_ellipse(data = ellipse_data, 
                                    geom = "polygon", 
                                    alpha = 0.2, 
                                    level = 0.95, 
                                    type = "norm", 
                                    show.legend = FALSE)
          }, error = function(e) {
            showNotification("Could not draw ellipses. Showing points only.", type = "warning")
          })
        }
      }
      
      # --- Labels (ggrepel) ---
      if (input$show_labels) {
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          gg <- gg + ggrepel::geom_text_repel(
            aes(label = Sample),
            size = 3.5,
            max.overlaps = Inf, 
            box.padding = 0.5,
            point.padding = 0.3,
            min.segment.length = 0, 
            show.legend = FALSE
          )
        } else {
          gg <- gg + geom_text(aes(label = Sample), vjust = -0.8, size = 3.5, show.legend = FALSE)
        }
      }
      
      return(gg)
    })
    
    
    # --- Output: Boxplots ---
    
    output$plot_before <- renderPlot({
      # INTELLIGENT PLOT LOGIC:
      # Only apply Log2 to "Before" plot if the method involves a Log2 transformation.
      # If we are doing 'already_log2' or 'linear_pass', we should show the raw input as-is.
      should_transform <- !(input$norm_method %in% c("already_log2", "linear_pass"))
      create_boxplot(raw_data(), "", log2_transform = should_transform)
    })
    
    output$ui_plot_before <- renderUI({
      req(raw_data())
      h <- get_plot_height(raw_data())
      div(style = "height: 500px; overflow-y: scroll; border: 1px solid #ddd; padding: 5px;",
          plotOutput(session$ns("plot_before"), height = paste0(h, "px"))
      )
    })
    
    output$plot_after <- renderPlot({
      req(normalized_data_reactive())
      log2_data_matrix <- normalized_data_reactive()$protein_data_log2
      create_boxplot(log2_data_matrix, "", log2_transform = FALSE)
    })
    
    output$ui_plot_after <- renderUI({
      req(normalized_data_reactive())
      h <- get_plot_height(normalized_data_reactive()$protein_data_log2)
      div(style = "height: 500px; overflow-y: scroll; border: 1px solid #ddd; padding: 5px;",
          plotOutput(session$ns("plot_after"), height = paste0(h, "px"))
      )
    })
    
    
    # --- Download Handler ---
    output$download_norm_data <- downloadHandler(
      filename = function() {
        paste0("normalized_data_", input$norm_method, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(loaded_data(), normalized_data_reactive())
        norm_prot_data <- normalized_data_reactive()$protein_data_log2
        sample_info <- loaded_data()$sample_info
        
        meta_df <- sample_info[, !colnames(sample_info) %in% c("Sample", "TotalProteins"), drop = FALSE]
        rownames(meta_df) <- sample_info$Sample
        meta_rows <- as.data.frame(t(meta_df), check.names = FALSE)
        meta_rows$Identifier <- rownames(meta_rows)
        
        norm_prot_data_df <- as.data.frame(norm_prot_data, check.names = FALSE) 
        norm_prot_data_df$Identifier <- rownames(norm_prot_data_df)
        
        sample_cols <- colnames(norm_prot_data) 
        meta_rows_final <- meta_rows[, c("Identifier", sample_cols)]
        norm_prot_data_final <- norm_prot_data_df[, c("Identifier", sample_cols)]
        
        output_df <- rbind(meta_rows_final, norm_prot_data_final)
        write.csv(output_df, file, row.names = FALSE)
      }
    )
    
    return(normalized_data_reactive)
    
  }) # end moduleServer
}