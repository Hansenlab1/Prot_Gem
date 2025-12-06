# --- Volcano Plot Module ---
# This file contains all UI and Server logic for
# the differential expression volcano plot.

library(shiny)
library(ggplot2)
library(dplyr)
library(ggrepel) # For non-overlapping labels
library(DT)      # For interactive tables

# --- 1. Module UI Function ---
volcano_ui <- function(id) {
  ns <- NS(id) # Namespace
  
  sidebarLayout(
    # --- Sidebar Panel for Inputs ---
    sidebarPanel(
      h4("Differential Expression"),
      p("Select two groups to compare. This will run a t-test for each protein."),
      
      # Group selection dropdowns
      selectInput(ns("group1"), "Select Group 1 (Control):",
                  choices = c("Loading..." = "")),
      
      selectInput(ns("group2"), "Select Group 2 (Treatment):",
                  choices = c("Loading..." = "")),
      
      tags$hr(),
      h5("Plot Controls"),
      
      # Cutoff inputs
      numericInput(ns("p_cutoff"), "P-Value Cutoff (e.g., 0.05):",
                   value = 0.05, min = 0, max = 1, step = 0.01),
      
      numericInput(ns("fc_cutoff"), "Log2 Fold Change Cutoff (e.g., 1):",
                   value = 1, min = 0, step = 0.5),
      
      tags$hr(),
      h5("Labeling"),
      
      # Labeling controls
      numericInput(ns("top_n"), "Label Top N Proteins (by p-value):",
                   value = 10, min = 0, step = 1),
      
      checkboxInput(ns("show_labels"), "Show Protein Labels", value = TRUE),
      checkboxInput(ns("show_lines"), "Show Cutoff Lines", value = TRUE),
      
      # Action button to run
      actionButton(ns("runAnalysis"), "Run Analysis")
    ),
    
    # --- Main Panel for Outputs ---
    mainPanel(
      tabsetPanel(
        type = "tabs",
        # --- Volcano Plot Tab ---
        tabPanel("Volcano Plot & Data",
                 h3("Volcano Plot"),
                 p("This plot shows protein significance (p-value) vs. fold change. Red points are significant based on your cutoffs."),
                 
                 # Plot Output
                 plotOutput(ns("volcanoPlot"), height = "600px"),
                 
                 tags$hr(),
                 h4("Download Plot Image"),
                 p("Adjust the dimensions below to control the aspect ratio of the downloaded file."),
                 
                 # --- Dimension Controls & Download Buttons ---
                 div(
                   style = "display: flex; align-items: flex-end; gap: 15px; background-color: #f9f9f9; padding: 15px; border-radius: 5px;",
                   
                   # Width Input
                   div(
                     style = "width: 100px;",
                     numericInput(ns("plot_width"), "Width (in):", value = 10, min = 2, max = 20)
                   ),
                   
                   # Height Input
                   div(
                     style = "width: 100px;",
                     numericInput(ns("plot_height"), "Height (in):", value = 8, min = 2, max = 20)
                   ),
                   
                   # Download Buttons
                   div(
                     downloadButton(ns("download_plot_pdf"), "Download PDF", class = "btn-default"),
                     downloadButton(ns("download_plot_svg"), "Download SVG", class = "btn-default")
                   )
                 ),
                 
                 # --- Download Data Buttons ---
                 tags$hr(),
                 h4("Download Results CSV"),
                 div(
                   style = "margin-top: 10px; display: flex; gap: 10px; margin-bottom: 20px;",
                   downloadButton(ns("download_all"), "Download All", class = "btn-default"),
                   downloadButton(ns("download_up"), "Download Up-regulated", class = "btn-primary"),
                   downloadButton(ns("download_down"), "Download Down-regulated", class = "btn-info")
                 ),
                 
                 # --- Results Table ---
                 tags$hr(),
                 h4("Top 25 Significant Proteins"),
                 p("Showing the top 25 most significant proteins, sorted by p-value."),
                 DT::DTOutput(ns("results_table"))
        ),
        # --- P-Value Histogram Tab ---
        tabPanel("P-Value Histogram",
                 h3("P-Value Distribution"),
                 p("This histogram shows the distribution of all p-values from the analysis. A 'healthy' distribution often has a spike near 0 and is flat near 1."),
                 plotOutput(ns("pvalue_histogram"), height = "500px")
        )
      )
    )
  )
}


# --- 2. Module Server Function ---
# This server takes the reactive 'normalized_data' from the main app
volcano_server <- function(id, normalized_data) {
  
  moduleServer(id, function(input, output, session) {
    
    # --- Reactive: Update Group Dropdowns ---
    observe({
      req(normalized_data()) 
      group_choices <- unique(normalized_data()$sample_info$Group)
      
      if (length(group_choices) >= 2) {
        updateSelectInput(session, "group1",
                          choices = group_choices,
                          selected = group_choices[1])
        updateSelectInput(session, "group2",
                          choices = group_choices,
                          selected = group_choices[2])
      } else {
        updateSelectInput(session, "group1", choices = c("Not enough groups" = ""))
        updateSelectInput(session, "group2", choices = c("Not enough groups" = ""))
      }
    })
    
    
    # --- Reactive: Run Analysis ---
    analysis_results <- eventReactive(input$runAnalysis, {
      
      req(normalized_data(), input$group1, input$group2)
      
      if (input$group1 == input$group2) {
        showNotification("Group 1 and Group 2 must be different.", type = "error")
        return(NULL)
      }
      
      prot_data <- normalized_data()$protein_data_log2
      sample_info <- normalized_data()$sample_info
      
      samples_g1 <- sample_info$Sample[sample_info$Group == input$group1]
      samples_g2 <- sample_info$Sample[sample_info$Group == input$group2]
      
      if (length(samples_g1) < 2 || length(samples_g2) < 2) {
        showNotification("Each group must have at least 2 samples for a t-test.", type = "error")
        return(NULL)
      }
      
      g1_data <- prot_data[, samples_g1, drop = FALSE]
      g2_data <- prot_data[, samples_g2, drop = FALSE]
      
      results_list <- lapply(rownames(prot_data), function(protein) {
        
        vec_g1 <- as.numeric(g1_data[protein, ])
        vec_g2 <- as.numeric(g2_data[protein, ])
        
        if (sum(!is.na(vec_g1)) < 2 || sum(!is.na(vec_g2)) < 2) {
          return(NULL)
        }
        
        ttest_res <- tryCatch({
          t.test(vec_g1, vec_g2)
        }, error = function(e) {
          return(list(p.value = NA, estimate = c(NA, NA)))
        })
        
        mean_g1 <- mean(vec_g1, na.rm = TRUE)
        mean_g2 <- mean(vec_g2, na.rm = TRUE)
        log2FC <- mean_g2 - mean_g1
        
        return(data.frame(
          Protein = protein,
          log2FC = log2FC,
          pvalue = ttest_res$p.value
        ))
      })
      
      results_df <- do.call(rbind, results_list)
      results_df <- na.omit(results_df) 
      
      results_df$negLog10P <- -log10(results_df$pvalue)
      
      results_df$significant <- "Not Significant"
      
      up_sig <- (results_df$log2FC > input$fc_cutoff) & (results_df$pvalue < input$p_cutoff)
      down_sig <- (results_df$log2FC < -input$fc_cutoff) & (results_df$pvalue < input$p_cutoff)
      
      results_df$significant[up_sig] <- "Up-regulated"
      results_df$significant[down_sig] <- "Down-regulated"
      
      return(results_df)
    })
    
    
    # --- Reactive: Create Volcano Plot Object ---
    volcano_plot_obj <- reactive({
      req(analysis_results())
      
      df <- analysis_results()
      
      color_map <- c(
        "Up-regulated" = "red",
        "Down-regulated" = "blue",
        "Not Significant" = "grey"
      )
      
      label_data <- df %>%
        arrange(pvalue) %>%
        head(input$top_n)
      
      gg <- ggplot(df, aes(x = log2FC, y = negLog10P, color = significant)) +
        geom_point(alpha = 0.6) +
        scale_color_manual(values = color_map) +
        
        # --- FIXED: Explicitly control Y-axis padding ---
        scale_y_continuous(
          limits = c(-0.1, NA), 
          expand = expansion(mult = c(0, 0.05)) # 0% padding bottom, 5% top
        ) +
        
        theme_bw(base_size = 14) +
        labs(
          title = paste("Volcano Plot:", input$group2, "vs.", input$group1),
          x = "log2(Fold Change)",
          y = "-log10(P-Value)",
          color = "Significance"
        )
      
      if (input$show_lines) {
        gg <- gg +
          geom_hline(yintercept = -log10(input$p_cutoff), linetype = "dashed") +
          geom_vline(xintercept = c(-input$fc_cutoff, input$fc_cutoff), linetype = "dashed")
      }
      
      if (input$show_labels && nrow(label_data) > 0) {
        gg <- gg + ggrepel::geom_text_repel(
          data = label_data,
          aes(label = Protein),
          box.padding = 0.5,
          point.padding = 0.5,
          segment.color = 'grey50',
          max.overlaps = 20
        )
      }
      
      return(gg)
    })
    
    
    # --- Output: Render Volcano Plot ---
    output$volcanoPlot <- renderPlot({
      volcano_plot_obj()
    })
    
    
    # --- Output: Download Handlers (Plot) ---
    output$download_plot_pdf <- downloadHandler(
      filename = function() {
        paste0("Volcano_", input$group2, "_vs_", input$group1, ".pdf")
      },
      content = function(file) {
        req(volcano_plot_obj())
        ggsave(file, plot = volcano_plot_obj(), device = "pdf", 
               width = input$plot_width, height = input$plot_height)
      }
    )
    
    output$download_plot_svg <- downloadHandler(
      filename = function() {
        paste0("Volcano_", input$group2, "_vs_", input$group1, ".svg")
      },
      content = function(file) {
        req(volcano_plot_obj())
        ggsave(file, plot = volcano_plot_obj(), device = "svg", 
               width = input$plot_width, height = input$plot_height)
      }
    )
    
    
    # --- Output: P-Value Histogram ---
    output$pvalue_histogram <- renderPlot({
      req(analysis_results())
      
      df <- analysis_results()
      
      ggplot(df, aes(x = pvalue)) +
        geom_histogram(bins = 30, fill = "lightblue", color = "black", boundary = 0) +
        theme_bw(base_size = 14) +
        labs(
          title = "P-Value Distribution",
          x = "P-Value",
          y = "Count"
        )
    })
    
    # --- Output: Results Table ---
    output$results_table <- DT::renderDT({
      req(analysis_results())
      
      df <- analysis_results() %>%
        arrange(pvalue) %>%
        head(25) %>%
        mutate(across(where(is.numeric), ~ round(., 4)))
      
      datatable(df,
                options = list(pageLength = 25, scrollX = TRUE),
                rownames = FALSE)
    })
    
    # --- Output: Download Handlers (CSV) ---
    output$download_all <- downloadHandler(
      filename = function() {
        paste0("volcano_results_all_", input$group2, "_vs_", input$group1, ".csv")
      },
      content = function(file) {
        req(analysis_results())
        write.csv(analysis_results(), file, row.names = FALSE)
      }
    )
    
    output$download_up <- downloadHandler(
      filename = function() {
        paste0("volcano_results_UP_", input$group2, "_vs_", input$group1, ".csv")
      },
      content = function(file) {
        req(analysis_results())
        up_df <- analysis_results() %>%
          filter(significant == "Up-regulated")
        write.csv(up_df, file, row.names = FALSE)
      }
    )
    
    output$download_down <- downloadHandler(
      filename = function() {
        paste0("volcano_results_DOWN_", input$group2, "_vs_", input$group1, ".csv")
      },
      content = function(file) {
        req(analysis_results())
        down_df <- analysis_results() %>%
          filter(significant == "Down-regulated")
        write.csv(down_df, file, row.names = FALSE)
      }
    )
    
    return(analysis_results)
    
  }) # end moduleServer
}