# --- Marker Analysis Module ---
# VERSION 4:
# - Fixes bug in "Full Results Table" (Aggregates signatures).
# - Adds "Dynamic Category Analysis" (Multi-select, separated table).
# - NEW: Adds a global "Data Source for Analysis" switch (Raw vs. Normalized)
#   for validation and testing. This switch controls all analysis buttons
#   and the main data download button.

library(shiny)
library(DT)
library(ggplot2)
library(plotly)
library(pheatmap)
library(shinycssloaders)
library(reshape2) 

marker_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h3("Analysis Controls"),
      fileInput(ns("markers_file"), "Upload Markers File (.csv)", multiple=F, accept=c(".csv",".txt")),
      fileInput(ns("annotator_file"), "Upload Annotator File (.csv)", multiple=F, accept=c(".csv",".txt")),
      
      # --- NEW VALIDATION FEATURE ---
      radioButtons(ns("data_source"), "Data Source for Analysis:",
                   choices = c("Normalized Data" = "norm", "Raw Data" = "raw"),
                   selected = "norm"),
      # --- END NEW FEATURE ---
      
      tags$hr(),
      h4("Single Comparison"),
      uiOutput(ns("group_var_ui")),
      uiOutput(ns("group_level_ui")),
      actionButton(ns("run_analysis"), "Run Single Analysis", class="btn-primary"),
      tags$hr(),
      h4("Meta-Matrix Controls"),
      uiOutput(ns("multi_group_ui")),
      actionButton(ns("run_multi_analysis"), "Run Meta-Matrix", class="btn-info"),
      tags$hr(),
      
      h4("Dynamic Category Analysis"),
      p("Run analysis on all unique values from one or more annotator columns."),
      uiOutput(ns("annotator_col_ui")), # This will now be a selectizeInput
      actionButton(ns("run_dynamic_analysis"), "Run Dynamic Analysis", class="btn-warning"),
      tags$hr(),
      
      selectInput(ns("padj_method"), "P-Value Correction", choices=c("FDR (BH)"="BH", "Bonferroni"="bonferroni", "None"="none"), selected="BH"),
      numericInput(ns("fdr_threshold"), "Significance (FDR) Threshold", value=0.05, min=0.001, max=1, step=0.01),
      tags$hr(),
      downloadButton(ns("download_data"), "Download Data + Signatures", class="btn-success"),
      width=3
    ),
    mainPanel(
      h3("Analysis Results"),
      tabsetPanel(id=ns("results_tabs"),
                  tabPanel("Signature List", DT::dataTableOutput(ns("signature_summary_table")) %>% withSpinner()),
                  tabPanel("Full Results Table", DT::dataTableOutput(ns("full_results_table")) %>% withSpinner()),
                  tabPanel("Volcano Plot", plotlyOutput(ns("volcano_plot"), height="600px") %>% withSpinner()),
                  tabPanel("Heatmap", plotOutput(ns("heatmap"), height="800px") %>% withSpinner()),
                  tabPanel("Feature Viewer",
                           tags$br(),
                           fluidRow(column(4, radioButtons(ns("feature_type"), "Plot Type:", choices=c("Gene","Signature"), inline=T)), column(4, uiOutput(ns("feature_selector_ui")))),
                           tags$hr(),
                           fluidRow(column(6, h4("Raw Data"), plotOutput(ns("plot_raw"), height="500px") %>% withSpinner()), column(6, h4("Normalized Data"), plotOutput(ns("plot_norm"), height="500px") %>% withSpinner()))
                  ),
                  tabPanel("Meta-Matrix",
                           tags$br(),
                           radioButtons(ns("matrix_val"), "Show in Matrix:", choices=c("Log2 Fold Change"="logfc", "P-Value"="pval", "FDR"="fdr"), inline=T),
                           plotOutput(ns("meta_heatmap"), height="800px") %>% withSpinner()
                  )
      )
    )
  )
}

marker_server <- function(id, loaded_data, normalized_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns 
    annotator_df <- reactiveVal(NULL)
    signature_defs <- reactiveVal(NULL) 
    analysis_results <- reactiveVal(NULL)
    multi_results <- reactiveVal(NULL)
    
    use_data <- reactive({
      req(loaded_data(), normalized_data())
      s_info <- loaded_data()$sample_info; n_data <- normalized_data(); n_samples <- nrow(s_info); p_data <- NULL
      if (is.list(n_data) && !is.data.frame(n_data)) {
        if ("protein_data" %in% names(n_data) && ncol(as.data.frame(n_data$protein_data))==n_samples) p_data <- n_data$protein_data
        else if ("normalized" %in% names(n_data) && ncol(as.data.frame(n_data$normalized))==n_samples) p_data <- n_data$normalized
        else { for (i in seq_along(n_data)) { if (length(dim(n_data[[i]]))==2 && ncol(n_data[[i]])==n_samples) { p_data <- n_data[[i]]; break } } }
      } else if (is.data.frame(n_data) || is.matrix(n_data)) { if (ncol(n_data)==n_samples) p_data <- n_data }
      if (is.null(p_data)) stop("Could not find matching data matrix.")
      p_data <- as.data.frame(p_data)
      s_info$Sample <- trimws(as.character(s_info$Sample))
      colnames(p_data) <- trimws(as.character(colnames(p_data)))
      common <- intersect(s_info$Sample, colnames(p_data))
      if (length(common)==0) return(NULL)
      list(sample_info=s_info[s_info$Sample %in% common,], raw_data=as.data.frame(loaded_data()$protein_data)[,common,drop=F], norm_data=p_data[,common,drop=F])
    })
    
    # --- NEW REACTIVE (V4) ---
    # This reactive returns the correct data matrix (raw or norm)
    # based on the user's selection.
    analysis_data <- reactive({
      req(use_data())
      if (input$data_source == "raw") {
        return(use_data()$raw_data)
      } else {
        return(use_data()$norm_data)
      }
    })
    # --- END NEW REACTIVE ---
    
    observeEvent(input$annotator_file, { req(input$annotator_file); tryCatch({ df <- read.csv(input$annotator_file$datapath, stringsAsFactors=F); if(ncol(df)<2) df <- read.csv(input$annotator_file$datapath, stringsAsFactors=F, sep=";"); cols <- colnames(df); gene_col <- cols[grep("^gene_?symbol$", cols, ignore.case=T)][1]; if(is.na(gene_col)) { showNotification("Annotator needs Gene_Symbol column", type="error"); annotator_df(NULL) } else { colnames(df)[colnames(df)==gene_col] <- "Gene_Symbol"; annotator_df(df); showNotification(paste("Loaded Annotator:", nrow(df)), type="message") } }, error=function(e) showNotification(paste("Annotator Error:", e$message), type="error")) })
    
    observe({
      req(input$markers_file, annotator_df())
      id_notif <- showNotification("Parsing markers file...", duration=NULL, closeButton=F); on.exit(removeNotification(id_notif), add=T)
      tryCatch({
        markers_df <- read.csv(input$markers_file$datapath, stringsAsFactors=F)
        if(ncol(markers_df)<3) markers_df <- read.csv(input$markers_file$datapath, stringsAsFactors=F, sep=";")
        annot <- annotator_df(); defs <- list(); cols <- colnames(markers_df)
        sig_idx <- grep("level1", cols, ignore.case=T)[1]
        if(is.na(sig_idx)) { pos <- grep("signature", cols, ignore.case=T); pos <- pos[!grepl("parent", cols[pos], ignore.case=T)]; if(length(pos)>0) sig_idx <- pos[1] }
        rule_idx <- grep("GENEs|gene", cols, ignore.case=T)[1]
        calc_idx <- setdiff(grep("calculation|calc", cols, ignore.case=T), rule_idx)[1]
        if(is.na(sig_idx)||is.na(rule_idx)||is.na(calc_idx)) return(); 
        
        for(i in 1:nrow(markers_df)) {
          row <- markers_df[i,]; sig <- trimws(row[[sig_idx]]); rule <- row[[rule_idx]]; calc <- tolower(trimws(row[[calc_idx]]))
          if(is.na(sig)||sig=="") next; key <- tolower(sig)
          if(!is.na(calc)&&calc=="sum all") {
            genes <- character(0)
            if(grepl("^Division:", rule, ignore.case=T)) { val <- trimws(sub("^Division:","",rule,ignore.case=T)); genes <- annot$Gene_Symbol[tolower(trimws(annot$Division))==tolower(val)] }
            else if(grepl("^Category:", rule, ignore.case=T)) { val <- trimws(sub("^Category:","",rule,ignore.case=T)); genes <- annot$Gene_Symbol[tolower(trimws(annot$Category))==tolower(val)] }
            else if(grepl("^Subcategory:", rule, ignore.case=T)) { val <- trimws(sub("^Subcategory:","",rule,ignore.case=T)); genes <- annot$Gene_Symbol[tolower(trimws(annot$Subcategory))==tolower(val)] }
            else { genes <- unique(unlist(strsplit(rule, "[,;]\\s*"))); genes <- genes[!is.na(genes)&genes!=""] }
            
            defs[[key]] <- list(type="sum", genes=genes, name=sig, valid=(length(genes)>0), 
                                note=if(length(genes)>0) "OK" else "No matching proteins found")
            
          } else if(!is.na(calc)&&grepl("/",calc)) {
            parts <- trimws(unlist(strsplit(calc,"/")))
            if(length(parts)==2) defs[[key]] <- list(type="ratio", num=parts[1], denom=parts[2], name=sig, valid=TRUE, note="OK")
          }
        }
        
        signature_defs(defs)
        showNotification(paste("Parsed", length(defs), "signatures."), type="message")
      }, error=function(e) showNotification(paste("Marker Error:", e$message), type="error"))
    })
    
    output$group_var_ui <- renderUI({ req(use_data()); cols <- colnames(use_data()$sample_info); selectInput(ns("group_var"), "Single: Grouping Variable", choices=setdiff(cols, c("Sample","TotalProteins")), selected="Group") })
    output$group_level_ui <- renderUI({ req(use_data(), input$group_var); levels <- unique(use_data()$sample_info[[input$group_var]]); levels <- levels[!is.na(levels)&levels!=""]; if(length(levels)<2) return(NULL); tagList(selectInput(ns("group1"),"Control Group",choices=levels,selected=levels[1]), selectInput(ns("group2"),"Case Group",choices=levels,selected=if(length(levels)>1) levels[2] else levels[1])) })
    output$multi_group_ui <- renderUI({ req(use_data()); cols <- colnames(use_data()$sample_info); valid <- cols[sapply(cols, function(c) { l <- length(unique(use_data()$sample_info[[c]])); l>=2 && l<=15 })]; selectizeInput(ns("meta_cols"), "Select Variables:", choices=setdiff(valid, c("Sample","TotalProteins")), multiple=T, options=list(placeholder="Select metadata...")) })
    
    output$annotator_col_ui <- renderUI({
      req(annotator_df())
      cols <- colnames(annotator_df())
      valid_cols <- setdiff(cols, "Gene_Symbol")
      selectizeInput(ns("dynamic_col"), 
                     "Group By Annotator Column(s):", 
                     choices = valid_cols, 
                     multiple = TRUE,
                     options = list(placeholder = 'Select one or more columns...'))
    })
    
    # --- UPDATED (V4) ---
    # Now uses analysis_data() instead of use_data()$norm_data
    observeEvent(input$run_analysis, { 
      req(use_data(), analysis_data(), signature_defs(), input$group_var, input$group1, input$group2)
      if(input$group1==input$group2){showNotification("Groups must be different",type="warning");return()}
      withProgress(message=paste('Running analysis on', input$data_source, 'data...'), value=0.2, { 
        tryCatch({ 
          res <- run_signature_analysis(analysis_data(), use_data()$sample_info, signature_defs(), input$group_var, input$group1, input$group2, input$padj_method, input$fdr_threshold, annotator_df())
          analysis_results(res)
          updateTabsetPanel(session, "results_tabs", selected="Signature List") 
        }, error=function(e) showNotification(paste("Analysis failed:", e$message), type="error")) 
      }) 
    })
    
    # --- UPDATED (V4) ---
    # Now uses analysis_data() instead of use_data()$norm_data
    observeEvent(input$run_multi_analysis, { 
      req(use_data(), analysis_data(), signature_defs(), input$meta_cols)
      withProgress(message=paste('Running Meta-Matrix on', input$data_source, 'data...'), value=0, { 
        all_comps <- list()
        s_info <- use_data()$sample_info
        for(col in input$meta_cols) { 
          levels <- unique(s_info[[col]])
          levels <- levels[!is.na(levels)&levels!=""]
          if(length(levels)<2) next
          ref <- levels[1]
          for(i in 2:length(levels)) { 
            case <- levels[i]
            comp <- paste0(col, ": ", case, " vs ", ref)
            incProgress(1/(length(input$meta_cols)*(length(levels)-1)), detail=comp)
            tryCatch({ 
              res <- run_signature_analysis(analysis_data(), s_info, signature_defs(), col, ref, case, input$padj_method, input$fdr_threshold, NULL)
              if(!is.null(res$full_results)) { 
                df <- res$full_results
                df$Comparison <- comp
                all_comps[[comp]] <- df 
              } 
            }, error=function(e) print(paste("Skipped", comp, ":", e$message))) 
          } 
        }
        if(length(all_comps)>0) { 
          multi_results(do.call(rbind, all_comps))
          updateTabsetPanel(session, "results_tabs", selected="Meta-Matrix") 
        } else { 
          showNotification("No valid comparisons generated.", type="warning") 
        } 
      }) 
    })
    
    # --- UPDATED (V4) ---
    # Now uses analysis_data() instead of use_data()$norm_data
    observeEvent(input$run_dynamic_analysis, {
      req(use_data(), analysis_data(), annotator_df(), input$dynamic_col, input$group_var, input$group1, input$group2)
      
      if(input$group1 == input$group2){
        showNotification("Control and Case groups must be different", type="warning")
        return()
      }
      
      withProgress(message=paste('Running Dynamic Analysis on', input$data_source, 'data...'), value = 0.2, {
        tryCatch({
          
          p_mat <- analysis_data() # <-- UPDATED (V4)
          s_info <- use_data()$sample_info
          annot <- annotator_df()
          
          dynamic_scores <- list()
          all_mappings <- list() 
          
          for (group_col in input$dynamic_col) {
            if (!group_col %in% colnames(annot)) next 
            mapping <- annot[, c("Gene_Symbol", group_col)]
            mapping <- mapping[!is.na(mapping[[group_col]]) & mapping[[group_col]] != "", ]
            unique_categories <- unique(mapping[[group_col]])
            
            for(cat in unique_categories) {
              prots_in_cat <- mapping$Gene_Symbol[mapping[[group_col]] == cat]
              mems <- intersect(prots_in_cat, rownames(p_mat))
              
              if(length(mems) > 0) {
                score_key <- paste(group_col, cat, sep = "::") 
                if (score_key %in% names(dynamic_scores)) {
                  score_key <- paste(score_key, "2", sep="::")
                }
                dynamic_scores[[score_key]] <- colSums(p_mat[mems, , drop=F], na.rm=T)
                all_mappings[[score_key]] <- data.frame(Column = group_col, Signature = cat, stringsAsFactors = F)
              }
            }
          }
          
          if(length(dynamic_scores) == 0) {
            stop("No matching proteins found for any categories in the selected annotator column(s).")
          }
          
          score_mat <- do.call(rbind, dynamic_scores)
          
          g_vals <- s_info[[input$group_var]]
          g1_s <- s_info$Sample[g_vals == input$group1]
          g2_s <- s_info$Sample[g_vals == input$group2]
          if(length(g1_s) < 2 || length(g2_s) < 2) stop("Need >= 2 samples per group.")
          
          sig_stats <- list()
          for(score_key in rownames(score_mat)) {
            scores <- score_mat[score_key, ]
            s_g1 <- scores[g1_s]; s_g2 <- scores[g2_s]
            map_info <- all_mappings[[score_key]]
            
            sig_stats[[score_key]] <- data.frame(
              Annotator_Column = map_info$Column, 
              Signature = map_info$Signature,
              Group1_Mean = mean(s_g1, na.rm=T),
              Group2_Mean = mean(s_g2, na.rm=T),
              Log2_FC = log2(mean(s_g2, na.rm=T) + 1e-12) - log2(mean(s_g1, na.rm=T) + 1e-12),
              P_Value = tryCatch(t.test(s_g1, s_g2)$p.value, error=function(e) NA),
              Notes = "OK (Dynamic)",
              stringsAsFactors = F
            )
          }
          
          sig_full <- do.call(rbind, sig_stats)
          rownames(sig_full) <- NULL
          
          if(!is.null(sig_full) && !all(is.na(sig_full$P_Value))) {
            sig_full$FDR <- p.adjust(sig_full$P_Value, method = input$padj_method)
          } else {
            sig_full$FDR <- NA
          }
          
          analysis_results(list(
            significant = sig_full[which(sig_full$FDR < input$fdr_threshold), ],
            full_results = sig_full,
            all_proteins = data.frame(Protein=c("Dynamic analysis does not populate this table."), stringsAsFactors=F)
          ))
          
          updateTabsetPanel(session, "results_tabs", selected = "Signature List")
          
        }, error = function(e) {
          showNotification(paste("Dynamic analysis failed:", e$message), type="error")
        })
      })
    })
    
    output$meta_heatmap <- renderPlot({ req(multi_results(), input$matrix_val); df <- multi_results(); val_col <- switch(input$matrix_val, "logfc"="Log2_FC", "pval"="P_Value", "fdr"="FDR"); mat <- dcast(df, Signature~Comparison, value.var=val_col); rownames(mat) <- mat$Signature; mat$Signature <- NULL; mat <- as.matrix(mat); mat[is.na(mat)] <- if(input$matrix_val=="logfc") 0 else 1; if(input$matrix_val=="logfc") { limit <- max(abs(mat),na.rm=T); breaks <- seq(-limit,limit,length.out=100); colors <- colorRampPalette(c("navy","white","firebrick"))(100) } else { breaks <- c(seq(0,0.05,length.out=50), seq(0.051,1,length.out=50)); colors <- colorRampPalette(c("firebrick","white","grey90"))(100) }; pheatmap(mat, cluster_rows=T, cluster_cols=F, color=colors, breaks=breaks, display_numbers=T, number_color="black", fontsize_number=8, main=paste("Meta-Matrix:", switch(input$matrix_val, "logfc"="Log2 FC", "pval"="P-Value", "fdr"="FDR"))) })
    output$feature_selector_ui <- renderUI({ req(use_data()); if(input$feature_type=="Gene") { selectizeInput(ns("gene_select"),"Select Gene:",choices=NULL,multiple=F,options=list(placeholder="Type gene...")) } else { req(signature_defs()); selectInput(ns("sig_select"),"Select Signature:",choices=unique(sapply(signature_defs(),`[[`,"name"))) } })
    observe({ req(use_data(), input$feature_type=="Gene"); updateSelectizeInput(session, "gene_select", choices=rownames(use_data()$norm_data), server=T) })
    create_feature_plot <- function(data_mat, s_info, feat, grp, g1, g2, pre) { if(!feat%in%rownames(data_mat)) return(NULL); df <- data.frame(Sample=colnames(data_mat), Value=as.numeric(data_mat[feat,]), Group=s_info[[grp]]); df <- df[df$Group%in%c(g1,g2),]; stats <- aggregate(Value~Group,df,function(x) c(mean=mean(x,na.rm=T),se=sd(x,na.rm=T)/sqrt(length(na.omit(x))))); stats <- do.call(data.frame,stats); names(stats)<-c("Group","mean","se"); ggplot(df,aes(x=Group,y=Value,fill=Group))+geom_bar(data=stats,aes(y=mean),stat="identity",alpha=0.6,width=0.7)+geom_errorbar(data=stats,aes(y=mean,ymin=mean-se,ymax=mean+se),width=0.2,size=0.8)+geom_jitter(width=0.2,size=3,alpha=0.8,color="black")+theme_minimal(base_size=16)+labs(title=paste(pre,feat),y="Intensity/Score",x="")+theme(legend.position="none",plot.title=element_text(hjust=0.5)) }
    
    # These plots are NOT updated by the new switch. They are
    # explicitly to show the original raw vs. norm data.
    output$plot_raw <- renderPlot({ req(use_data(),input$group_var,input$group1,input$group2); dat <- use_data()$raw_data; feat <- if(input$feature_type=="Gene") input$gene_select else { req(input$sig_select); dat <- calculate_scores_only(dat, signature_defs()); input$sig_select }; validate(need(!is.null(feat)&&feat%in%rownames(dat),"Feature not found")); create_feature_plot(dat, use_data()$sample_info, feat, input$group_var, input$group1, input$group2, "Raw:") })
    output$plot_norm <- renderPlot({ req(use_data(),input$group_var,input$group1,input$group2); dat <- use_data()$norm_data; feat <- if(input$feature_type=="Gene") input$gene_select else { req(input$sig_select); dat <- calculate_scores_only(dat, signature_defs()); input$sig_select }; validate(need(!is.null(feat)&&feat%in%rownames(dat),"Feature not found")); create_feature_plot(dat, use_data()$sample_info, feat, input$group_var, input$group1, input$group2, "Normalized:") })
    
    # --- UPDATED (V4) ---
    # Download now respects the analysis_data() reactive and updates filename.
    output$download_data <- downloadHandler(
      filename=function(){
        paste0("Data_with_Signatures_", input$data_source, "_", Sys.Date(), ".csv")
      }, 
      content=function(file){ 
        req(use_data(), analysis_data(), signature_defs())
        withProgress(message='Preparing download...',{ 
          current_data <- analysis_data() # <-- UPDATED (V4)
          scores <- calculate_scores_only(current_data, signature_defs())
          combined_data <- rbind(current_data, scores)
          write.csv(combined_data, file) 
        }) 
      })
    
    output$signature_summary_table <- DT::renderDataTable({
      req(analysis_results()$full_results)
      df <- analysis_results()$full_results
      
      datatable(
        df,
        rownames = FALSE,
        extensions = 'Buttons',
        options = list(
          pageLength = 25,
          ordering  = FALSE,
          scrollX   = TRUE,
          dom       = 'Bfrtip',
          buttons   = list(
            list(
              extend        = 'copy',
              title         = NULL,
              exportOptions = list(modifier = list(page = 'all'))
            ),
            list(
              extend        = 'csv',
              title         = NULL,
              filename      = "Signature_List",
              exportOptions = list(modifier = list(page = 'all'))
            ),
            list(
              extend        = 'excel',
              title         = NULL,
              filename      = "Signature_List",
              exportOptions = list(modifier = list(page = 'all'))
            )
          )
        )
      ) %>%
        formatSignif(
          columns = c("Group1_Mean","Group2_Mean","Log2_FC","P_Value","FDR"),
          digits  = 3
        )
    })
    
    output$full_results_table <- DT::renderDataTable({
      req(analysis_results()$all_proteins)
      
      df <- analysis_results()$all_proteins
      
      cols_to_format <- c("Log2_FC", "P_Value", "FDR")
      
      dt <- datatable(
        df,
        rownames   = FALSE,
        extensions = 'Buttons',
        options    = list(
          scrollX    = TRUE,
          pageLength = 25,
          dom        = 'BfIrtip',
          buttons    = list(
            list(
              extend        = 'copy',
              title         = NULL,
              exportOptions = list(modifier = list(page = 'all'))
            ),
            list(
              extend        = 'csv',
              title         = NULL,
              filename      = "Protein_Level_Results",
              exportOptions = list(modifier = list(page = 'all'))
            ),
            list(
              extend        = 'excel',
              title         = NULL,
              filename      = "Protein_Level_Results",
              exportOptions = list(modifier = list(page = 'all'))
            )
          )
        )
      )
      
      if (all(cols_to_format %in% colnames(df))) {
        dt <- dt %>% formatSignif(columns = cols_to_format, digits = 3)
      }
      
      dt
    })
    
    
    output$volcano_plot <- renderPlotly({ req(analysis_results()$full_results); df <- analysis_results()$full_results; df$Significant <- df$FDR < input$fdr_threshold; plot_ly(data=df, x=~Log2_FC, y=~-log10(P_Value), type='scatter', mode='markers', color=~Significant, colors=c("grey","red"), text=~paste("Sig:", Signature, "<br>FDR:", format(FDR,digits=3)), hoverinfo="text") %>% layout(title="Signature Volcano Plot", xaxis=list(title="Log2 FC"), yaxis=list(title="-log10(P-Value)")) })
    output$heatmap <- renderPlot({ req(analysis_results()$significant, use_data(), input$group_var); sig_res <- analysis_results()$significant; defs <- signature_defs(); valid_sigs_lc <- intersect(tolower(sig_res$Signature), names(defs)); sum_sigs_lc <- valid_sigs_lc[sapply(defs[valid_sigs_lc], function(x) x$type=="sum" && isTRUE(x$valid))]; if(length(sum_sigs_lc)==0) return(NULL); s_data <- use_data(); mat <- data.matrix(s_data$norm_data); sig_prots <- unique(unlist(lapply(sum_sigs_lc, function(key) defs[[key]]$genes))); sig_prots <- intersect(sig_prots, rownames(mat)); if(length(sig_prots)<2) return(NULL); g_vals <- as.character(s_data$sample_info[[input$group_var]]); all_samples <- as.character(s_data$sample_info$Sample); g1s <- all_samples[g_vals==input$group1]; g2s <- all_samples[g_vals==input$group2]; plot_mat <- mat[sig_prots, c(g1s,g2s), drop=F]; plot_mat <- plot_mat[apply(plot_mat,1,var,na.rm=T)>0,,drop=F]; if(nrow(plot_mat)<2) return(NULL); mat_scaled <- t(scale(t(log2(plot_mat+1e-9)))); annot_col <- data.frame(Var=factor(rep(c(input$group1,input$group2), times=c(length(g1s),length(g2s))))); colnames(annot_col) <- input$group_var; rownames(annot_col) <- colnames(mat_scaled); pheatmap(mat_scaled, cluster_cols=F, cluster_rows=T, annotation_col=annot_col, show_rownames=F, show_colnames=F, scale="none", main="Significant Signatures Heatmap") })
  }) 
}

calculate_scores_only <- function(protein_data, sig_defs) {
  p_mat <- data.matrix(protein_data); scores <- list(); scores[["total protein signal"]] <- colSums(p_mat, na.rm=T)
  for(key in names(sig_defs)) { def <- sig_defs[[key]]; if(def$type=="sum" && isTRUE(def$valid)) { mems <- intersect(def$genes, rownames(p_mat)); if(length(mems)>0) scores[[key]] <- colSums(p_mat[mems,,drop=F], na.rm=T) } }
  for(key in names(sig_defs)) { def <- sig_defs[[key]]; if(def$type=="ratio") { num <- tolower(def$num); den <- tolower(def$denom); if(!is.null(scores[[num]]) && !is.null(scores[[den]])) scores[[key]] <- scores[[num]]/scores[[den]] } }
  if(length(scores)>0) { score_mat <- do.call(rbind, scores); rownames(score_mat) <- sapply(rownames(score_mat), function(k) { if(k %in% names(sig_defs)) sig_defs[[k]]$name else k }); ordered_names <- c("total protein signal", sapply(sig_defs, `[[`, "name")); common_names <- intersect(ordered_names, rownames(score_mat)); return(score_mat[common_names, , drop=FALSE]) }
  return(NULL)
}

run_signature_analysis <- function(protein_data, sample_info, sig_defs, group_col, g1, g2, padj, fdr, annot_df) {
  p_mat <- data.matrix(protein_data); g_vals <- sample_info[[group_col]]; g1_s <- sample_info$Sample[g_vals==g1]; g2_s <- sample_info$Sample[g_vals==g2]
  if(length(g1_s)<2 || length(g2_s)<2) stop("Need >= 2 samples per group.")
  p_vals <- unlist(apply(p_mat, 1, function(r) tryCatch(t.test(r[g1_s], r[g2_s])$p.value, error=function(e) NA)))
  l2fc <- log2(rowMeans(p_mat[,g2_s,drop=F],na.rm=T)+1e-9) - log2(rowMeans(p_mat[,g1_s,drop=F],na.rm=T)+1e-9)
  prot_res <- data.frame(Protein=rownames(p_mat), Log2_FC=l2fc, P_Value=p_vals, FDR=p.adjust(p_vals, method=padj), stringsAsFactors=F)
  
  # Note: This merge only runs if annot_df is provided (i.e., in "Run Single Analysis")
  if(!is.null(annot_df)) {
    # Get all relevant annotator columns
    annot_cols <- c("Gene_Symbol", "Category", "Division", "Subcategory") # Add any others
    annot_cols <- intersect(annot_cols, colnames(annot_df)) # Find ones that actually exist
    prot_res <- merge(prot_res, annot_df[,annot_cols], by.x="Protein", by.y="Gene_Symbol", all.x=T)
  }
  
  score_mat <- calculate_scores_only(protein_data, sig_defs)
  sig_stats <- list()
  
  name_to_key <- setNames(names(sig_defs), sapply(sig_defs, `[[`, "name"))
  
  for(key in names(sig_defs)) {
    def <- sig_defs[[key]]
    if(key == "total protein signal") next
    
    if(isTRUE(def$valid)) {
      if(!is.null(score_mat) && def$name %in% rownames(score_mat)) {
        scores <- score_mat[def$name, ]
        s_g1 <- scores[g1_s]; s_g2 <- scores[g2_s]
        sig_stats[[length(sig_stats)+1]] <- data.frame(Signature=def$name, Group1_Mean=mean(s_g1,na.rm=T), Group2_Mean=mean(s_g2,na.rm=T), 
                                                       Log2_FC=log2(mean(s_g2,na.rm=T)+1e-12)-log2(mean(s_g1,na.rm=T)+1e-12), 
                                                       P_Value=tryCatch(t.test(s_g1, s_g2)$p.value, error=function(e) NA), 
                                                       Notes="OK", stringsAsFactors=F)
      }
    } else {
      sig_stats[[length(sig_stats)+1]] <- data.frame(Signature=def$name, Group1_Mean=0, Group2_Mean=0, Log2_FC=0, P_Value=NA, Notes=def$note, stringsAsFactors=F)
    }
  }
  
  sig_full <- do.call(rbind, sig_stats)
  if(!is.null(sig_full) && !all(is.na(sig_full$P_Value))) sig_full$FDR <- p.adjust(sig_full$P_Value, method=padj) else sig_full$FDR <- NA
  
  sig_mappings <- list(); for(key in names(sig_defs)) { if(sig_defs[[key]]$type=="sum" && isTRUE(sig_defs[[key]]$valid)) { mems <- intersect(sig_defs[[key]]$genes, rownames(p_mat)); if(length(mems)>0) sig_mappings[[length(sig_mappings)+1]] <- data.frame(Protein=mems, Signature=sig_defs[[key]]$name, stringsAsFactors=F) } }
  
  sig_map_df <- do.call(rbind, sig_mappings)
  if(!is.null(sig_map_df) && nrow(sig_map_df) > 0) {
    agg_sigs <- aggregate(Signature ~ Protein, data = sig_map_df, FUN = function(x) paste(unique(x), collapse = "; "))
    prot_res <- merge(prot_res, agg_sigs, by = "Protein", all.x = T)
  } else {
    prot_res$Signature <- NA
  }
  
  list(significant=sig_full[which(sig_full$FDR < fdr),], full_results=sig_full, all_proteins=prot_res[!is.na(prot_res$Signature),])
}