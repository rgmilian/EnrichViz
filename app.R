# app.R
# Heatmap, Bar Plot, Chord Diagram, and Boxplot visualization of enriched Pathways or Functions
#v.1.0.2
# Copyright 2026 RGM
# MIT License — see README.md for full license text
# ─────────────────────────────────────────────────────────

library(shiny)
library(pheatmap)
library(tidyverse)
library(circlize)

# ── Increase max upload size to 500 MB ──────────────────────────────────────
options(shiny.maxRequestSize = 500 * 1024^2)

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════
ui <- fluidPage(
  
  tags$style(HTML("

    /* ── Stop the whole page from scrolling ── */
    html, body {
      height: 100%;
      overflow: hidden;
    }

    /* ── Sidebar: scrolls independently, stays fixed in place ── */
    .col-sm-3 {
      position : fixed;
      top       : 60px;          /* below the title bar */
      left      : 0;
      height    : calc(100vh - 60px);
      overflow-y: auto;
      overflow-x: hidden;
      padding   : 10px 15px;
      z-index   : 100;
    }

    /* ── Main panel: also fixed, sits to the right of sidebar ── */
    .col-sm-9 {
      position  : fixed;
      top       : 60px;
      left      : 25%;           /* matches sidebar width=3 (3/12 = 25%) */
      width     : 75%;
      height    : calc(100vh - 60px);
      overflow-y: auto;
      overflow-x: hidden;
      padding   : 10px 15px;
    }

    /* ── Title panel stays at the top ── */
    .container-fluid > h2 {
      position  : fixed;
      top       : 0;
      left      : 0;
      width     : 100%;
      background: white;
      z-index   : 200;
      padding   : 10px 20px 5px 20px;
      margin    : 0;
      border-bottom: 1px solid #ddd;
    }

    /* ── Scrollbar styling (optional, cleaner look) ── */
    .col-sm-3::-webkit-scrollbar,
    .col-sm-9::-webkit-scrollbar { width: 6px; }

    .col-sm-3::-webkit-scrollbar-thumb,
    .col-sm-9::-webkit-scrollbar-thumb {
      background   : #cccccc;
      border-radius: 3px;
    }

  ")),
  
  titlePanel("Heatmap, Bar Plot and Chord Diagram of Normalized Data by Pathway or Function"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      
      # ── File uploads ──────────────────────────────────────────────────────
      h4("Upload Data Files"),
      fileInput("file_norm_prot",  "Normalized Data (.csv)",               accept = ".csv"),
      fileInput("file_annotation", "Metadata or Sample Annotation (.csv)", accept = ".csv"),
      fileInput("file_ipa_funct",  "Enriched Pathway or Functions (.csv)", accept = ".csv"),
      
      hr(),
      
      # ── Column-role mapping: IPA Functions file ────────────────────────────
      h4("Select Pathway / Functions — Column"),
      helpText("Category = col 1. Molecules = last column by default."),
      uiOutput("ui_col_category"),
      uiOutput("ui_col_molecules"),
      
      hr(),
      
      # ── Column-role mapping: Normalized counts / abundance file ────────────
      h4("Normalized counts or protein abundance — Column"),
      helpText("Gene/Protein ID column (e.g. 'Gene.Symbol')."),
      uiOutput("ui_col_gene"),
      
      hr(),
      
      # ── Column-role mapping: Annotation file ──────────────────────────────
      h4("Sample Annotation — Column"),
      uiOutput("ui_col_sample_id"),
      uiOutput("ui_col_group"),
      
      hr(),
      
      # ── Bar plot settings ──────────────────────────────────────────────────
      h4("Bar Plot Settings"),
      uiOutput("ui_col_pvalue"),
      uiOutput("ui_pvalue_transform_info"),
      numericInput("bar_top_n",
                   "Show top N pathways",
                   value = 20, min = 5, max = 200),
      textInput("bar_color",
                "Bar fill colour (hex or name)",
                value = "steelblue"),
      numericInput("bar_plot_height",
                   "Bar plot height (px)",
                   value = 600, min = 300, max = 2000),
      downloadButton("download_barplot", "Download Bar Plot (.png)"),
      
      hr(),
      
      # ── Chord diagram settings ─────────────────────────────────────────────
      h4("Chord Diagram Settings"),
      
      numericInput("chord_top_n",
                   "Top N pathways to include",
                   value = 10, min = 2, max = 50),
      
      numericInput("chord_min_freq",
                   "Minimum times a molecule must appear (across pathways)",
                   value = 1, min = 1, max = 20),
      
      helpText("Increase 'minimum appearances' to reduce clutter by keeping",
               "only molecules shared across multiple pathways."),
      
      # ── Inner circle radius ────────────────────────────────────────────────
      sliderInput("chord_inner_radius",
                  "Inner circle size (smaller = more space for labels)",
                  min   = 0.1,
                  max   = 0.7,
                  value = 0.4,
                  step  = 0.05),
      
      helpText("Reduce this value to shrink the central hole and give more",
               "radial space to sector labels."),
      
      # ── Label font sizes ───────────────────────────────────────────────────
      numericInput("chord_label_size_path",
                   "Pathway label font size",
                   value = 0.55, min = 0.3, max = 1.5, step = 0.05),
      
      numericInput("chord_label_size_mol",
                   "Molecule label font size",
                   value = 0.70, min = 0.3, max = 1.5, step = 0.05),
      
      textInput("chord_pathway_color",
                "Pathway sector colour (hex or name)",
                value = "tomato"),
      
      textInput("chord_molecule_color",
                "Molecule sector colour (hex or name)",
                value = "steelblue"),
      
      numericInput("chord_plot_height",
                   "Chord plot height (px)",
                   value = 750, min = 400, max = 2000),
      
      downloadButton("download_chord", "Download Chord Diagram (.png)"),
      
      hr(),
      
      # ── Heatmap settings ───────────────────────────────────────────────────
      h4("Heatmap Settings"),
      
      # ── Category selector moved here ──────────────────────────────────────
      uiOutput("category_selector"),
      
      hr(),
      
      uiOutput("ui_group_colors"),
      numericInput("plot_height",
                   "Heatmap height (px)",
                   value = 700, min = 300, max = 2000),
      
      # ── Heatmap download buttons ───────────────────────────────────────────
      h4("Download Heatmaps"),
      downloadButton("download_single", "Download Selected Heatmap (.png)"),
      br(), br(),
      downloadButton("download_all", "Download ALL Heatmaps (.zip)"),
      
      hr(),
      
      # ── Boxplot settings ───────────────────────────────────────────────────
      h4("Boxplot Settings"),
      uiOutput("ui_boxplot_protein"),
      helpText("Group colours are shared with Heatmap Settings above."),
      numericInput("boxplot_height",
                   "Boxplot height (px)",
                   value = 500, min = 300, max = 2000),
      downloadButton("download_boxplot_sidebar", "Download Boxplot (.png)")
    ),
    
    mainPanel(
      width = 9,
      
      tabsetPanel(
        id = "main_tabs",
        
        # ── Tab 1: Bar Plot ──────────────────────────────────────────────────
        tabPanel(
          "Bar Plot",
          br(),
          verbatimTextOutput("barplot_status"),
          plotOutput("barplot", height = "600px"),
          br(),
          downloadButton("download_barplot_main", "\u2b07 Download Bar Plot (.png)")
        ),
        
        # ── Tab 2: Chord Diagram ─────────────────────────────────────────────
        tabPanel(
          "Chord Diagram",
          br(),
          verbatimTextOutput("chord_status"),
          plotOutput("chord_plot", height = "750px"),
          br(),
          downloadButton("download_chord_main", "\u2b07 Download Chord Diagram (.png)")
        ),
        
        # ── Tab 3: Heatmap ───────────────────────────────────────────────────
        tabPanel(
          "Heatmap",
          br(),
          verbatimTextOutput("status"),
          plotOutput("heatmap_plot", height = "700px"),
          br(),
          downloadButton("download_single_main", "\u2b07 Download Selected Heatmap (.png)"),
          br(), br(),
          downloadButton("download_all_main", "\u2b07 Download ALL Heatmaps (.zip)")
        ),
        
        # ── Tab 4: Boxplot ───────────────────────────────────────────────────
        tabPanel(
          "Boxplot",
          br(),
          verbatimTextOutput("boxplot_status"),
          plotOutput("boxplot_plot", height = "500px"),
          br(),
          downloadButton("download_boxplot_main", "\u2b07 Download Boxplot (.png)")
        )
      )
    )
  )
)

# ══════════════════════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {
  
  # ── 1. Load uploaded data ───────────────────────────────────────────────
  norm_prot  <- reactive({ req(input$file_norm_prot);  read.csv(input$file_norm_prot$datapath,  check.names = FALSE) })
  ipa_funct  <- reactive({ req(input$file_ipa_funct);  read.csv(input$file_ipa_funct$datapath,  check.names = FALSE) })
  annotation <- reactive({ req(input$file_annotation); read.csv(input$file_annotation$datapath, check.names = FALSE) })
  
  # ── 2. Column-selector widgets ──────────────────────────────────────────
  
  output$ui_col_category <- renderUI({
    req(ipa_funct())
    cols <- colnames(ipa_funct())
    selectInput("col_category",
                label   = paste0("Category column (detected: '", cols[1], "')"),
                choices = cols, selected = cols[1])
  })
  
  output$ui_col_molecules <- renderUI({
    req(ipa_funct())
    cols    <- colnames(ipa_funct())
    default <- cols[length(cols)]
    selectInput("col_molecules",
                label   = paste0("Molecules column (detected: '", default, "')"),
                choices = cols, selected = default)
  })
  
  output$ui_col_gene <- renderUI({
    req(norm_prot())
    cols    <- colnames(norm_prot())
    default <- if (length(cols) >= 3) cols[3] else cols[1]
    selectInput("col_gene",
                label   = paste0("Gene/Protein ID column (detected: '", default, "')"),
                choices = cols, selected = default)
  })
  
  output$ui_col_sample_id <- renderUI({
    req(annotation())
    cols <- colnames(annotation())
    selectInput("col_sample_id",
                label   = paste0("Sample ID column (detected: '", cols[1], "')"),
                choices = cols, selected = cols[1])
  })
  
  output$ui_col_group <- renderUI({
    req(annotation())
    cols    <- colnames(annotation())
    default <- if (length(cols) >= 2) cols[2] else cols[1]
    selectInput("col_group",
                label   = paste0("Group column (detected: '", default, "')"),
                choices = cols, selected = default)
  })
  
  # ── 3. P-value column selector ──────────────────────────────────────────
  output$ui_col_pvalue <- renderUI({
    req(ipa_funct())
    cols <- colnames(ipa_funct())
    
    keywords <- c("fdr", "adj", "p\\.val", "pval", "p-val", "p_val",
                  "p\\.adjust", "padj", "log\\(p", "-log", "bh", "bonf")
    
    detected <- cols[sapply(cols, function(x)
      any(sapply(keywords, function(k)
        grepl(k, x, ignore.case = TRUE, perl = TRUE)))
    )]
    
    default <- if (length(detected) > 0) detected[1] else cols[2]
    
    selectInput("col_pvalue",
                label   = "P-value or -log10(p-value) column",
                choices = cols,
                selected = default)
  })
  
  # ── 4. Transformation info ──────────────────────────────────────────────
  pvalue_transform <- reactive({
    req(ipa_funct(), input$col_pvalue)
    raw <- suppressWarnings(as.numeric(ipa_funct()[[input$col_pvalue]]))
    raw <- raw[!is.na(raw)]
    if (length(raw) == 0) return("unknown")
    if (all(raw >= 0 & raw <= 1)) "raw" else "log10"
  })
  
  output$ui_pvalue_transform_info <- renderUI({
    req(pvalue_transform())
    if (pvalue_transform() == "raw") {
      helpText("\u2714 Raw p-values detected (0-1). Will apply -log10 transformation.")
    } else if (pvalue_transform() == "log10") {
      helpText("\u2714 Already -log10 transformed values detected. Will use as-is.")
    } else {
      helpText("\u26a0 Could not determine transformation. Check column selection.")
    }
  })
  
  # ── 5. Dynamic colour inputs: one per group ─────────────────────────────
  output$ui_group_colors <- renderUI({
    req(annotation(), input$col_group)
    
    groups <- unique(annotation()[[input$col_group]])
    
    default_colors <- c("purple", "darkorange", "steelblue", "forestgreen",
                        "firebrick", "gold", "hotpink", "cyan4", "sienna", "gray40")
    
    color_inputs <- lapply(seq_along(groups), function(i) {
      textInput(
        inputId = paste0("color_group_", i),
        label   = paste0("Colour — ", groups[i], " (hex or name)"),
        value   = default_colors[i]
      )
    })
    
    tagList(color_inputs)
  })
  
  # ── 6. Sample columns: auto-detected via annotation ─────────────────────
  sample_cols_resolved <- reactive({
    req(norm_prot(), annotation(), input$col_sample_id)
    ann_samples <- annotation()[[input$col_sample_id]]
    intersect(ann_samples, colnames(norm_prot()))
  })
  
  # ── 7. Category dropdown (heatmap) ──────────────────────────────────────
  output$category_selector <- renderUI({
    req(ipa_funct(), input$col_category)
    choices <- unique(ipa_funct()[[input$col_category]])
    selectInput("selected_category",
                label    = paste0("Select: ", input$col_category),
                choices  = choices,
                selected = choices[1])
  })
  
  # ── 8. Gene list from selected category ─────────────────────────────────
  genes_interest <- reactive({
    req(ipa_funct(), input$selected_category, input$col_category, input$col_molecules)
    ipa_funct() %>%
      filter(.data[[input$col_category]] == input$selected_category) %>%
      pull(.data[[input$col_molecules]]) %>%
      strsplit(",") %>%
      unlist() %>%
      trimws() %>%
      unique()
  })
  
  # ── 9. Heatmap status ───────────────────────────────────────────────────
  output$status <- renderText({
    req(genes_interest(), norm_prot(), input$col_gene, sample_cols_resolved())
    matched <- sum(norm_prot()[[input$col_gene]] %in% genes_interest())
    paste0(
      "Category         : ", input$selected_category,                      "\n",
      "Genes in list    : ", length(genes_interest()),                      "\n",
      "Proteins matched : ", matched,                                       "\n",
      "Sample columns   : ", paste(sample_cols_resolved(), collapse = ", ")
    )
  })
  
  # ══════════════════════════════════════════════════════════════════════════
  # BAR PLOT
  # ══════════════════════════════════════════════════════════════════════════
  
  barplot_data <- reactive({
    req(ipa_funct(), input$col_category, input$col_pvalue, pvalue_transform())
    
    df        <- ipa_funct()
    raw_vals  <- suppressWarnings(as.numeric(df[[input$col_pvalue]]))
    transform <- pvalue_transform()
    
    df %>%
      mutate(
        .raw       = raw_vals,
        .neg_log10 = if (transform == "raw") -log10(.raw) else .raw
      ) %>%
      filter(!is.na(.neg_log10), is.finite(.neg_log10)) %>%
      arrange(desc(.neg_log10)) %>%
      slice_head(n = input$bar_top_n) %>%
      mutate(
        .category = factor(.data[[input$col_category]],
                           levels = rev(.data[[input$col_category]]))
      )
  })
  
  output$barplot_status <- renderText({
    req(ipa_funct(), input$col_pvalue, pvalue_transform())
    transform <- pvalue_transform()
    paste0(
      "P-value column   : ", input$col_pvalue,                            "\n",
      "Values detected  : ", if (transform == "raw")
        "raw p-values (0-1) -> -log10 applied"
      else
        "already -log10 transformed -> used as-is",                       "\n",
      "Total pathways   : ", nrow(ipa_funct()),                           "\n",
      "Showing top      : ", min(input$bar_top_n, nrow(ipa_funct()))
    )
  })
  
  make_barplot <- function() {
    req(barplot_data())
    df        <- barplot_data()
    transform <- pvalue_transform()
    bar_color <- if (is.null(input$bar_color) || input$bar_color == "") "steelblue" else input$bar_color
    x_label   <- if (transform == "raw") paste0("-log10(", input$col_pvalue, ")")
    else input$col_pvalue
    
    ggplot(df, aes(x = .neg_log10, y = .category)) +
      geom_col(fill = bar_color, width = 0.7) +
      labs(
        title = paste0("Top ", nrow(df), " Enriched Pathways / Functions"),
        x     = x_label,
        y     = NULL
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.y        = element_text(size = 10),
        axis.text.x        = element_text(size = 10),
        plot.title         = element_text(face = "bold", size = 13),
        panel.grid.major.y = element_blank()
      )
  }
  
  output$barplot <- renderPlot({
    req(barplot_data())
    make_barplot()
  }, height = function() input$bar_plot_height)
  
  # ── Download bar plot (sidebar button) ────────────────────────────────────
  output$download_barplot <- downloadHandler(
    filename = function() "Barplot_Enriched_Pathways.png",
    content  = function(file) {
      ggsave(file, plot = make_barplot(),
             width = 10, height = 8, dpi = 300, units = "in")
    }
  )
  
  # ── Download bar plot (main panel button) ─────────────────────────────────
  output$download_barplot_main <- downloadHandler(
    filename = function() "Barplot_Enriched_Pathways.png",
    content  = function(file) {
      ggsave(file, plot = make_barplot(),
             width = 10, height = 8, dpi = 300, units = "in")
    }
  )
  
  # ══════════════════════════════════════════════════════════════════════════
  # CHORD DIAGRAM
  # ══════════════════════════════════════════════════════════════════════════
  
  # ── 10. Build pathway-molecule edge table ──────────────────────────────
  chord_data <- reactive({
    req(ipa_funct(), input$col_category, input$col_molecules,
        input$col_pvalue, pvalue_transform())
    
    df        <- ipa_funct()
    raw_vals  <- suppressWarnings(as.numeric(df[[input$col_pvalue]]))
    transform <- pvalue_transform()
    
    df <- df %>%
      mutate(.neg_log10 = if (transform == "raw") -log10(raw_vals) else raw_vals) %>%
      filter(!is.na(.neg_log10), is.finite(.neg_log10)) %>%
      arrange(desc(.neg_log10)) %>%
      slice_head(n = input$chord_top_n)
    
    edges <- df %>%
      select(pathway   = all_of(input$col_category),
             molecules = all_of(input$col_molecules)) %>%
      mutate(molecules = strsplit(as.character(molecules), ",")) %>%
      unnest(molecules) %>%
      mutate(molecules = trimws(molecules)) %>%
      filter(molecules != "", !is.na(molecules))
    
    mol_freq <- edges %>%
      count(molecules, name = "freq") %>%
      filter(freq >= input$chord_min_freq)
    
    edges %>% filter(molecules %in% mol_freq$molecules)
  })
  
  # ── 11. Chord status ────────────────────────────────────────────────────
  output$chord_status <- renderText({
    req(chord_data())
    cd <- chord_data()
    paste0(
      "Pathways shown   : ", n_distinct(cd$pathway),  "\n",
      "Molecules shown  : ", n_distinct(cd$molecules), "\n",
      "Total connections: ", nrow(cd),                 "\n",
      "Min appearances  : ", input$chord_min_freq,     "\n",
      "Inner radius     : ", input$chord_inner_radius
    )
  })
  
  # ── 12. Core chord diagram function ────────────────────────────────────
  make_chord <- function() {
    cd <- chord_data()
    
    validate(
      need(nrow(cd) > 0,
           "No pathway-molecule connections to display. Try lowering 'Minimum appearances'."),
      need(n_distinct(cd$pathway) >= 2,
           "At least 2 pathways are needed for a chord diagram.")
    )
    
    mat    <- as.matrix(table(cd$pathway, cd$molecules))
    n_path <- nrow(mat)
    n_mol  <- ncol(mat)
    
    pathway_color  <- if (is.null(input$chord_pathway_color)  || input$chord_pathway_color  == "") "tomato"    else input$chord_pathway_color
    molecule_color <- if (is.null(input$chord_molecule_color) || input$chord_molecule_color == "") "steelblue" else input$chord_molecule_color
    
    path_cols <- colorRampPalette(c(pathway_color,  "lightyellow"))(n_path)
    mol_cols  <- colorRampPalette(c(molecule_color, "lightyellow"))(n_mol)
    
    grid_colors <- c(
      setNames(path_cols, rownames(mat)),
      setNames(mol_cols,  colnames(mat))
    )
    
    circos.clear()
    
    par(mar = c(2, 2, 3, 2))
    
    circos.par(
      "gap.after"    = 2,
      "start.degree" = 90
    )
    
    chordDiagram(
      mat,
      h                 = input$chord_inner_radius,
      grid.col          = grid_colors,
      transparency      = 0.4,
      annotationTrack   = "grid",
      preAllocateTracks = list(track.height = 0.08),
      directional       = FALSE
    )
    
    # ── Labels: clockwise perpendicular orientation ──────────────────────
    circos.trackPlotRegion(
      track.index = 1,
      panel.fun   = function(x, y) {
        
        sector_name <- get.cell.meta.data("sector.index")
        xlim        <- get.cell.meta.data("xlim")
        ylim        <- get.cell.meta.data("ylim")
        
        is_pathway  <- sector_name %in% rownames(mat)
        
        circos.text(
          mean(xlim),
          ylim[1] + 0.1,
          sector_name,
          facing     = "clockwise",
          niceFacing = TRUE,
          adj        = c(0, 0.5),
          cex        = if (is_pathway) input$chord_label_size_path
          else           input$chord_label_size_mol
        )
      },
      bg.border = NA
    )
    
    title(
      paste0("Chord Diagram — Top ", n_distinct(cd$pathway),
             " Pathways & Molecules"),
      cex.main = 1.1
    )
    
    legend(
      "bottomleft",
      legend = c("Pathways / Functions", "Molecules"),
      fill   = c(pathway_color, molecule_color),
      border = NA,
      bty    = "n",
      cex    = 0.85
    )
    
    circos.clear()
  }
  
  # ── 13. Render chord diagram ────────────────────────────────────────────
  output$chord_plot <- renderPlot({
    req(chord_data())
    make_chord()
  }, height = function() input$chord_plot_height)
  
  # ── Download chord diagram (sidebar button) ────────────────────────────
  output$download_chord <- downloadHandler(
    filename = function() "Chord_Diagram_Pathways_Molecules.png",
    content  = function(file) {
      png(file, width = 10, height = 10, units = "in", res = 300)
      make_chord()
      dev.off()
    }
  )
  
  # ── Download chord diagram (main panel button) ─────────────────────────
  output$download_chord_main <- downloadHandler(
    filename = function() "Chord_Diagram_Pathways_Molecules.png",
    content  = function(file) {
      png(file, width = 10, height = 10, units = "in", res = 300)
      make_chord()
      dev.off()
    }
  )
  
  # ══════════════════════════════════════════════════════════════════════════
  # HEATMAP
  # ══════════════════════════════════════════════════════════════════════════
  
  make_heatmap <- function(gene_filter, title_text) {
    
    req(input$col_gene, input$col_sample_id, input$col_group)
    
    data        <- norm_prot()
    ann         <- annotation()
    sample_cols <- sample_cols_resolved()
    
    if (!is.null(gene_filter)) {
      data <- data[data[[input$col_gene]] %in% gene_filter, ]
    }
    
    validate(
      need(nrow(data) > 1,
           "Fewer than 2 proteins matched for this category — cannot draw heatmap. Try another.")
    )
    
    mat           <- as.matrix(data[, sample_cols, drop = FALSE])
    rownames(mat) <- data[[input$col_gene]]
    mat_scaled    <- t(scale(t(mat)))
    
    col_annotation           <- data.frame(ann[[input$col_group]],
                                           row.names   = ann[[input$col_sample_id]],
                                           check.names = FALSE)
    colnames(col_annotation) <- input$col_group
    
    groups   <- unique(ann[[input$col_group]])
    grp_cols <- setNames(
      sapply(seq_along(groups), function(i) {
        val <- input[[paste0("color_group_", i)]]
        if (is.null(val) || val == "") "grey" else val
      }),
      groups
    )
    ann_colors        <- list(grp_cols)
    names(ann_colors) <- input$col_group
    
    pheatmap(
      mat_scaled,
      annotation_col    = col_annotation,
      annotation_colors = ann_colors,
      cluster_rows      = TRUE,
      cluster_cols      = FALSE,
      show_rownames     = nrow(mat) <= 100,
      show_colnames     = TRUE,
      scale             = "none",
      color             = colorRampPalette(c("darkblue", "white", "darkred"))(100),
      border_color      = NA,
      fontsize_row      = 8,
      fontsize_col      = 10,
      main              = title_text
    )
  }
  
  output$heatmap_plot <- renderPlot({
    req(norm_prot(), annotation(), genes_interest(), input$selected_category)
    make_heatmap(
      gene_filter = genes_interest(),
      title_text  = paste0(input$selected_category, "\n(Z-score, by Group)")
    )
  }, height = function() input$plot_height)
  
  # ── Download single heatmap (sidebar button) ───────────────────────────
  output$download_single <- downloadHandler(
    filename = function() {
      clean <- gsub("[^A-Za-z0-9]", "_", input$selected_category)
      paste0("Heatmap_", clean, ".png")
    },
    content = function(file) {
      png(file, width = 8, height = 10, units = "in", res = 300)
      make_heatmap(
        gene_filter = genes_interest(),
        title_text  = paste0(input$selected_category, "\n(Z-score, by Group)")
      )
      dev.off()
    }
  )
  
  # ── Download single heatmap (main panel button) ────────────────────────
  output$download_single_main <- downloadHandler(
    filename = function() {
      clean <- gsub("[^A-Za-z0-9]", "_", input$selected_category)
      paste0("Heatmap_", clean, ".png")
    },
    content = function(file) {
      png(file, width = 8, height = 10, units = "in", res = 300)
      make_heatmap(
        gene_filter = genes_interest(),
        title_text  = paste0(input$selected_category, "\n(Z-score, by Group)")
      )
      dev.off()
    }
  )
  
  # ── Download all heatmaps (sidebar button) ─────────────────────────────
  output$download_all <- downloadHandler(
    filename = function() "All_Heatmaps.zip",
    content  = function(file) {
      req(input$col_category, input$col_molecules, input$col_gene)
      
      tmp_dir <- tempdir()
      saved   <- c()
      
      for (cat_val in unique(ipa_funct()[[input$col_category]])) {
        
        genes <- ipa_funct() %>%
          filter(.data[[input$col_category]] == cat_val) %>%
          pull(.data[[input$col_molecules]]) %>%
          strsplit(",") %>%
          unlist() %>%
          trimws() %>%
          unique()
        
        if (length(genes) == 0) next
        
        data <- norm_prot()[norm_prot()[[input$col_gene]] %in% genes, ]
        if (nrow(data) <= 1) next
        
        clean    <- gsub("[^A-Za-z0-9]", "_", cat_val)
        png_path <- file.path(tmp_dir, paste0("Heatmap_", clean, ".png"))
        
        png(png_path, width = 8, height = 10, units = "in", res = 300)
        make_heatmap(gene_filter = genes,
                     title_text  = paste0(cat_val, "\n(Z-score, by Group)"))
        dev.off()
        
        saved <- c(saved, png_path)
      }
      
      zip(zipfile = file, files = saved, flags = "-j")
    }
  )
  
  # ── Download all heatmaps (main panel button) ──────────────────────────
  output$download_all_main <- downloadHandler(
    filename = function() "All_Heatmaps.zip",
    content  = function(file) {
      req(input$col_category, input$col_molecules, input$col_gene)
      
      tmp_dir <- tempdir()
      saved   <- c()
      
      for (cat_val in unique(ipa_funct()[[input$col_category]])) {
        
        genes <- ipa_funct() %>%
          filter(.data[[input$col_category]] == cat_val) %>%
          pull(.data[[input$col_molecules]]) %>%
          strsplit(",") %>%
          unlist() %>%
          trimws() %>%
          unique()
        
        if (length(genes) == 0) next
        
        data <- norm_prot()[norm_prot()[[input$col_gene]] %in% genes, ]
        if (nrow(data) <= 1) next
        
        clean    <- gsub("[^A-Za-z0-9]", "_", cat_val)
        png_path <- file.path(tmp_dir, paste0("Heatmap_", clean, ".png"))
        
        png(png_path, width = 8, height = 10, units = "in", res = 300)
        make_heatmap(gene_filter = genes,
                     title_text  = paste0(cat_val, "\n(Z-score, by Group)"))
        dev.off()
        
        saved <- c(saved, png_path)
      }
      
      zip(zipfile = file, files = saved, flags = "-j")
    }
  )
  
  # ══════════════════════════════════════════════════════════════════════════
  # BOXPLOT
  # ══════════════════════════════════════════════════════════════════════════
  
  # ── All unique molecules across all pathways ────────────────────────────
  all_molecules <- reactive({
    req(ipa_funct(), input$col_molecules)
    ipa_funct()[[input$col_molecules]] %>%
      strsplit(",") %>%
      unlist() %>%
      trimws() %>%
      unique() %>%
      sort()
  })
  
  # ── Boxplot protein selector ────────────────────────────────────────────
  output$ui_boxplot_protein <- renderUI({
    req(all_molecules())
    mols <- all_molecules()
    selectInput("boxplot_protein",
                label   = "Select protein / gene for boxplot",
                choices = mols,
                selected = mols[1])
  })
  
  # ── Boxplot data ────────────────────────────────────────────────────────
  boxplot_data <- reactive({
    req(input$boxplot_protein, norm_prot(), annotation(),
        input$col_gene, input$col_sample_id, input$col_group,
        sample_cols_resolved())
    
    protein     <- input$boxplot_protein
    data        <- norm_prot()
    ann         <- annotation()
    sample_cols <- sample_cols_resolved()
    
    row_match <- data[[input$col_gene]] == protein
    prot_row  <- data[row_match, , drop = FALSE]
    
    validate(
      need(nrow(prot_row) >= 1,
           paste0("Protein '", protein, "' was not found in the normalized data file."))
    )
    
    # Use the first matching row if duplicates exist
    prot_row <- prot_row[1, ]
    
    vals <- as.numeric(prot_row[, sample_cols, drop = TRUE])
    
    df <- data.frame(
      sample = sample_cols,
      value  = vals,
      stringsAsFactors = FALSE
    )
    
    # Merge with annotation to get group
    ann_sub <- ann[, c(input$col_sample_id, input$col_group), drop = FALSE]
    colnames(ann_sub) <- c("sample", "group")
    
    df <- merge(df, ann_sub, by = "sample", all.x = TRUE)
    df <- df[!is.na(df$value), ]
    df
  })
  
  # ── Group colours for boxplot (reuse heatmap colour inputs) ────────────
  boxplot_group_colors <- reactive({
    req(annotation(), input$col_group)
    groups <- unique(annotation()[[input$col_group]])
    setNames(
      sapply(seq_along(groups), function(i) {
        val <- input[[paste0("color_group_", i)]]
        if (is.null(val) || val == "") "grey" else val
      }),
      groups
    )
  })
  
  # ── Boxplot status line ─────────────────────────────────────────────────
  output$boxplot_status <- renderText({
    req(boxplot_data(), input$boxplot_protein)
    df <- boxplot_data()
    paste0(
      "Selected protein : ", input$boxplot_protein, "\n",
      "Samples found    : ", nrow(df),               "\n",
      "Groups           : ", paste(unique(df$group), collapse = ", ")
    )
  })
  
  # ── Core boxplot function ───────────────────────────────────────────────
  make_boxplot <- function() {
    req(boxplot_data(), boxplot_group_colors())
    
    df         <- boxplot_data()
    grp_colors <- boxplot_group_colors()
    protein    <- input$boxplot_protein
    
    # Keep only groups present in the data
    present_groups <- unique(df$group)
    grp_colors     <- grp_colors[names(grp_colors) %in% present_groups]
    
    df$group <- factor(df$group, levels = names(grp_colors))
    
    ggplot(df, aes(x = group, y = value, fill = group, colour = group)) +
      geom_boxplot(
        alpha    = 0.5,
        outlier.shape = NA,
        width    = 0.5
      ) +
      geom_jitter(
        width = 0.15,
        size  = 2.5,
        alpha = 0.85
      ) +
      scale_fill_manual(values   = grp_colors) +
      scale_colour_manual(values = grp_colors) +
      labs(
        title = paste0("Normalized Abundance — ", protein),
        x     = NULL,
        y     = "Normalized abundance / counts"
      ) +
      theme_bw(base_size = 13) +
      theme(
        legend.position    = "none",
        plot.title         = element_text(face = "bold", size = 14),
        axis.text.x        = element_text(size = 12),
        axis.text.y        = element_text(size = 11),
        panel.grid.major.x = element_blank()
      )
  }
  
  # ── Render boxplot ──────────────────────────────────────────────────────
  output$boxplot_plot <- renderPlot({
    req(boxplot_data())
    make_boxplot()
  }, height = function() input$boxplot_height)
  
  # ── Download boxplot (sidebar button) ──────────────────────────────────
  output$download_boxplot_sidebar <- downloadHandler(
    filename = function() {
      clean <- gsub("[^A-Za-z0-9]", "_", input$boxplot_protein)
      paste0("Boxplot_", clean, ".png")
    },
    content = function(file) {
      ggsave(file, plot = make_boxplot(),
             width = 8, height = 6, dpi = 300, units = "in")
    }
  )
  
  # ── Download boxplot (main panel button) ───────────────────────────────
  output$download_boxplot_main <- downloadHandler(
    filename = function() {
      clean <- gsub("[^A-Za-z0-9]", "_", input$boxplot_protein)
      paste0("Boxplot_", clean, ".png")
    },
    content = function(file) {
      ggsave(file, plot = make_boxplot(),
             width = 8, height = 6, dpi = 300, units = "in")
    }
  )
  
}

# ══════════════════════════════════════════════════════════════════════════════
shinyApp(ui = ui, server = server)