library(shiny)
library(DT)
library(ggplot2)
library(gridExtra)

# ============================================================================
# CONSTANTS
# ============================================================================
METRIC_COLS <- c("RPBZ", "MQBZ", "BQBZ", "MQSBZ", "SCBZ", "MQ0F", "MQ", "QUAL", "DP")
RAINBOW     <- rainbow(length(METRIC_COLS))

SAMPLE_HIST_COLS <- c("sample_DP", "SP")
SAMPLE_AD_COLS   <- c("ADREF_pct", "ADALT_pct")
SAMPLE_PL_COLS   <- c("PL00", "PL01", "PL11")
SAMPLE_FILTER_COLS <- c(SAMPLE_HIST_COLS, SAMPLE_AD_COLS)
SAMPLE_ALL_COLS    <- c(SAMPLE_HIST_COLS, SAMPLE_AD_COLS, SAMPLE_PL_COLS)
SAMPLE_RAINBOW     <- rainbow(length(SAMPLE_ALL_COLS))

GT_CLASSES <- c("0/0", "0/1", "1/1")
GT_IDS     <- c("00", "01", "11")
GT_BCF     <- c("00" = "RR", "01" = "het", "11" = "AA")

SAMPLE_LABELS <- c(
  "sample_DP"  = "Sample DP",
  "SP"         = "SP",
  "ADREF_pct"  = "ADREF %",
  "ADALT_pct"  = "ADALT %",
  "PL00"       = "PL(0/0)",
  "PL01"       = "PL(0/1)",
  "PL11"       = "PL(1/1)"
)

match_gt <- function(gt_vals, gt_class) {
  gt_vals %in% c(gt_class, gsub("/", "|", gt_class, fixed = TRUE))
}

# ============================================================================
# UI
# ============================================================================
ui <- fluidPage(
  tags$head(tags$script(HTML("
    function copyCmd(elId, btnId) {
      var el = document.getElementById(elId);
      if (!el) return;
      var text = el.innerText || el.textContent;
      text = text.replace(/^\\s+|\\s+$/g, '');
      if (navigator.clipboard && window.isSecureContext) {
        navigator.clipboard.writeText(text).then(function() {
          flashBtn(btnId);
        }, function() {
          fallbackCopy(text, btnId);
        });
      } else {
        fallbackCopy(text, btnId);
      }
    }
    function fallbackCopy(text, btnId) {
      var ta = document.createElement('textarea');
      ta.value = text;
      ta.style.position = 'fixed';
      ta.style.top = '0';
      ta.style.left = '0';
      ta.style.width = '2em';
      ta.style.height = '2em';
      ta.style.padding = '0';
      ta.style.border = 'none';
      ta.style.outline = 'none';
      ta.style.boxShadow = 'none';
      ta.style.background = 'transparent';
      document.body.appendChild(ta);
      ta.focus();
      ta.select();
      try {
        document.execCommand('copy');
        flashBtn(btnId);
      } catch(e) {
        alert('Copy failed. Please select and copy manually.');
      }
      document.body.removeChild(ta);
    }
    function flashBtn(btnId) {
      var btn = document.getElementById(btnId);
      if (!btn) return;
      var orig = btn.innerHTML;
      btn.innerHTML = '\\u2713 Copied!';
      btn.style.color = '#27ae60';
      setTimeout(function() { btn.innerHTML = orig; btn.style.color = ''; }, 2000);
    }
  "))),

  titlePanel("BCFtools VCF Metrics Viewer"),

  tabsetPanel(
    id = "main_tabs",

    # ================================================================
    # TAB 1 - Site-Level Metrics
    # ================================================================
    tabPanel(
      "Site-Level Metrics",
      br(),
      fluidRow(
        column(
          width = 3,
          wellPanel(
            h5(tags$b("Pre-filtering Options")),
            checkboxInput("biallelic_only", "Keep only biallelic sites", FALSE),
            radioButtons("variant_type", "Variant type to keep:",
                         choices  = c("All variants" = "all",
                                      "SNPs only"    = "snps",
                                      "INDELs only"  = "indels"),
                         selected = "all", inline = TRUE)
          ),
          wellPanel(
            h5(tags$b("BCFtools Filter Command")),
            verbatimTextOutput("bcf_command"),
            tags$button(
              id = "copy_filter_btn", class = "btn btn-default btn-sm",
              onclick = "copyCmd('bcf_command', 'copy_filter_btn')",
              tags$i(class = "fa fa-clipboard"), " Copy"
            )
          ),
          h4("Data Preview"),
          DTOutput("data_table"),
          br(),
          downloadButton("download_png", "Export Plot as PNG",
                          style = "margin-top:12px;")
        ),
        column(
          width = 9,
          wellPanel(
            fluidRow(
              column(
                4,
                textInput("vcf_name", "VCF Filename",
                          value = "", placeholder = "e.g. snp_final.vcf.gz")
              ),
              column(
                8,
                h5(tags$b("Step 1: Run these commands in your terminal")),
                verbatimTextOutput("site_query_cmd"),
                tags$button(
                  id = "copy_site_query_btn", class = "btn btn-default btn-sm",
                  onclick = "copyCmd('site_query_cmd', 'copy_site_query_btn')",
                  tags$i(class = "fa fa-clipboard"), " Copy"
                )
              )
            ),
            fileInput("file", "Step 2: Upload the resulting subsampled TSV file",
                      accept = c(".tsv", ".txt", ".tab"))
          ),
          uiOutput("hist_grid_ui")
        )
      )
    ),

    # ================================================================
    # TAB 2 - Per-Sample Metrics
    # ================================================================
    tabPanel(
      "Per-Sample Metrics",
      br(),
      fluidRow(
        column(
          width = 3,
          wellPanel(
            h5(tags$b("BCFtools Filter Command")),
            verbatimTextOutput("sample_bcf_command"),
            tags$button(
              id = "copy_sample_filter_btn", class = "btn btn-default btn-sm",
              onclick = "copyCmd('sample_bcf_command', 'copy_sample_filter_btn')",
              tags$i(class = "fa fa-clipboard"), " Copy"
            )
          ),
          h4("Data Preview"),
          DTOutput("sample_data_table"),
          br(),
          downloadButton("download_sample_png", "Export Plots as PNG",
                          style = "margin-top:12px;")
        ),
        column(
          width = 9,
          wellPanel(
            fluidRow(
              column(
                4,
                textInput("sample_vcf_name", "VCF Filename",
                          value = "", placeholder = "e.g. snp_final.vcf.gz")
              ),
              column(
                8,
                h5(tags$b("Step 1: Run these commands in your terminal")),
                verbatimTextOutput("sample_query_cmd"),
                tags$button(
                  id = "copy_sample_query_btn", class = "btn btn-default btn-sm",
                  onclick = "copyCmd('sample_query_cmd', 'copy_sample_query_btn')",
                  tags$i(class = "fa fa-clipboard"), " Copy"
                )
              )
            ),
            fileInput("sample_file", "Step 2: Upload the resulting subsampled TSV file",
                      accept = c(".tsv", ".txt", ".tab"))
          ),
          wellPanel(
            fluidRow(
              column(6,
                textInput("pct_missing", "Max % Missing Across Samples",
                          value = "", placeholder = "e.g. 20")
              ),
              column(6,
                tags$div(style = "margin-top: 25px; font-size: 12px; color: #555;",
                  uiOutput("pct_missing_hint")
                )
              )
            )
          ),
          uiOutput("sample_grid_ui")
        )
      )
    ),

    # ================================================================
    # TAB 3 - PL Viewer (view-only, no filtering)
    # ================================================================
    tabPanel(
      "PL Viewer",
      br(),
      conditionalPanel(
        condition = "!output.has_sample_data",
        tags$div(
          style = "color:#888; font-style:italic; padding:40px; text-align:center; font-size:16px;",
          "Upload per-sample data on the ",
          tags$b("Per-Sample Metrics"),
          " tab to view PL distributions."
        )
      ),
      conditionalPanel(
        condition = "output.has_sample_data",
        fluidRow(
          column(4, plotOutput("pl_plot_PL00", height = "400px")),
          column(4, plotOutput("pl_plot_PL01", height = "400px")),
          column(4, plotOutput("pl_plot_PL11", height = "400px"))
        ),
        br(),
        downloadButton("download_pl_png", "Export PL Plots as PNG",
                        style = "margin-top:4px; margin-bottom:12px;")
      )
    ),

    # ================================================================
    # TAB 4 - Filter Summary
    # ================================================================
    tabPanel(
      "Filter Summary",
      br(),
      fluidRow(
        column(
          width = 6,
          h3("Site-Level Filter Summary"),
          hr(),
          conditionalPanel(
            condition = "output.has_site_summary",
            DTOutput("filter_summary_table"),
            br(),
            downloadButton("download_summary",
                            "Export Site-Level Summary as TSV",
                            style = "margin-top:4px; margin-bottom:12px;")
          ),
          conditionalPanel(
            condition = "!output.has_site_summary",
            tags$div(
              style = "color:#888; font-style:italic; padding:30px;",
              "No site-level filters have been set. ",
              "Set min/max thresholds on the ",
              tags$b("Site-Level Metrics"),
              " tab to see the filter summary here."
            )
          )
        ),
        column(
          width = 6,
          h3("Per-Sample Filter Summary"),
          hr(),
          conditionalPanel(
            condition = "output.has_sample_summary",
            DTOutput("sample_filter_summary_table"),
            br(),
            downloadButton("download_sample_summary",
                            "Export Per-Sample Summary as TSV",
                            style = "margin-top:4px; margin-bottom:12px;")
          ),
          conditionalPanel(
            condition = "!output.has_sample_summary",
            tags$div(
              style = "color:#888; font-style:italic; padding:30px;",
              "No per-sample filters have been set. ",
              "Set min/max thresholds on the ",
              tags$b("Per-Sample Metrics"),
              " tab to see the filter summary here."
            )
          )
        )
      )
    )
  )
)

# ============================================================================
# SERVER
# ============================================================================
server <- function(input, output, session) {

  # ==========================================================================
  # Site-level query command
  # ==========================================================================
  output$site_query_cmd <- renderText({
    vcf <- input$vcf_name
    if (is.null(vcf) || vcf == "") return("Enter VCF filename above to generate the query command.")
    base <- sub("\\.(vcf\\.gz|vcf|bcf)$", "", vcf)
    tsv_name  <- paste0(base, "_info_metrics.tsv")
    sub_name  <- paste0("subsample_", base, "_info_metrics.tsv")
    paste0(
      "bcftools query \\\n",
      "  -H \\\n",
      "  -e 'ALT=\".\"' \\\n",
      "  -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/DP\\t%INFO/DP4\\t%INFO/RPBZ\\t%INFO/MQBZ\\t%INFO/BQBZ\\t%INFO/MQSBZ\\t%INFO/SCBZ\\t%INFO/MQ0F\\t%INFO/AD\\t%INFO/MQ\\n' \\\n",
      "  ", vcf, " > ", tsv_name, "\n\n",
      "head -n 1 ", tsv_name, " > ", sub_name, "\n",
      "tail -n +2 ", tsv_name, " | shuf -n 10000 >> ", sub_name
    )
  })

  # ==========================================================================
  # Per-sample query command
  # ==========================================================================
  output$sample_query_cmd <- renderText({
    vcf <- input$sample_vcf_name
    if (is.null(vcf) || vcf == "") return("Enter VCF filename above to generate the query command.")
    base <- sub("\\.(vcf\\.gz|vcf|bcf)$", "", vcf)
    tsv_name  <- paste0(base, "_info_metrics_persample.tsv")
    sub_name  <- paste0("subsample_", base, "_info_metrics_persample.tsv")
    paste0(
      "bcftools query \\\n",
      "  -H \\\n",
      "  -e 'ALT=\".\"' \\\n",
      "  -f '[%CHROM\\t%POS\\t%SAMPLE\\t%GT\\t%PL\\t%DP\\t%SP\\t%AD\\n]' \\\n",
      "  ", vcf, " > ", tsv_name, "\n\n",
      "head -n 1 ", tsv_name, " > ", sub_name, "\n",
      "tail -n +2 ", tsv_name, " | shuf -n 10000 >> ", sub_name
    )
  })

  # ==========================================================================
  # TAB 1 - Site-level data
  # ==========================================================================
  vcf_data <- reactive({
    req(input$file)
    raw <- readLines(input$file$datapath)
    header_idx  <- max(grep("^#", raw))
    header_line <- sub("^#\\s*", "", raw[header_idx])
    col_names   <- unlist(strsplit(header_line, "\t"))
    col_names   <- gsub("\\[\\d+\\]", "", col_names)
    col_names   <- gsub("^%", "", col_names)
    col_names   <- gsub("^INFO/", "", col_names)
    col_names   <- trimws(col_names)
    data_lines <- raw[(header_idx + 1):length(raw)]
    data_lines <- data_lines[data_lines != ""]
    df <- read.delim(
      textConnection(paste(data_lines, collapse = "\n")),
      header = FALSE, sep = "\t",
      stringsAsFactors = FALSE, na.strings = c(".", "NA", "")
    )
    if (ncol(df) == length(col_names)) colnames(df) <- col_names
    else colnames(df) <- paste0("V", seq_len(ncol(df)))
    for (col in METRIC_COLS) {
      if (col %in% colnames(df)) df[[col]] <- as.numeric(df[[col]])
    }
    df
  })

  original_snp_count <- reactive({ req(vcf_data()); nrow(vcf_data()) })

  # ---- Site histogram grid UI ----
  output$hist_grid_ui <- renderUI({
    req(vcf_data())
    available <- METRIC_COLS[METRIC_COLS %in% colnames(vcf_data())]
    rows <- list()
    for (start in seq(1, length(available), by = 3)) {
      end  <- min(start + 2, length(available))
      cols <- lapply(start:end, function(i) {
        metric <- available[i]
        column(4,
          plotOutput(paste0("plot_", metric), height = "220px"),
          fluidRow(
            column(6, textInput(paste0("min_", metric), paste("Min", metric),
                                value = "", placeholder = "e.g. -2.5")),
            column(6, textInput(paste0("max_", metric), paste("Max", metric),
                                value = "", placeholder = "e.g. 2.5"))
          ),
          uiOutput(paste0("counts_", metric)), br()
        )
      })
      rows <- c(rows, list(fluidRow(cols)))
    }
    do.call(tagList, rows)
  })

  # ---- Render site-level plots + counts ----
  observe({
    req(vcf_data())
    df <- vcf_data()
    available <- METRIC_COLS[METRIC_COLS %in% colnames(df)]
    lapply(seq_along(available), function(i) {
      metric <- available[i]
      colour_idx <- which(METRIC_COLS == metric)
      output[[paste0("plot_", metric)]] <- renderPlot({
        plot_df <- df[!is.na(df[[metric]]), , drop = FALSE]
        p <- ggplot(plot_df, aes(x = .data[[metric]])) +
          geom_histogram(bins = 50, fill = RAINBOW[colour_idx],
                         color = "grey25", alpha = 0.85) +
          labs(title = metric, x = metric, y = "Count") +
          theme_minimal(base_size = 11) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        min_val <- suppressWarnings(as.numeric(input[[paste0("min_", metric)]]))
        if (!is.null(min_val) && !is.na(min_val))
          p <- p + geom_vline(xintercept = min_val, linetype = "dashed",
                              color = "black", linewidth = 0.9) +
            annotate("label", x = min_val, y = Inf, label = as.character(min_val),
                     vjust = 1.5, size = 3.4, fontface = "bold",
                     fill = "white", color = "black",
                     label.padding = unit(0.25, "lines"))
        max_val <- suppressWarnings(as.numeric(input[[paste0("max_", metric)]]))
        if (!is.null(max_val) && !is.na(max_val))
          p <- p + geom_vline(xintercept = max_val, linetype = "dashed",
                              color = "black", linewidth = 0.9) +
            annotate("label", x = max_val, y = Inf, label = as.character(max_val),
                     vjust = 1.5, size = 3.4, fontface = "bold",
                     fill = "white", color = "black",
                     label.padding = unit(0.25, "lines"))
        p
      })

      output[[paste0("counts_", metric)]] <- renderUI({
        total_snps <- original_snp_count()
        vals <- df[[metric]]; vals <- vals[!is.na(vals)]
        min_val <- suppressWarnings(as.numeric(input[[paste0("min_", metric)]]))
        max_val <- suppressWarnings(as.numeric(input[[paste0("max_", metric)]]))
        has_min <- !is.null(min_val) && !is.na(min_val)
        has_max <- !is.null(max_val) && !is.na(max_val)
        if (!has_min && !has_max) return(NULL)
        below_n <- if (has_min) sum(vals < min_val) else NA
        above_n <- if (has_max) sum(vals > max_val) else NA
        pct_fmt <- function(n) sprintf("%.1f%%", 100 * n / total_snps)
        if (has_min && has_max) {
          within_n <- sum(vals >= min_val & vals <= max_val)
          removed_n <- nrow(df) - within_n
          tags$div(style = "background:#f7f7f7; border:1px solid #ddd; border-radius:4px; padding:6px 10px; font-size:12px;",
            tags$b("Total SNPs: "), format(total_snps, big.mark = ","), tags$br(),
            tags$span(style = "color:#c0392b;",
              tags$b("Below min: "), format(below_n, big.mark = ","),
              paste0(" (", pct_fmt(below_n), ")")), tags$br(),
            tags$span(style = "color:#c0392b;",
              tags$b("Above max: "), format(above_n, big.mark = ","),
              paste0(" (", pct_fmt(above_n), ")")), tags$br(),
            tags$span(style = "color:#27ae60;",
              tags$b("In range: "), format(within_n, big.mark = ","),
              paste0(" (", pct_fmt(within_n), ")")), tags$br(),
            tags$span(style = "color:#e67e22;",
              tags$b("Removed: "), format(removed_n, big.mark = ","),
              paste0(" (", pct_fmt(removed_n), ")"))
          )
        } else if (has_min) {
          kept_n <- sum(vals >= min_val)
          tags$div(style = "background:#f7f7f7; border:1px solid #ddd; border-radius:4px; padding:6px 10px; font-size:12px;",
            tags$b("Total SNPs: "), format(total_snps, big.mark = ","), tags$br(),
            tags$span(style = "color:#c0392b;",
              tags$b("Below min: "), format(below_n, big.mark = ","),
              paste0(" (", pct_fmt(below_n), ")")), tags$br(),
            tags$span(style = "color:#27ae60;",
              tags$b("Kept (\u2265 min): "), format(kept_n, big.mark = ","),
              paste0(" (", pct_fmt(kept_n), ")"))
          )
        } else {
          kept_n <- sum(vals <= max_val)
          tags$div(style = "background:#f7f7f7; border:1px solid #ddd; border-radius:4px; padding:6px 10px; font-size:12px;",
            tags$b("Total SNPs: "), format(total_snps, big.mark = ","), tags$br(),
            tags$span(style = "color:#c0392b;",
              tags$b("Above max: "), format(above_n, big.mark = ","),
              paste0(" (", pct_fmt(above_n), ")")), tags$br(),
            tags$span(style = "color:#27ae60;",
              tags$b("Kept (\u2264 max): "), format(kept_n, big.mark = ","),
              paste0(" (", pct_fmt(kept_n), ")"))
          )
        }
      })
    })
  })

  output$data_table <- renderDT({
    req(vcf_data())
    datatable(vcf_data(), options = list(scrollX = TRUE, pageLength = 5),
              class = "compact cell-border stripe")
  })

  # ---- Site-level BCFtools filter command ----
  site_filter_info <- reactive({
    df <- tryCatch(vcf_data(), error = function(e) NULL)
    available <- if (!is.null(df)) METRIC_COLS[METRIC_COLS %in% colnames(df)]
                 else METRIC_COLS
    filters <- list(); name_parts <- character(0)
    for (metric in available) {
      min_val <- suppressWarnings(as.numeric(input[[paste0("min_", metric)]]))
      max_val <- suppressWarnings(as.numeric(input[[paste0("max_", metric)]]))
      has_min <- !is.null(min_val) && !is.na(min_val)
      has_max <- !is.null(max_val) && !is.na(max_val)
      if (!has_min && !has_max) next
      prefix <- if (metric %in% c("QUAL", "DP")) "" else "INFO/"
      if (has_min) filters <- c(filters, paste0(prefix, metric, " >= ", min_val))
      if (has_max) filters <- c(filters, paste0(prefix, metric, " <= ", max_val))
      tag <- metric
      if (has_min && has_max) tag <- paste0(tag, "min", min_val, "max", max_val)
      else if (has_min)       tag <- paste0(tag, "min", min_val)
      else                    tag <- paste0(tag, "max", max_val)
      name_parts <- c(name_parts, tag)
    }
    list(filters = filters, name_parts = name_parts)
  })

  output$bcf_command <- renderText({
    info <- site_filter_info()
    vcf  <- input$vcf_name
    if (vcf == "") return("Set thresholds and VCF name above.")

    base <- sub("\\.(vcf\\.gz|vcf|bcf)$", "", vcf)
    ext  <- sub(paste0("^", gsub("\\.", "\\\\.", base)), "", vcf)
    flag <- if (grepl("\\.vcf\\.gz$", vcf)) "-Oz" else if (grepl("\\.bcf$", vcf)) "-Ob" else "-Ov"

    view_flags <- c()
    view_name_parts <- c()
    if (isTRUE(input$biallelic_only)) {
      view_flags <- c(view_flags, "-m2 -M2")
      view_name_parts <- c(view_name_parts, "biallelic")
    }
    if (!is.null(input$variant_type) && input$variant_type != "all") {
      view_flags <- c(view_flags, paste0("-v ", input$variant_type))
      view_name_parts <- c(view_name_parts, input$variant_type)
    }

    has_view    <- length(view_flags) > 0
    has_filter  <- length(info$filters) > 0

    if (!has_view && !has_filter) return("Set thresholds or pre-filtering options above.")

    all_name_parts <- c(view_name_parts, info$name_parts)
    out <- paste0(base, "_", paste(all_name_parts, collapse = "_"), ext)

    if (has_view && has_filter) {
      expr <- paste(info$filters, collapse = " && ")
      cmd <- paste0(
        "bcftools view ", paste(view_flags, collapse = " "), " ", vcf,
        " | \\\nbcftools filter \\\n  -i '", expr, "' \\\n  ",
        flag, " -o ", out
      )
    } else if (has_view) {
      cmd <- paste0(
        "bcftools view ", paste(view_flags, collapse = " "), " \\\n  ",
        vcf, " \\\n  ", flag, " -o ", out
      )
    } else {
      expr <- paste(info$filters, collapse = " && ")
      cmd <- paste0(
        "bcftools filter \\\n  -i '", expr, "' \\\n  ", vcf,
        " \\\n  ", flag, " -o ", out
      )
    }

    paste0(cmd,
           " && \\\nbcftools index ", out,
           " && \\\necho \"Filtered SNP count:\" && \\\nbcftools view ", out,
           " | grep -v '^#' | wc -l")
  })

  # ---- Site-level filter summary ----
  site_summary_df <- reactive({
    rows <- list()
    if (isTRUE(input$biallelic_only)) {
      rows <- c(rows, list(data.frame(
        Metric = "Biallelic", Range = "biallelic only",
        Total = NA_integer_, Kept = NA_integer_, Removed = NA_integer_,
        Pct_Kept = "-", Pct_Removed = "-",
        stringsAsFactors = FALSE
      )))
    }
    if (!is.null(input$variant_type) && input$variant_type != "all") {
      rows <- c(rows, list(data.frame(
        Metric = "Variant Type", Range = input$variant_type,
        Total = NA_integer_, Kept = NA_integer_, Removed = NA_integer_,
        Pct_Kept = "-", Pct_Removed = "-",
        stringsAsFactors = FALSE
      )))
    }

    df <- tryCatch(vcf_data(), error = function(e) NULL)
    if (!is.null(df)) {
      total_snps <- nrow(df)
      available <- METRIC_COLS[METRIC_COLS %in% colnames(df)]
      for (metric in available) {
        min_val <- suppressWarnings(as.numeric(input[[paste0("min_", metric)]]))
        max_val <- suppressWarnings(as.numeric(input[[paste0("max_", metric)]]))
        has_min <- !is.null(min_val) && !is.na(min_val)
        has_max <- !is.null(max_val) && !is.na(max_val)
        if (!has_min && !has_max) next
        vals <- df[[metric]]
        range_str <- paste0(
          if (has_min) paste0("min=", min_val) else "",
          if (has_min && has_max) ", " else "",
          if (has_max) paste0("max=", max_val) else ""
        )
        if (has_min && has_max) kept_n <- sum(vals >= min_val & vals <= max_val, na.rm = TRUE)
        else if (has_min)       kept_n <- sum(vals >= min_val, na.rm = TRUE)
        else                    kept_n <- sum(vals <= max_val, na.rm = TRUE)
        removed_n <- total_snps - kept_n
        rows <- c(rows, list(data.frame(
          Metric = metric, Range = range_str, Total = total_snps,
          Kept = kept_n, Removed = removed_n,
          Pct_Kept = sprintf("%.1f%%", 100 * kept_n / total_snps),
          Pct_Removed = sprintf("%.1f%%", 100 * removed_n / total_snps),
          stringsAsFactors = FALSE
        )))
      }
    }
    if (length(rows) == 0) return(data.frame(Metric = character(0)))
    do.call(rbind, rows)
  })

  output$has_site_summary <- reactive({
    df <- site_summary_df()
    !is.null(df) && nrow(df) > 0
  })
  outputOptions(output, "has_site_summary", suspendWhenHidden = FALSE)

  output$filter_summary_table <- renderDT({
    req(nrow(site_summary_df()) > 0)
    datatable(site_summary_df(), options = list(dom = "t", paging = FALSE),
              rownames = FALSE, class = "compact cell-border stripe")
  })

  output$download_summary <- downloadHandler(
    filename = function() paste0("filter_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv"),
    content  = function(file) write.table(site_summary_df(), file, sep = "\t", row.names = FALSE, quote = FALSE)
  )

  output$download_png <- downloadHandler(
    filename = function() paste0("site_metrics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content  = function(file) {
      df <- vcf_data()
      available <- METRIC_COLS[METRIC_COLS %in% colnames(df)]
      plots <- lapply(seq_along(available), function(i) {
        metric <- available[i]
        colour_idx <- which(METRIC_COLS == metric)
        plot_df <- df[!is.na(df[[metric]]), , drop = FALSE]
        p <- ggplot(plot_df, aes(x = .data[[metric]])) +
          geom_histogram(bins = 50, fill = RAINBOW[colour_idx], color = "grey25", alpha = 0.85) +
          labs(title = metric, x = metric, y = "Count") +
          theme_minimal(base_size = 11) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        min_val <- suppressWarnings(as.numeric(input[[paste0("min_", metric)]]))
        if (!is.null(min_val) && !is.na(min_val))
          p <- p + geom_vline(xintercept = min_val, linetype = "dashed", color = "black", linewidth = 0.9) +
            annotate("label", x = min_val, y = Inf, label = as.character(min_val),
                     vjust = 1.5, size = 3.4, fontface = "bold", fill = "white",
                     label.padding = unit(0.25, "lines"))
        max_val <- suppressWarnings(as.numeric(input[[paste0("max_", metric)]]))
        if (!is.null(max_val) && !is.na(max_val))
          p <- p + geom_vline(xintercept = max_val, linetype = "dashed", color = "black", linewidth = 0.9) +
            annotate("label", x = max_val, y = Inf, label = as.character(max_val),
                     vjust = 1.5, size = 3.4, fontface = "bold", fill = "white",
                     label.padding = unit(0.25, "lines"))
        p
      })
      png(file, width = 1800, height = 1500, res = 150)
      do.call(grid.arrange, c(plots, ncol = 3))
      dev.off()
    }
  )

  # ==========================================================================
  # TAB 2 - Per-sample data
  # ==========================================================================
  sample_data <- reactive({
    req(input$sample_file)
    raw <- readLines(input$sample_file$datapath)
    header_idx  <- max(grep("^#", raw))
    header_line <- sub("^#\\s*", "", raw[header_idx])
    col_names   <- unlist(strsplit(header_line, "\t"))
    col_names   <- gsub("\\[\\d+\\]", "", col_names)
    col_names   <- gsub("^%", "", col_names)
    col_names   <- gsub("^FORMAT/", "", col_names)
    col_names   <- gsub("^INFO/", "", col_names)
    col_names   <- trimws(col_names)
    data_lines <- raw[(header_idx + 1):length(raw)]
    data_lines <- data_lines[data_lines != ""]
    df <- read.delim(
      textConnection(paste(data_lines, collapse = "\n")),
      header = FALSE, sep = "\t",
      stringsAsFactors = FALSE, na.strings = c(".", "NA", "")
    )
    if (ncol(df) == length(col_names)) colnames(df) <- col_names
    else colnames(df) <- paste0("V", seq_len(ncol(df)))
    if ("DP" %in% colnames(df)) colnames(df)[colnames(df) == "DP"] <- "sample_DP"

    # Split PL
    if ("PL" %in% colnames(df)) {
      pl_split <- strsplit(as.character(df$PL), ",")
      df$PL00 <- suppressWarnings(as.numeric(sapply(pl_split, function(x) if (length(x) >= 1) x[1] else NA)))
      df$PL01 <- suppressWarnings(as.numeric(sapply(pl_split, function(x) if (length(x) >= 2) x[2] else NA)))
      df$PL11 <- suppressWarnings(as.numeric(sapply(pl_split, function(x) if (length(x) >= 3) x[3] else NA)))
    }

    # Split AD
    if ("AD" %in% colnames(df)) {
      ad_split <- strsplit(as.character(df$AD), ",")
      df$ADREF <- suppressWarnings(as.numeric(sapply(ad_split, function(x) if (length(x) >= 1) x[1] else NA)))
      df$ADALT <- suppressWarnings(as.numeric(sapply(ad_split, function(x) if (length(x) >= 2) x[2] else NA)))
      ad_total <- df$ADREF + df$ADALT
      ad_total[ad_total == 0] <- NA
      df$ADREF_pct <- df$ADREF / ad_total
      df$ADALT_pct <- df$ADALT / ad_total
    }

    for (col in c("sample_DP", "SP", "ADREF", "ADALT", "ADREF_pct", "ADALT_pct",
                  "PL00", "PL01", "PL11")) {
      if (col %in% colnames(df)) df[[col]] <- as.numeric(df[[col]])
    }
    df
  })

  original_sample_count <- reactive({ req(sample_data()); nrow(sample_data()) })

  output$has_sample_data <- reactive({
    !is.null(tryCatch(sample_data(), error = function(e) NULL))
  })
  outputOptions(output, "has_sample_data", suspendWhenHidden = FALSE)

  # ---- Percent missing hint ----
  output$pct_missing_hint <- renderUI({
    pct <- suppressWarnings(as.numeric(input$pct_missing))
    if (is.null(pct) || is.na(pct)) return(NULL)
    frac <- pct / 100
    tags$span(style = "color:#2980b9;",
      tags$b("F_MISSING < ", sprintf("%.2f", frac)),
      " \u2014 sites with \u2265 ", pct, "% missing genotypes will be removed."
    )
  })

  # ---- Per-sample grid UI (2 columns wide) ----
  output$sample_grid_ui <- renderUI({
    req(sample_data())
    df <- sample_data()
    available <- SAMPLE_FILTER_COLS[SAMPLE_FILTER_COLS %in% colnames(df)]
    rows <- list()

    for (start in seq(1, length(available), by = 2)) {
      end  <- min(start + 1, length(available))
      cols <- lapply(start:end, function(i) {
        metric <- available[i]
        is_ad  <- metric %in% SAMPLE_AD_COLS

        if (is_ad) {
          gt_inputs <- lapply(seq_along(GT_CLASSES), function(g) {
            gt_id <- GT_IDS[g]
            gt_label <- GT_CLASSES[g]
            fluidRow(
              column(12, tags$label(paste0("GT = ", gt_label), style = "font-size:11px; margin-bottom:2px;")),
              column(6, textInput(paste0("smin_", metric, "_", gt_id),
                                  paste("Min"), value = "", placeholder = "e.g. 0.2")),
              column(6, textInput(paste0("smax_", metric, "_", gt_id),
                                  paste("Max"), value = "", placeholder = "e.g. 0.8"))
            )
          })
          column(6,
            plotOutput(paste0("splot_", metric), height = "280px"),
            gt_inputs,
            uiOutput(paste0("scounts_", metric)), br()
          )
        } else {
          column(6,
            plotOutput(paste0("splot_", metric), height = "260px"),
            fluidRow(
              column(6, textInput(paste0("smin_", metric), paste("Min", SAMPLE_LABELS[metric]),
                                  value = "", placeholder = "e.g. 0")),
              column(6, textInput(paste0("smax_", metric), paste("Max", SAMPLE_LABELS[metric]),
                                  value = "", placeholder = "e.g. 100"))
            ),
            uiOutput(paste0("scounts_", metric)), br()
          )
        }
      })
      rows <- c(rows, list(fluidRow(cols)))
    }
    do.call(tagList, rows)
  })

  # ---- Render per-sample plots and counts ----
  observe({
    req(sample_data())
    df <- sample_data()
    available <- SAMPLE_FILTER_COLS[SAMPLE_FILTER_COLS %in% colnames(df)]

    lapply(seq_along(available), function(i) {
      metric  <- available[i]
      col_idx <- which(SAMPLE_ALL_COLS == metric)
      fill_col <- SAMPLE_RAINBOW[col_idx]
      is_ad   <- metric %in% SAMPLE_AD_COLS

      output[[paste0("splot_", metric)]] <- renderPlot({
        plot_df <- df[!is.na(df[[metric]]), , drop = FALSE]

        if (is_ad) {
          plot_df$GT <- as.factor(plot_df$GT)
          p <- ggplot(plot_df, aes(x = GT, y = .data[[metric]])) +
            geom_boxplot(fill = fill_col, color = "grey25", alpha = 0.85,
                         outlier.size = 0.5, outlier.alpha = 0.3) +
            scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                               limits = c(0, 1)) +
            labs(title = paste(SAMPLE_LABELS[metric], "by Genotype"),
                 x = "Genotype (GT)", y = SAMPLE_LABELS[metric]) +
            theme_minimal(base_size = 11) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))

          gt_levels <- levels(plot_df$GT)
          for (g in seq_along(GT_CLASSES)) {
            gt_id <- GT_IDS[g]
            gt_class <- GT_CLASSES[g]
            gt_match <- which(gt_levels %in% c(gt_class, gsub("/", "|", gt_class, fixed = TRUE)))
            if (length(gt_match) == 0) next
            x_pos <- gt_match[1]
            min_val <- suppressWarnings(as.numeric(input[[paste0("smin_", metric, "_", gt_id)]]))
            max_val <- suppressWarnings(as.numeric(input[[paste0("smax_", metric, "_", gt_id)]]))
            if (!is.null(min_val) && !is.na(min_val))
              p <- p + annotate("segment", x = x_pos - 0.4, xend = x_pos + 0.4,
                                y = min_val, yend = min_val,
                                linetype = "dashed", color = "red", linewidth = 0.9)
            if (!is.null(max_val) && !is.na(max_val))
              p <- p + annotate("segment", x = x_pos - 0.4, xend = x_pos + 0.4,
                                y = max_val, yend = max_val,
                                linetype = "dashed", color = "red", linewidth = 0.9)
          }
          p
        } else {
          p <- ggplot(plot_df, aes(x = .data[[metric]])) +
            geom_histogram(bins = 50, fill = fill_col, color = "grey25", alpha = 0.85) +
            labs(title = SAMPLE_LABELS[metric], x = SAMPLE_LABELS[metric], y = "Count") +
            theme_minimal(base_size = 11) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
          min_val <- suppressWarnings(as.numeric(input[[paste0("smin_", metric)]]))
          if (!is.null(min_val) && !is.na(min_val))
            p <- p + geom_vline(xintercept = min_val, linetype = "dashed",
                                color = "black", linewidth = 0.9) +
              annotate("label", x = min_val, y = Inf, label = as.character(min_val),
                       vjust = 1.5, size = 3.4, fontface = "bold",
                       fill = "white", color = "black",
                       label.padding = unit(0.25, "lines"))
          max_val <- suppressWarnings(as.numeric(input[[paste0("smax_", metric)]]))
          if (!is.null(max_val) && !is.na(max_val))
            p <- p + geom_vline(xintercept = max_val, linetype = "dashed",
                                color = "black", linewidth = 0.9) +
              annotate("label", x = max_val, y = Inf, label = as.character(max_val),
                       vjust = 1.5, size = 3.4, fontface = "bold",
                       fill = "white", color = "black",
                       label.padding = unit(0.25, "lines"))
          p
        }
      })

      # --- Counts display ---
      output[[paste0("scounts_", metric)]] <- renderUI({
        total <- original_sample_count()
        pct_fmt <- function(n) sprintf("%.1f%%", 100 * n / total)

        if (is_ad) {
          count_blocks <- lapply(seq_along(GT_CLASSES), function(g) {
            gt_id <- GT_IDS[g]; gt_class <- GT_CLASSES[g]
            gt_mask <- match_gt(df$GT, gt_class)
            vals <- df[[metric]][gt_mask & !is.na(df[[metric]])]
            gt_total <- sum(gt_mask, na.rm = TRUE)
            min_val <- suppressWarnings(as.numeric(input[[paste0("smin_", metric, "_", gt_id)]]))
            max_val <- suppressWarnings(as.numeric(input[[paste0("smax_", metric, "_", gt_id)]]))
            has_min <- !is.null(min_val) && !is.na(min_val)
            has_max <- !is.null(max_val) && !is.na(max_val)
            if (!has_min && !has_max) return(NULL)
            if (gt_total == 0) return(NULL)
            gt_pct <- function(n) sprintf("%.1f%%", 100 * n / gt_total)
            if (has_min && has_max) {
              outside_n <- sum(vals < min_val | vals > max_val)
              within_n <- gt_total - outside_n
            } else if (has_min) {
              outside_n <- sum(vals < min_val)
              within_n <- gt_total - outside_n
            } else {
              outside_n <- sum(vals > max_val)
              within_n <- gt_total - outside_n
            }
            tags$div(style = "font-size:11px; margin-bottom:4px;",
              tags$b(paste0("GT=", gt_class, ": ")),
              tags$span(style = "color:#27ae60;",
                format(within_n, big.mark = ","), " kept (",
                gt_pct(within_n), ")"),
              " | ",
              tags$span(style = "color:#c0392b;",
                format(outside_n, big.mark = ","), " set to missing (",
                gt_pct(outside_n), ")")
            )
          })
          count_blocks <- Filter(Negate(is.null), count_blocks)
          if (length(count_blocks) == 0) return(NULL)
          tags$div(style = "background:#f7f7f7; border:1px solid #ddd; border-radius:4px; padding:6px 10px;",
            count_blocks
          )
        } else {
          vals <- df[[metric]]; vals <- vals[!is.na(vals)]
          min_val <- suppressWarnings(as.numeric(input[[paste0("smin_", metric)]]))
          max_val <- suppressWarnings(as.numeric(input[[paste0("smax_", metric)]]))
          has_min <- !is.null(min_val) && !is.na(min_val)
          has_max <- !is.null(max_val) && !is.na(max_val)
          if (!has_min && !has_max) return(NULL)
          below_n <- if (has_min) sum(vals < min_val) else NA
          above_n <- if (has_max) sum(vals > max_val) else NA
          if (has_min && has_max) {
            within_n  <- sum(vals >= min_val & vals <= max_val)
            removed_n <- total - within_n
            tags$div(style = "background:#f7f7f7; border:1px solid #ddd; border-radius:4px; padding:6px 10px; font-size:12px;",
              tags$b("Total genotypes: "), format(total, big.mark = ","), tags$br(),
              tags$span(style = "color:#c0392b;",
                tags$b("Below min: "), format(below_n, big.mark = ","),
                paste0(" (", pct_fmt(below_n), ")")), tags$br(),
              tags$span(style = "color:#c0392b;",
                tags$b("Above max: "), format(above_n, big.mark = ","),
                paste0(" (", pct_fmt(above_n), ")")), tags$br(),
              tags$span(style = "color:#27ae60;",
                tags$b("In range: "), format(within_n, big.mark = ","),
                paste0(" (", pct_fmt(within_n), ")")), tags$br(),
              tags$span(style = "color:#e67e22;",
                tags$b("Removed: "), format(removed_n, big.mark = ","),
                paste0(" (", pct_fmt(removed_n), ")"))
            )
          } else if (has_min) {
            kept_n <- sum(vals >= min_val)
            tags$div(style = "background:#f7f7f7; border:1px solid #ddd; border-radius:4px; padding:6px 10px; font-size:12px;",
              tags$b("Total genotypes: "), format(total, big.mark = ","), tags$br(),
              tags$span(style = "color:#c0392b;",
                tags$b("Below min: "), format(below_n, big.mark = ","),
                paste0(" (", pct_fmt(below_n), ")")), tags$br(),
              tags$span(style = "color:#27ae60;",
                tags$b("Kept (\u2265 min): "), format(kept_n, big.mark = ","),
                paste0(" (", pct_fmt(kept_n), ")"))
            )
          } else {
            kept_n <- sum(vals <= max_val)
            tags$div(style = "background:#f7f7f7; border:1px solid #ddd; border-radius:4px; padding:6px 10px; font-size:12px;",
              tags$b("Total genotypes: "), format(total, big.mark = ","), tags$br(),
              tags$span(style = "color:#c0392b;",
                tags$b("Above max: "), format(above_n, big.mark = ","),
                paste0(" (", pct_fmt(above_n), ")")), tags$br(),
              tags$span(style = "color:#27ae60;",
                tags$b("Kept (\u2264 max): "), format(kept_n, big.mark = ","),
                paste0(" (", pct_fmt(kept_n), ")"))
            )
          }
        }
      })
    })
  })

  # ---- Per-sample data preview ----
  output$sample_data_table <- renderDT({
    req(sample_data())
    show_cols <- c("CHROM", "POS", "SAMPLE", "GT",
                   "sample_DP", "SP", "ADREF", "ADALT", "ADREF_pct", "ADALT_pct",
                   "PL00", "PL01", "PL11")
    show_cols <- show_cols[show_cols %in% colnames(sample_data())]
    datatable(sample_data()[, show_cols, drop = FALSE],
              options = list(scrollX = TRUE, pageLength = 5),
              class = "compact cell-border stripe")
  })

  # ---- Per-sample BCFtools command ----
  sample_filter_info <- reactive({
    df <- tryCatch(sample_data(), error = function(e) NULL)

    std_filters <- list(); std_name_parts <- character(0)
    for (metric in SAMPLE_HIST_COLS) {
      min_val <- suppressWarnings(as.numeric(input[[paste0("smin_", metric)]]))
      max_val <- suppressWarnings(as.numeric(input[[paste0("smax_", metric)]]))
      has_min <- !is.null(min_val) && !is.na(min_val)
      has_max <- !is.null(max_val) && !is.na(max_val)
      if (!has_min && !has_max) next
      fmt_field <- if (metric == "sample_DP") "FORMAT/DP" else paste0("FORMAT/", metric)
      if (has_min) std_filters <- c(std_filters, paste0(fmt_field, " >= ", min_val))
      if (has_max) std_filters <- c(std_filters, paste0(fmt_field, " <= ", max_val))
      tag <- metric
      if (has_min && has_max) tag <- paste0(tag, "min", min_val, "max", max_val)
      else if (has_min)       tag <- paste0(tag, "min", min_val)
      else                    tag <- paste0(tag, "max", max_val)
      std_name_parts <- c(std_name_parts, tag)
    }

    dp_min <- suppressWarnings(as.numeric(input[["smin_sample_DP"]]))
    dp_max <- suppressWarnings(as.numeric(input[["smax_sample_DP"]]))
    has_dp_filter <- (!is.null(dp_min) && !is.na(dp_min)) || (!is.null(dp_max) && !is.na(dp_max))

    ad_setgt_cmds <- list(); ad_name_parts <- character(0)
    for (metric in SAMPLE_AD_COLS) {
      ad_idx <- if (metric == "ADREF_pct") 0 else 1
      for (g in seq_along(GT_CLASSES)) {
        gt_id <- GT_IDS[g]; gt_class <- GT_CLASSES[g]
        gt_bcf <- GT_BCF[gt_id]
        min_val <- suppressWarnings(as.numeric(input[[paste0("smin_", metric, "_", gt_id)]]))
        max_val <- suppressWarnings(as.numeric(input[[paste0("smax_", metric, "_", gt_id)]]))
        has_min <- !is.null(min_val) && !is.na(min_val)
        has_max <- !is.null(max_val) && !is.na(max_val)
        if (!has_min && !has_max) next

        conditions <- c()
        if (has_min) conditions <- c(conditions,
          paste0("FMT/AD[*:", ad_idx, "] < ", min_val,
                 " * (FMT/AD[*:0] + FMT/AD[*:1])"))
        if (has_max) conditions <- c(conditions,
          paste0("FMT/AD[*:", ad_idx, "] > ", max_val,
                 " * (FMT/AD[*:0] + FMT/AD[*:1])"))

        expr <- paste0('GT="', gt_bcf, '" && (', paste(conditions, collapse = " || "), ")")
        ad_setgt_cmds <- c(ad_setgt_cmds, list(list(expr = expr)))

        tag <- paste0(metric, "_", gt_class)
        if (has_min && has_max) tag <- paste0(tag, "min", min_val, "max", max_val)
        else if (has_min)       tag <- paste0(tag, "min", min_val)
        else                    tag <- paste0(tag, "max", max_val)
        ad_name_parts <- c(ad_name_parts, tag)
      }
    }

    pct_missing <- suppressWarnings(as.numeric(input$pct_missing))
    has_missing <- !is.null(pct_missing) && !is.na(pct_missing) && pct_missing > 0
    missing_name <- character(0)
    if (has_missing) missing_name <- paste0("Fmiss", pct_missing)

    list(
      std_filters = std_filters, std_name_parts = std_name_parts,
      ad_setgt_cmds = ad_setgt_cmds, ad_name_parts = ad_name_parts,
      has_dp_filter = has_dp_filter,
      has_missing = has_missing,
      pct_missing = if (has_missing) pct_missing else NULL,
      missing_name = missing_name
    )
  })

  output$sample_bcf_command <- renderText({
    info <- sample_filter_info()
    vcf  <- input$sample_vcf_name
    if (vcf == "") return("Set thresholds and VCF name above.")

    has_std <- length(info$std_filters) > 0
    has_ad  <- length(info$ad_setgt_cmds) > 0
    has_dp  <- info$has_dp_filter
    has_miss <- info$has_missing

    if (!has_std && !has_ad && !has_miss) return("Set thresholds above.")

    base <- sub("\\.(vcf\\.gz|vcf|bcf)$", "", vcf)
    ext  <- sub(paste0("^", gsub("\\.", "\\\\.", base)), "", vcf)
    flag <- if (grepl("\\.vcf\\.gz$", vcf)) "-Oz" else if (grepl("\\.bcf$", vcf)) "-Ob" else "-Ov"
    all_name_parts <- c(info$std_name_parts, info$ad_name_parts, info$missing_name)
    out <- paste0(base, "_", paste(all_name_parts, collapse = "_"), ext)

    stages <- c()
    current_input <- vcf

    if (has_std) {
      expr <- paste(info$std_filters, collapse = " && ")
      if (has_dp) {
        std_cmd <- paste0("bcftools filter \\\n  --set-GTs . \\\n  -i '", expr, "'")
      } else {
        std_cmd <- paste0("bcftools filter \\\n  -i '", expr, "'")
      }
      if (length(stages) == 0) {
        std_cmd <- paste0(std_cmd, " \\\n  ", current_input)
      }
      stages <- c(stages, std_cmd)
    }

    if (has_ad) {
      all_exprs <- sapply(info$ad_setgt_cmds, function(x) x$expr)
      combined_expr <- paste0("(", all_exprs, ")")
      combined_expr <- paste(combined_expr, collapse = " || ")
      ad_cmd <- paste0("bcftools +setGT \\\n  -- -t q -n . \\\n  -i '", combined_expr, "'")
      if (length(stages) == 0) {
        ad_cmd <- paste0("bcftools +setGT ", current_input, " \\\n  -- -t q -n . \\\n  -i '", combined_expr, "'")
      }
      stages <- c(stages, ad_cmd)
    }

    if (has_miss) {
      frac <- info$pct_missing / 100
      miss_cmd <- paste0("bcftools view \\\n  -i 'F_MISSING < ", sprintf("%.2f", frac), "'")
      if (length(stages) == 0) {
        miss_cmd <- paste0(miss_cmd, " \\\n  ", current_input)
      }
      stages <- c(stages, miss_cmd)
    }

    if (length(stages) == 1) {
      pipeline <- paste0(stages[[1]], " \\\n  ", flag, " -o ", out)
    } else {
      pipeline <- paste(stages, collapse = " | \\\n")
      pipeline <- paste0(pipeline, " \\\n  ", flag, " -o ", out)
    }

    paste0(pipeline,
           " && \\\nbcftools index ", out,
           " && \\\necho \"Filtered SNP count:\" && \\\nbcftools view ", out,
           " | grep -v '^#' | wc -l")
  })

  # ---- Per-sample filter summary ----
  sample_summary_df <- reactive({
    df <- tryCatch(sample_data(), error = function(e) NULL)
    rows <- list()

    if (!is.null(df)) {
      total <- nrow(df)

      for (metric in SAMPLE_HIST_COLS) {
        if (!(metric %in% colnames(df))) next
        min_val <- suppressWarnings(as.numeric(input[[paste0("smin_", metric)]]))
        max_val <- suppressWarnings(as.numeric(input[[paste0("smax_", metric)]]))
        has_min <- !is.null(min_val) && !is.na(min_val)
        has_max <- !is.null(max_val) && !is.na(max_val)
        if (!has_min && !has_max) next
        vals <- df[[metric]]
        range_str <- paste0(
          if (has_min) paste0("min=", min_val) else "",
          if (has_min && has_max) ", " else "",
          if (has_max) paste0("max=", max_val) else ""
        )
        if (has_min && has_max) kept_n <- sum(vals >= min_val & vals <= max_val, na.rm = TRUE)
        else if (has_min)       kept_n <- sum(vals >= min_val, na.rm = TRUE)
        else                    kept_n <- sum(vals <= max_val, na.rm = TRUE)
        removed_n <- total - kept_n
        rows <- c(rows, list(data.frame(
          Metric = SAMPLE_LABELS[metric], GT_Class = "All", Range = range_str,
          Total = total, Kept = kept_n, Removed = removed_n,
          Pct_Kept = sprintf("%.1f%%", 100 * kept_n / total),
          Pct_Removed = sprintf("%.1f%%", 100 * removed_n / total),
          stringsAsFactors = FALSE
        )))
      }

      for (metric in SAMPLE_AD_COLS) {
        if (!(metric %in% colnames(df))) next
        for (g in seq_along(GT_CLASSES)) {
          gt_id <- GT_IDS[g]; gt_class <- GT_CLASSES[g]
          min_val <- suppressWarnings(as.numeric(input[[paste0("smin_", metric, "_", gt_id)]]))
          max_val <- suppressWarnings(as.numeric(input[[paste0("smax_", metric, "_", gt_id)]]))
          has_min <- !is.null(min_val) && !is.na(min_val)
          has_max <- !is.null(max_val) && !is.na(max_val)
          if (!has_min && !has_max) next
          gt_mask <- match_gt(df$GT, gt_class)
          vals <- df[[metric]][gt_mask]
          gt_total <- sum(gt_mask, na.rm = TRUE)
          range_str <- paste0(
            if (has_min) paste0("min=", min_val) else "",
            if (has_min && has_max) ", " else "",
            if (has_max) paste0("max=", max_val) else ""
          )
          if (has_min && has_max) kept_n <- sum(vals >= min_val & vals <= max_val, na.rm = TRUE)
          else if (has_min)       kept_n <- sum(vals >= min_val, na.rm = TRUE)
          else                    kept_n <- sum(vals <= max_val, na.rm = TRUE)
          removed_n <- gt_total - kept_n
          rows <- c(rows, list(data.frame(
            Metric = SAMPLE_LABELS[metric], GT_Class = gt_class, Range = range_str,
            Total = gt_total, Kept = kept_n, Removed = removed_n,
            Pct_Kept = if (gt_total > 0) sprintf("%.1f%%", 100 * kept_n / gt_total) else "-",
            Pct_Removed = if (gt_total > 0) sprintf("%.1f%%", 100 * removed_n / gt_total) else "-",
            stringsAsFactors = FALSE
          )))
        }
      }

      pct_missing <- suppressWarnings(as.numeric(input$pct_missing))
      if (!is.null(pct_missing) && !is.na(pct_missing) && pct_missing > 0) {
        rows <- c(rows, list(data.frame(
          Metric = "F_MISSING", GT_Class = "All",
          Range = paste0("< ", pct_missing, "%"),
          Total = NA_integer_, Kept = NA_integer_, Removed = NA_integer_,
          Pct_Kept = "-", Pct_Removed = "-",
          stringsAsFactors = FALSE
        )))
      }
    }
    if (length(rows) == 0) return(data.frame(Metric = character(0)))
    do.call(rbind, rows)
  })

  output$has_sample_summary <- reactive({
    df <- sample_summary_df()
    !is.null(df) && nrow(df) > 0
  })
  outputOptions(output, "has_sample_summary", suspendWhenHidden = FALSE)

  output$sample_filter_summary_table <- renderDT({
    req(nrow(sample_summary_df()) > 0)
    datatable(sample_summary_df(), options = list(dom = "t", paging = FALSE),
              rownames = FALSE, class = "compact cell-border stripe")
  })

  output$download_sample_summary <- downloadHandler(
    filename = function() paste0("sample_filter_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv"),
    content  = function(file) write.table(sample_summary_df(), file, sep = "\t", row.names = FALSE, quote = FALSE)
  )

  # ---- Per-sample PNG export (2 columns) ----
  output$download_sample_png <- downloadHandler(
    filename = function() paste0("sample_metrics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content  = function(file) {
      df <- sample_data()
      available <- SAMPLE_FILTER_COLS[SAMPLE_FILTER_COLS %in% colnames(df)]
      plots <- lapply(seq_along(available), function(i) {
        metric  <- available[i]
        col_idx <- which(SAMPLE_ALL_COLS == metric)
        fill_col <- SAMPLE_RAINBOW[col_idx]
        is_ad   <- metric %in% SAMPLE_AD_COLS
        plot_df <- df[!is.na(df[[metric]]), , drop = FALSE]

        if (is_ad) {
          plot_df$GT <- as.factor(plot_df$GT)
          p <- ggplot(plot_df, aes(x = GT, y = .data[[metric]])) +
            geom_boxplot(fill = fill_col, color = "grey25", alpha = 0.85,
                         outlier.size = 0.5, outlier.alpha = 0.3) +
            scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                               limits = c(0, 1)) +
            labs(title = paste(SAMPLE_LABELS[metric], "by GT"),
                 x = "GT", y = SAMPLE_LABELS[metric]) +
            theme_minimal(base_size = 11) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
          gt_levels <- levels(plot_df$GT)
          for (g in seq_along(GT_CLASSES)) {
            gt_id <- GT_IDS[g]; gt_class <- GT_CLASSES[g]
            gt_match <- which(gt_levels %in% c(gt_class, gsub("/", "|", gt_class, fixed = TRUE)))
            if (length(gt_match) == 0) next
            x_pos <- gt_match[1]
            min_val <- suppressWarnings(as.numeric(input[[paste0("smin_", metric, "_", gt_id)]]))
            if (!is.null(min_val) && !is.na(min_val))
              p <- p + annotate("segment", x = x_pos - 0.4, xend = x_pos + 0.4,
                                y = min_val, yend = min_val,
                                linetype = "dashed", color = "red", linewidth = 0.9)
            max_val <- suppressWarnings(as.numeric(input[[paste0("smax_", metric, "_", gt_id)]]))
            if (!is.null(max_val) && !is.na(max_val))
              p <- p + annotate("segment", x = x_pos - 0.4, xend = x_pos + 0.4,
                                y = max_val, yend = max_val,
                                linetype = "dashed", color = "red", linewidth = 0.9)
          }
          p
        } else {
          p <- ggplot(plot_df, aes(x = .data[[metric]])) +
            geom_histogram(bins = 50, fill = fill_col, color = "grey25", alpha = 0.85) +
            labs(title = SAMPLE_LABELS[metric], x = SAMPLE_LABELS[metric], y = "Count") +
            theme_minimal(base_size = 11) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
          min_val <- suppressWarnings(as.numeric(input[[paste0("smin_", metric)]]))
          if (!is.null(min_val) && !is.na(min_val))
            p <- p + geom_vline(xintercept = min_val, linetype = "dashed",
                                color = "black", linewidth = 0.9) +
              annotate("label", x = min_val, y = Inf, label = as.character(min_val),
                       vjust = 1.5, size = 3.4, fontface = "bold", fill = "white",
                       label.padding = unit(0.25, "lines"))
          max_val <- suppressWarnings(as.numeric(input[[paste0("smax_", metric)]]))
          if (!is.null(max_val) && !is.na(max_val))
            p <- p + geom_vline(xintercept = max_val, linetype = "dashed",
                                color = "black", linewidth = 0.9) +
              annotate("label", x = max_val, y = Inf, label = as.character(max_val),
                       vjust = 1.5, size = 3.4, fontface = "bold", fill = "white",
                       label.padding = unit(0.25, "lines"))
          p
        }
      })
      n_plots <- length(plots)
      n_rows  <- ceiling(n_plots / 2)
      png(file, width = 1200, height = 600 * n_rows, res = 150)
      do.call(grid.arrange, c(plots, ncol = 2))
      dev.off()
    }
  )

  # ==========================================================================
  # TAB 3 - PL Viewer (view-only)
  # ==========================================================================
  observe({
    df <- tryCatch(sample_data(), error = function(e) NULL)
    if (is.null(df)) return()

    for (pl_col in SAMPLE_PL_COLS) {
      local({
        my_col <- pl_col
        col_idx <- which(SAMPLE_ALL_COLS == my_col)
        fill_col <- SAMPLE_RAINBOW[col_idx]
        output[[paste0("pl_plot_", my_col)]] <- renderPlot({
          req(my_col %in% colnames(df))
          plot_df <- df[!is.na(df[[my_col]]), , drop = FALSE]
          plot_df$GT <- as.factor(plot_df$GT)
          ggplot(plot_df, aes(x = GT, y = .data[[my_col]])) +
            geom_boxplot(fill = fill_col, color = "grey25", alpha = 0.85,
                         outlier.size = 0.5, outlier.alpha = 0.3) +
            labs(title = paste(my_col, "by Genotype"), x = "Genotype (GT)", y = my_col) +
            theme_minimal(base_size = 13) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        })
      })
    }
  })

  # ---- PL Viewer PNG export ----
  output$download_pl_png <- downloadHandler(
    filename = function() paste0("pl_plots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content  = function(file) {
      df <- sample_data()
      pl_available <- SAMPLE_PL_COLS[SAMPLE_PL_COLS %in% colnames(df)]
      plots <- lapply(pl_available, function(my_col) {
        col_idx  <- which(SAMPLE_ALL_COLS == my_col)
        fill_col <- SAMPLE_RAINBOW[col_idx]
        plot_df  <- df[!is.na(df[[my_col]]), , drop = FALSE]
        plot_df$GT <- as.factor(plot_df$GT)
        ggplot(plot_df, aes(x = GT, y = .data[[my_col]])) +
          geom_boxplot(fill = fill_col, color = "grey25", alpha = 0.85,
                       outlier.size = 0.5, outlier.alpha = 0.3) +
          labs(title = paste(my_col, "by Genotype"), x = "Genotype (GT)", y = my_col) +
          theme_minimal(base_size = 13) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      })
      png(file, width = 1800, height = 500, res = 150)
      do.call(grid.arrange, c(plots, ncol = length(plots)))
      dev.off()
    }
  )
}

shinyApp(ui, server)
