## Purity / Ploidy Fit Explorer
##
## Interactive Shiny app for assessing purity/ploidy CN fits.
##
## Given a balanced JaBbA gGraph and tumor coverage, computes per-segment
## coverage statistics (via Skilift::get_segstats) then lets the user drag
## intercept and slope sliders to see how integer-CN lines align with the
## segment-width-weighted coverage histogram, while back-calculating
## purity and ploidy.
##
## Run with:
##   Skilift::ppfit_explorer("/path/to/jabba.rds", "/path/to/cov.rds")
## or directly:
##   shiny::runApp(system.file("shiny/ppfit_explorer", package = "Skilift"))

suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(data.table)
  library(Skilift)
})

## Paths optionally pre-populated by Skilift::ppfit_explorer() via env vars
.init_jabba_path    <- Sys.getenv("SKILIFT_JABBA_PATH",     unset = "")
.init_coverage_path <- Sys.getenv("SKILIFT_COVERAGE_PATH",  unset = "")
.init_cov_field     <- Sys.getenv("SKILIFT_COVERAGE_FIELD", unset = "foreground")

# ── Helpers ───────────────────────────────────────────────────────────────────

## Compute coverage positions for integer CN states given purity/ploidy.
## Model: E[cov | CN=k] = intercept + slope * k
##   intercept = 2*(1-p) / (p*pl + 2*(1-p))   (coverage from normal cells)
##   slope     = p       / (p*pl + 2*(1-p))   (coverage gain per copy)
## This is the standard rel2abs transformation from skitools.
cn_cov_lines <- function(purity, ploidy, max_cn = 12) {
  denom     <- purity * ploidy + 2 * (1 - purity)
  intercept <- 2 * (1 - purity) / denom
  slope     <- purity / denom
  data.table(
    cn       = 0L:as.integer(max_cn),
    coverage = intercept + slope * (0:max_cn)
  )
}

## Convert normalized coverage to estimated CN
cov_to_cn <- function(coverage, purity, ploidy) {
  denom     <- purity * ploidy + 2 * (1 - purity)
  intercept <- 2 * (1 - purity) / denom
  slope     <- purity / denom
  (coverage - intercept) / slope
}

## Convert intercept / slope to purity / ploidy under the same rel2abs model.
## purity = 2*slope / (intercept + 2*slope)
## ploidy = (1 - intercept) / slope
infer_pp_from_line <- function(intercept, slope) {
  purity <- ifelse(
    slope <= 0,
    0,
    2 * slope / (intercept + 2 * slope)
  )

  ploidy <- ifelse(
    slope <= 0,
    NA_real_,
    (1 - intercept) / slope
  )

  list(
    purity = as.numeric(purity),
    ploidy = as.numeric(ploidy)
  )
}

## Extract purity / ploidy stored in a gGraph object
extract_pp <- function(jab) {
  p  <- if (!is.null(jab$purity))      jab$purity
        else if (!is.null(jab$meta$purity))  jab$meta$purity
        else NA_real_
  pl <- if (!is.null(jab$ploidy))      jab$ploidy
        else if (!is.null(jab$meta$ploidy))  jab$meta$ploidy
        else NA_real_
  list(purity = as.numeric(p), ploidy = as.numeric(pl))
}

AUTOSOME_NAMES <- c(as.character(1:22), paste0("chr", 1:22))

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- fluidPage(

      tags$head(tags$style(HTML("
        .sidebar-section { margin-bottom: 14px; }
        .ref-box  { background:#eaf4fb; border-left:3px solid #3498db;
                    padding:6px 10px; font-size:0.85em; border-radius:3px; }
        .ok-box   { background:#eafaf1; border-left:3px solid #27ae60;
                    padding:6px 10px; font-size:0.85em; border-radius:3px; }
        .warn-box { background:#fef9e7; border-left:3px solid #f39c12;
                    padding:6px 10px; font-size:0.85em; border-radius:3px; }
        .err-box  { background:#fdedec; border-left:3px solid #e74c3c;
                    padding:6px 10px; font-size:0.85em; border-radius:3px; }
        h5 { font-weight:600; color:#2c3e50; margin-top:10px; }
        .summary-bar { font-family:monospace; font-size:0.82em;
                       background:#f4f6f7; padding:6px 10px; border-radius:3px; }
    "))),

  titlePanel("Purity / Ploidy Fit Explorer"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      # ── Files ──────────────────────────────────────────────────────
      div(class = "sidebar-section",
          h5("Input Files"),
          textInput("jabba_path",     "Balanced JaBbA gGraph (.rds):",
                    value = .init_jabba_path,
                    placeholder = "/path/to/balanced.rds"),
          textInput("coverage_path",  "Tumor coverage (.rds):",
                    value = .init_coverage_path,
                    placeholder = "/path/to/cov.rds"),
          textInput("coverage_field", "Coverage field:",
                    value = .init_cov_field),
          actionButton("load_btn", "Load & Compute Segstats",
                       class = "btn-primary", width = "100%")
          ),

      # ── Status ─────────────────────────────────────────────────────
      uiOutput("status_ui"),
      hr(),

      # ── Intercept / slope sliders ─────────────────────────────────
      div(class = "sidebar-section",
          h5("Intercept / Slope"),
          uiOutput("jabba_ref_ui"),
          sliderInput("intercept", "Intercept:",
                      min = 0.00, max = 1.00, value = 0.50, step = 0.005),
          sliderInput("slope", "Slope:",
                      min = 0.001, max = 1.00, value = 0.25, step = 0.001),
          div(class = "ok-box", uiOutput("derived_pp_ui")),
          actionButton("reset_btn", "Reset to JaBbA values",
                       class = "btn-sm btn-default", width = "100%")
          ),
      hr(),

      # ── Display options ────────────────────────────────────────────
      div(class = "sidebar-section",
          h5("Display Options"),
          numericInput("max_cn",  "Max CN:",           value = 10, min = 1, max = 40),
          numericInput("n_bins",  "Histogram bins:",   value = 150, min = 20, max = 500),
          checkboxInput("autosomes_only",  "Autosomes only",         value = TRUE),
          checkboxInput("show_jabba_ref", "Show JaBbA reference lines", value = TRUE),
          radioButtons("view_mode", "X-axis mode:",
                       choices  = c("Normalized coverage" = "coverage",
                                    "Estimated CN"        = "cn"),
                       selected = "coverage"
                       )
          )
    ),

    mainPanel(
      width = 9,
      plotOutput("main_plot", height = "580px"),
      br(),
      uiOutput("summary_ui")
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  rv <- reactiveValues(
    segstats    = NULL,
    jabba_p     = NULL,
    jabba_pl    = NULL,
    jabba_i     = NULL,
    jabba_s     = NULL,
    status_msg  = "No data loaded.",
    status_type = "warn"   # "warn" | "ok" | "err" | "info"
  )

  current_pp <- reactive({
    infer_pp_from_line(input$intercept, input$slope)
  })

  # ── Load data ─────────────────────────────────────────────────────────────
  ## Core load logic extracted so it can be triggered by the button OR
  ## automatically on startup when paths are pre-supplied.
  do_load <- function(jabba_path, coverage_path, cov_field) {
    if (!nzchar(jabba_path) || !nzchar(coverage_path)) {
      rv$status_msg  <- "Please fill in both file paths."
      rv$status_type <- "warn"
      return()
    }

    rv$status_msg  <- "Loading… computing segstats (may take a minute)"
    rv$status_type <- "info"

    withCallingHandlers(
      tryCatch({
        if (!file.exists(jabba_path))
          stop("JaBbA file not found: ", jabba_path)
        if (!file.exists(coverage_path))
          stop("Coverage file not found: ", coverage_path)

        jab <- Skilift:::process_jabba(jabba_path)
        pp  <- extract_pp(jab)

        colnames_check <- c("mean", "sd", "var", "nbins_ok",
                            "nbins_nafrac", "raw_mean", "raw_var")
        has_segstats <- all(colnames_check %in% names(jab$nodes$dt))

        if (has_segstats) {
          segstats_dt <- jab$nodes$dt
          mu   <- segstats_dt$mean
          w    <- as.numeric(segstats_dt$width)
          mu[is.infinite(mu)] <- NA
          w[is.na(mu)] <- NA
          sw   <- sum(w, na.rm = TRUE)
          mutl <- sum(mu * w, na.rm = TRUE)
          segstats_dt[, mean := mu * (sw / mutl)]
        } else {

          segstats_dt <- Skilift::get_segstats(
                                    balanced_jabba_gg = jabba_path,
                                    tumor_coverage    = coverage_path,
                                    coverage_field    = cov_field
                                  )
        }

        rv$jabba_p  <- pp$purity
        rv$jabba_pl <- pp$ploidy
        rv$segstats <- segstats_dt

        if (!is.na(pp$purity) && !is.na(pp$ploidy)) {
          denom      <- pp$purity * pp$ploidy + 2 * (1 - pp$purity)
          rv$jabba_i <- 2 * (1 - pp$purity) / denom
          rv$jabba_s <- pp$purity / denom
        } else {
          rv$jabba_i <- NULL
          rv$jabba_s <- NULL
        }

        rv$status_msg <- sprintf(
          "Loaded %d segments.%s",
          nrow(segstats_dt),
          if (!is.na(pp$purity))
            sprintf("  JaBbA: purity=%.3f  ploidy=%.3f",
                    pp$purity, pp$ploidy)
          else ""
        )
        rv$status_type <- "ok"

        if (!is.null(rv$jabba_i) && !is.na(rv$jabba_i))
          updateSliderInput(session, "intercept", value = round(rv$jabba_i, 3))
        if (!is.null(rv$jabba_s) && !is.na(rv$jabba_s))
          updateSliderInput(session, "slope", value = round(rv$jabba_s, 3))

      }, error = function(e) {
        rv$status_msg  <- paste("Error:", conditionMessage(e))
        rv$status_type <- "err"
      }),
      message = function(m) invokeRestart("muffleMessage")
    )
  }

  ## Auto-load if paths were pre-supplied via ppfit_explorer()
  session$onFlushed(function() {
    if (nzchar(.init_jabba_path) && nzchar(.init_coverage_path))
      do_load(.init_jabba_path, .init_coverage_path, .init_cov_field)
  }, once = TRUE)

  observeEvent(input$load_btn, {
    do_load(
      jabba_path    = trimws(input$jabba_path),
      coverage_path = trimws(input$coverage_path),
      cov_field     = trimws(input$coverage_field)
    )
  })

  # ── Reset sliders to JaBbA-derived intercept / slope ─────────────────────
  observeEvent(input$reset_btn, {
    if (!is.null(rv$jabba_i) && !is.na(rv$jabba_i))
      updateSliderInput(session, "intercept", value = round(rv$jabba_i, 3))
    if (!is.null(rv$jabba_s) && !is.na(rv$jabba_s))
      updateSliderInput(session, "slope", value = round(rv$jabba_s, 3))
  })

  # ── Status UI ─────────────────────────────────────────────────────────────
  output$status_ui <- renderUI({
    cls <- switch(rv$status_type,
                  ok   = "ok-box",
                  err  = "err-box",
                  info = "ref-box",
                  "warn-box"
                  )
    div(class = cls, rv$status_msg)
  })

  # ── JaBbA reference badge ─────────────────────────────────────────────────
  output$jabba_ref_ui <- renderUI({
    p  <- rv$jabba_p
    pl <- rv$jabba_pl
    ji <- rv$jabba_i
    js <- rv$jabba_s
    if (!is.null(p) && !is.na(p)) {
      div(class = "ref-box",
          sprintf(
            "JaBbA fit — purity: %.3f   ploidy: %.3f   intercept: %.4f   slope: %.4f",
            p, pl, ji, js
          ))
    }
  })

  output$derived_pp_ui <- renderUI({
    pp <- current_pp()
    HTML(sprintf(
      "<b>Derived purity:</b> %.4f<br/><b>Derived ploidy:</b> %s",
      pp$purity,
      if (is.na(pp$ploidy)) "NA" else sprintf("%.4f", pp$ploidy)
    ))
  })

  # ── Filtered segstats ─────────────────────────────────────────────────────
  filtered_segs <- reactive({
    dt <- rv$segstats
    req(!is.null(dt))
    dt <- dt[!is.na(mean) & is.finite(mean) & mean > 0]
    if (isTRUE(input$autosomes_only))
      dt <- dt[seqnames %in% AUTOSOME_NAMES]
    dt
  })

  # ── Main plot ─────────────────────────────────────────────────────────────
  output$main_plot <- renderPlot({
    dt  <- filtered_segs()
    i   <- input$intercept
    s   <- input$slope
    pp  <- current_pp()
    p   <- pp$purity
    pl  <- pp$ploidy
    mcn <- input$max_cn

    validate(
      need(!is.na(pl) && is.finite(pl), "Current intercept/slope imply undefined ploidy."),
      need(p > 0, "Current intercept/slope imply zero purity; CN transform is undefined."),
      need(pl > 0, "Current intercept/slope imply non-positive ploidy."),
      need(i >= 0, "Intercept must be non-negative."),
      need(s > 0, "Slope must be > 0.")
    )

    # Segment weight column
    wt_col <- if ("nbins_ok" %in% names(dt)) "nbins_ok" else "width"
    wt_label <- if (wt_col == "nbins_ok") "bins (nbins_ok)" else "bp (width)"

    # Current fit lines
    lines <- data.table(
      cn = 0L:as.integer(mcn),
      coverage = i + s * (0:mcn)
    )

    # JaBbA reference lines (only when stored pp exists)
    jab_p  <- rv$jabba_p
    jab_pl <- rv$jabba_pl
    jab_i  <- rv$jabba_i
    jab_s  <- rv$jabba_s
    has_jabba_ref <- isTRUE(input$show_jabba_ref) &&
      !is.null(jab_i) && !is.na(jab_i) &&
      !(abs(jab_i - i) < 1e-6 && abs(jab_s - s) < 1e-6)
    jlines <- if (has_jabba_ref) {
                data.table(
                  cn = 0L:as.integer(mcn),
                  coverage = jab_i + jab_s * (0:mcn)
                )
              } else NULL

    if (input$view_mode == "coverage") {
      # ── Coverage-space plot ──────────────────────────────────────
      xlim_lo <- max(0, min(lines$coverage) - 0.15)
      xlim_hi <- max(lines$coverage) * 1.12

      g <- ggplot(dt, aes(x = mean, weight = .data[[wt_col]])) +
        geom_histogram(bins = input$n_bins,
                       fill = "#4e79a7", color = NA, alpha = 0.85)

      if (!is.null(jlines)) {
        g <- g +
          geom_vline(data = jlines,
                     aes(xintercept = coverage),
                     color = "grey60", linetype = "dotted",
                     linewidth = 0.5, alpha = 0.7,
                     inherit.aes = FALSE)
      }

      g <- g +
        geom_vline(data = lines,
                   aes(xintercept = coverage),
                   color = "#e15759", linetype = "dashed",
                   linewidth = 0.65,
                   inherit.aes = FALSE) +
        geom_text(data = lines,
                  aes(x = coverage, label = cn),
                  y = Inf, vjust = 1.4, hjust = -0.25,
                  size = 3.2, color = "#e15759", fontface = "bold",
                  inherit.aes = FALSE) +
        coord_cartesian(xlim = c(xlim_lo, xlim_hi)) +
        labs(
          title = sprintf(
            "Segment coverage distribution — intercept: %.4f   slope: %.4f",
            i, s),
          subtitle = paste0(
            sprintf("Derived purity: %.4f   ploidy: %.4f", p, pl),
            if (has_jabba_ref)
              sprintf("   |   Gray dotted = JaBbA fit (purity: %.3f, ploidy: %.3f)",
                      jab_p, jab_pl)
            else ""
          ),
          x = "Normalized segment mean coverage",
          y = paste("Total", wt_label, "(segment-width weighted)")
        ) +
        theme_bw(base_size = 14) +
        theme(
          plot.title    = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 10, color = "grey40")
        )

    } else {
      # ── CN-space plot ────────────────────────────────────────────
      dt_plot <- copy(dt)
      dt_plot[, est_cn := (mean - i) / s]
      dt_plot <- dt_plot[is.finite(est_cn) &
                         est_cn >= -0.5 & est_cn <= (mcn + 0.5)]

      g <- ggplot(dt_plot, aes(x = est_cn, weight = .data[[wt_col]])) +
        geom_histogram(bins = input$n_bins,
                       fill = "#4e79a7", color = NA, alpha = 0.85) +
        geom_vline(xintercept = 0:mcn,
                   color = "#e15759", linetype = "dashed",
                   linewidth = 0.65) +
        geom_text(data = data.table(cn = 0:mcn),
                  aes(x = cn, label = cn),
                  y = Inf, vjust = 1.4, hjust = -0.25,
                  size = 3.2, color = "#e15759", fontface = "bold",
                  inherit.aes = FALSE) +
        scale_x_continuous(breaks = 0:mcn, minor_breaks = NULL) +
        coord_cartesian(xlim = c(-0.5, mcn + 0.5)) +
        labs(
          title = sprintf(
            "Estimated CN distribution — intercept: %.4f   slope: %.4f",
            i, s),
          subtitle = sprintf(
            "Derived purity: %.4f   ploidy: %.4f   |   Peaks should align with integer CN lines for correct fit",
            p, pl
          ),
          x = "Estimated copy number",
          y = paste("Total", wt_label, "(segment-width weighted)")
        ) +
        theme_bw(base_size = 14) +
        theme(
          plot.title    = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 10, color = "grey40")
        )
    }

    g
  })

  # ── Summary bar ───────────────────────────────────────────────────────────
  output$summary_ui <- renderUI({
    dt <- filtered_segs()
    if (is.null(dt)) return(NULL)

    i  <- input$intercept
    s  <- input$slope
    pp <- current_pp()

    n_segs    <- nrow(dt)
    n_chroms  <- length(unique(dt$seqnames))
    wt_col    <- if ("nbins_ok" %in% names(dt)) "nbins_ok" else "width"
    wmean_cov <- weighted.mean(dt$mean, dt[[wt_col]], na.rm = TRUE)

    msg <- sprintf(
      "%d segments │ %d chromosomes │ weighted mean cov: %.4f  │  intercept (CN=0 cov): %.4f  │  slope: %.4f  │  purity: %.4f  │  ploidy: %s",
      n_segs, n_chroms, wmean_cov, i, s, pp$purity,
      if (is.na(pp$ploidy)) "NA" else sprintf("%.4f", pp$ploidy)
    )
    div(class = "summary-bar", msg)
  })
}

shinyApp(ui, server)
