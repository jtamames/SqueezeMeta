library(shiny)
library(shinyjs)
library(shinyFiles)
library(bslib)
library(SQMtools)
library(DT)

# ─────────────────────────────────────────────
#  TEMA CLARO
# ─────────────────────────────────────────────
sqm_theme <- bs_theme(
  version      = 5,
  bg           = "#f7f9fc",
  fg           = "#1a2a3a",
  primary      = "#1a6eb5",
  secondary    = "#e8eef5",
  success      = "#1a9e6e",
  info         = "#3b9ede",
  font_scale   = 0.92,
  base_font    = font_google("IBM Plex Sans"),
  heading_font = font_google("IBM Plex Mono")
)

# ─────────────────────────────────────────────
#  HELPERS UI
# ─────────────────────────────────────────────
stat_card <- function(id, label, icon_char = "◈") {
  card(
    class = "stat-card",
    card_body(
      class = "p-3",
      tags$div(class = "stat-icon", icon_char),
      tags$div(class = "stat-value", textOutput(id, inline = TRUE)),
      tags$div(class = "stat-label", label)
    )
  )
}

# ─────────────────────────────────────────────
#  CSS
# ─────────────────────────────────────────────
custom_css <- "
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;500&display=swap');

:root {
  --blue:   #1a6eb5;
  --teal:   #1a9e6e;
  --bg:     #f7f9fc;
  --panel:  #eef2f7;
  --card:   #ffffff;
  --border: #d0dae6;
  --muted:  #7a90a8;
  --text:   #1a2a3a;
}

body { background: var(--bg); color: var(--text); }

/* Navbar */
.navbar {
  background: linear-gradient(90deg, #0e4a82 0%, #1a6eb5 100%) !important;
  border-bottom: 3px solid #1a9e6e !important;
  padding: 0.5rem 1.5rem !important;
}
.navbar-brand {
  font-family: 'IBM Plex Mono', monospace !important;
  font-weight: 600;
  color: #ffffff !important;
  letter-spacing: 0.05em;
  font-size: 1.1rem;
}
.nav-link { color: rgba(255,255,255,0.65) !important; font-size: 0.85rem; }
.nav-link:hover { color: #ffffff !important; }
.nav-link.active { color: #ffffff !important; border-bottom: 2px solid #1a9e6e; }

/* Cards */
.card {
  background: var(--card) !important;
  border: 1px solid var(--border) !important;
  border-radius: 8px !important;
  box-shadow: 0 1px 4px rgba(26,42,58,0.07) !important;
}
.card-header {
  background: var(--panel) !important;
  border-bottom: 1px solid var(--border) !important;
  font-family: 'IBM Plex Mono', monospace;
  font-size: 0.78rem;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: var(--blue) !important;
  padding: 0.6rem 1rem !important;
}

/* Stat cards */
.stat-card { border-left: 3px solid var(--blue) !important; }
.stat-icon { font-size: 1.4rem; color: var(--blue); margin-bottom: 4px; }
.stat-value {
  font-family: 'IBM Plex Mono', monospace;
  font-size: 1.6rem;
  font-weight: 600;
  line-height: 1;
  color: var(--text);
}
.stat-label {
  font-size: 0.72rem;
  color: var(--muted);
  text-transform: uppercase;
  letter-spacing: 0.06em;
  margin-top: 4px;
}

/* Sidebar */
.bslib-sidebar-layout > .sidebar {
  background: var(--panel) !important;
  border-right: 1px solid var(--border) !important;
}

/* Inputs */
.form-control, .form-select {
  background: #ffffff !important;
  border: 1px solid var(--border) !important;
  color: var(--text) !important;
  font-size: 0.85rem !important;
}
.form-control:focus, .form-select:focus {
  border-color: var(--blue) !important;
  box-shadow: 0 0 0 2px rgba(26,110,181,0.12) !important;
}
.form-label {
  font-size: 0.78rem;
  color: var(--muted);
  text-transform: uppercase;
  letter-spacing: 0.05em;
  margin-bottom: 3px;
}

/* Buttons */
.btn-primary {
  background: var(--blue) !important;
  border-color: var(--blue) !important;
  color: #ffffff !important;
  font-weight: 600;
  font-size: 0.82rem;
  letter-spacing: 0.04em;
}
.btn-primary:hover { background: #1558a0 !important; }
.btn-outline-secondary {
  border-color: var(--border) !important;
  color: var(--muted) !important;
  font-size: 0.82rem;
  background: #ffffff !important;
}
.btn-outline-secondary:hover {
  border-color: var(--blue) !important;
  color: var(--blue) !important;
}

/* Project badge */
.project-badge {
  display: inline-block;
  background: rgba(26,110,181,0.08);
  border: 1px solid rgba(26,110,181,0.25);
  color: var(--blue);
  border-radius: 4px;
  padding: 2px 8px;
  font-size: 0.75rem;
  font-family: 'IBM Plex Mono', monospace;
  margin: 2px;
}

/* Section divider */
.section-divider {
  border: none;
  border-top: 1px solid var(--border);
  margin: 10px 0;
}

/* DT table */
.dataTables_wrapper { color: var(--text) !important; }
table.dataTable { background: #ffffff !important; color: var(--text) !important; }
table.dataTable thead th {
  background: var(--panel) !important;
  color: var(--blue) !important;
  border-bottom: 1px solid var(--border) !important;
  font-family: 'IBM Plex Mono', monospace;
  font-size: 0.75rem;
  letter-spacing: 0.05em;
}
table.dataTable tbody tr:hover { background: #eef5fc !important; }
.dataTables_filter input, .dataTables_length select {
  background: #ffffff !important;
  border: 1px solid var(--border) !important;
  color: var(--text) !important;
}
.dataTables_info, .dataTables_paginate { color: var(--muted) !important; font-size: 0.8rem; }
.paginate_button.current {
  background: var(--blue) !important;
  color: #ffffff !important;
  border-radius: 4px;
}

/* shinyFiles button */
.btn-default {
  background: #ffffff !important;
  border: 1px solid var(--border) !important;
  color: var(--muted) !important;
  font-size: 0.82rem !important;
}
.btn-default:hover { border-color: var(--blue) !important; color: var(--blue) !important; }

/* Scrollbar */
::-webkit-scrollbar { width: 5px; height: 5px; }
::-webkit-scrollbar-track { background: var(--bg); }
::-webkit-scrollbar-thumb { background: var(--border); border-radius: 3px; }
::-webkit-scrollbar-thumb:hover { background: var(--blue); }
"

# ─────────────────────────────────────────────
#  UI
# ─────────────────────────────────────────────
ui <- page_navbar(
  title = "SQM Viewer",
  theme = sqm_theme,

  navbar_options = navbar_options(theme = "dark", bg = "#0e4a82"),

  header = tagList(
    useShinyjs(),
    tags$head(tags$style(HTML(custom_css)))
  ),

  # ══════════════════════════════════════════
  #  TAB 1 — PROYECTO
  # ══════════════════════════════════════════
  nav_panel(
    "Proyecto",

    layout_sidebar(
      sidebar = sidebar(
        width = 280, open = TRUE,

        tags$div(class = "form-label mt-1", "Directorio de tablas"),
        shinyDirButton(
          "dir_A", "Seleccionar directorio", "Elige el directorio de tablas",
          multiple = FALSE, class = "btn-default w-100 mb-1"
        ),
        tags$div(
          style = "font-size:0.72rem; color:var(--muted); word-break:break-all; margin-bottom:8px;",
          textOutput("path_A", inline = TRUE)
        ),

        tags$div(class = "form-label", "Tipo de carga"),
        radioButtons(
          "load_type_A", NULL,
          choices = c("Completo (loadSQM)" = "full", "Ligero (loadSQMlite)" = "lite"),
          inline  = TRUE
        ),

        actionButton("load_A", "Cargar proyecto", class = "btn-primary w-100 mb-2"),

        hr(class = "section-divider"),
        uiOutput("project_status_ui")
      ),

      fluidRow(
        column(3, stat_card("stat_contigs", "Contigs",  "⬡")),
        column(3, stat_card("stat_orfs",    "ORFs",     "◈")),
        column(3, stat_card("stat_samples", "Muestras", "◉")),
        column(3, stat_card("stat_taxa",    "Taxones",  "✦"))
      ),

      br(),

      fluidRow(
        column(12,
          card(
            card_header("Muestras"),
            card_body(uiOutput("samples_badges"))
          )
        )
      )
    )
  ),

  # ══════════════════════════════════════════
  #  TAB 2 — GRÁFICOS
  # ══════════════════════════════════════════
  nav_panel(
    "Gráficos",

    layout_sidebar(
      sidebar = sidebar(
        width = 270,

        tags$div(class = "form-label mt-1", "Tipo de gráfico"),
        selectInput("plot_type", NULL,
          choices = c(
            "Taxonomía (barras)"  = "taxonomy_bar",
            "Taxonomía (pie)"     = "taxonomy_pie",
            "Funciones COG"       = "func_cog",
            "Funciones KEGG"      = "func_kegg",
            "Funciones PFAM"      = "func_pfam",
            "Binning"             = "bins"
          )
        ),

        hr(class = "section-divider"),
        uiOutput("plot_controls_ui"),
        hr(class = "section-divider"),

        actionButton("do_plot", "Generar gráfico", class = "btn-primary w-100 mt-1"),
        br(), br(),
        downloadButton("download_plot", "Descargar PNG", class = "btn-outline-secondary w-100")
      ),

      card(
        card_header(
          div(
            style = "display:flex; justify-content:space-between; align-items:center;",
            span("Visualización"),
            uiOutput("plot_status_badge")
          )
        ),
        card_body(
          class = "p-2",
          plotOutput("sqm_plot", height = "560px")
        )
      )
    )
  ),

  # ══════════════════════════════════════════
  #  TAB 3 — TABLAS
  # ══════════════════════════════════════════
  nav_panel(
    "Tablas",

    layout_sidebar(
      sidebar = sidebar(
        width = 270,

        tags$div(class = "form-label mt-1", "Tabla"),
        selectInput("table_type", NULL,
          choices = c(
            "Contigs"             = "contigs",
            "ORFs"                = "orfs",
            "Taxonomía — phylum"  = "tax_phylum",
            "Taxonomía — genus"   = "tax_genus",
            "Funciones — COG"     = "fun_cog",
            "Funciones — KEGG"    = "fun_kegg"
          )
        ),

        hr(class = "section-divider"),
        tags$div(class = "form-label", "Filtrar muestras"),
        uiOutput("table_sample_filter"),

        hr(class = "section-divider"),
        downloadButton("download_table", "Descargar CSV", class = "btn-outline-secondary w-100")
      ),

      card(
        card_header("Datos"),
        card_body(class = "p-2", DTOutput("data_table"))
      )
    )
  )
)

# ─────────────────────────────────────────────
#  SERVER
# ─────────────────────────────────────────────
server <- function(input, output, session) {

  roots    <- c(home = normalizePath("~"), root = "/")
  sqm_A    <- reactiveVal(NULL)
  status_A <- reactiveVal("idle")

  # ── Dir chooser ───────────────────────────
  shinyDirChoose(input, "dir_A", roots = roots)

  path_A <- reactive({
    req(input$dir_A)
    parseDirPath(roots, input$dir_A)
  })

  output$path_A <- renderText({
    tryCatch(path_A(), error = function(e) "")
  })

  # ── Cargar proyecto ───────────────────────
  observeEvent(input$load_A, {
    req(path_A())
    status_A("loading")
    tryCatch({
      proj <- if (input$load_type_A == "full") loadSQM(path_A()) else loadSQMlite(path_A())
      sqm_A(proj)
      status_A("ready")
    }, error = function(e) {
      status_A("error")
      showNotification(paste("Error:", e$message), type = "error", duration = 8)
    })
  })

  # ── Status ────────────────────────────────
  output$project_status_ui <- renderUI({
    s   <- status_A()
    col <- switch(s, idle = "#7a90a8", loading = "#3b9ede", ready = "#1a9e6e", error = "#c0392b")
    ico <- switch(s, idle = "○", loading = "◌", ready = "●", error = "✕")
    tags$div(
      style = "font-size:0.8rem;",
      tags$span(style = paste0("color:", col, "; margin-right:5px;"), ico),
      tags$span(style = "color:#7a90a8;", "Estado: "),
      tags$span(style = paste0("color:", col, "; font-weight:600;"), toupper(s))
    )
  })

  # ── Stats ─────────────────────────────────
  output$stat_contigs <- renderText({
    req(sqm_A())
    tryCatch(format(nrow(sqm_A()$contigs$table), big.mark = ","), error = function(e) "—")
  })
  output$stat_orfs <- renderText({
    req(sqm_A())
    tryCatch(format(nrow(sqm_A()$orfs$table), big.mark = ","), error = function(e) "—")
  })
  output$stat_samples <- renderText({
    req(sqm_A())
    tryCatch(as.character(length(sqm_A()$samples)), error = function(e) "—")
  })
  output$stat_taxa <- renderText({
    req(sqm_A())
    tryCatch(format(nrow(sqm_A()$taxa$phylum$abund), big.mark = ","), error = function(e) "—")
  })

  output$samples_badges <- renderUI({
    req(sqm_A())
    samples <- tryCatch(sqm_A()$samples, error = function(e) NULL)
    req(samples)
    tagList(lapply(samples, function(s) tags$span(class = "project-badge", s)))
  })

  # ── Controles dinámicos de gráfico ────────
  output$plot_controls_ui <- renderUI({
    pt <- input$plot_type
    if (is.null(pt)) return(NULL)

    rank_choices <- c(
      "Phylum" = "phylum", "Class" = "class", "Order" = "order",
      "Family" = "family", "Genus" = "genus", "Species" = "species"
    )

    if (pt %in% c("taxonomy_bar", "taxonomy_pie")) {
      tagList(
        tags$div(class = "form-label", "Rango taxonómico"),
        selectInput("tax_rank", NULL, choices = rank_choices),
        tags$div(class = "form-label", "Nº taxones"),
        numericInput("n_taxa", NULL, value = 15, min = 5, max = 50, step = 5)
      )
    } else if (pt %in% c("func_cog", "func_kegg", "func_pfam")) {
      tagList(
        tags$div(class = "form-label", "Nº funciones"),
        numericInput("n_funcs", NULL, value = 20, min = 5, max = 60, step = 5)
      )
    } else {
      NULL
    }
  })

  # ── Generar gráfico ───────────────────────
  plot_reactive <- eventReactive(input$do_plot, {
    req(sqm_A())
    proj <- sqm_A()
    pt   <- input$plot_type

    if      (pt == "taxonomy_bar") plotTaxonomy(proj,  rank = input$tax_rank, N = input$n_taxa)
    else if (pt == "taxonomy_pie") plotTaxonomy(proj,  rank = input$tax_rank, N = input$n_taxa, type = "pie")
    else if (pt == "func_cog")     plotFunctions(proj, fun_level = "COG",  N = input$n_funcs)
    else if (pt == "func_kegg")    plotFunctions(proj, fun_level = "KEGG", N = input$n_funcs)
    else if (pt == "func_pfam")    plotFunctions(proj, fun_level = "PFAM", N = input$n_funcs)
    else if (pt == "bins")         plotBins(proj)
  })

  output$sqm_plot <- renderPlot({
    plot_reactive()
  }, bg = "#ffffff")

  output$plot_status_badge <- renderUI({
    if (is.null(sqm_A())) {
      tags$span(
        class = "badge",
        style = "background:#eef2f7; color:#7a90a8; font-size:0.72rem; border:1px solid #d0dae6;",
        "Sin proyecto"
      )
    } else {
      tags$span(
        class = "badge",
        style = "background:rgba(26,158,110,0.1); color:#1a9e6e; font-size:0.72rem; border:1px solid rgba(26,158,110,0.3);",
        "● Listo"
      )
    }
  })

  output$download_plot <- downloadHandler(
    filename = function() paste0("sqm_plot_", Sys.Date(), ".png"),
    content  = function(file) {
      png(file, width = 1400, height = 900, res = 150, bg = "#ffffff")
      print(plot_reactive())
      dev.off()
    }
  )

  # ── Filtro de muestras ────────────────────
  output$table_sample_filter <- renderUI({
    req(sqm_A())
    samples <- tryCatch(sqm_A()$samples, error = function(e) NULL)
    req(samples)
    checkboxGroupInput("selected_samples", NULL, choices = samples, selected = samples)
  })

  # ── Datos de tabla ────────────────────────
  get_table_data <- reactive({
    req(sqm_A())
    proj <- sqm_A()
    tt   <- input$table_type
    smp  <- input$selected_samples

    tryCatch({
      if      (tt == "contigs")    as.data.frame(proj$contigs$table)
      else if (tt == "orfs")       as.data.frame(proj$orfs$table)
      else if (tt == "tax_phylum") {
        d <- as.data.frame(proj$taxa$phylum$abund)
        if (!is.null(smp)) d[, colnames(d) %in% smp, drop = FALSE] else d
      }
      else if (tt == "tax_genus") {
        d <- as.data.frame(proj$taxa$genus$abund)
        if (!is.null(smp)) d[, colnames(d) %in% smp, drop = FALSE] else d
      }
      else if (tt == "fun_cog") {
        d <- as.data.frame(proj$functions$COG$abund)
        if (!is.null(smp)) d[, colnames(d) %in% smp, drop = FALSE] else d
      }
      else if (tt == "fun_kegg") {
        d <- as.data.frame(proj$functions$KEGG$abund)
        if (!is.null(smp)) d[, colnames(d) %in% smp, drop = FALSE] else d
      }
    }, error = function(e) NULL)
  })

  output$data_table <- renderDT({
    df <- get_table_data()
    req(df)
    datatable(
      df,
      options = list(
        pageLength = 20,
        scrollX    = TRUE,
        dom        = "lfrtip"
      ),
      class = "compact hover stripe"
    )
  })

  output$download_table <- downloadHandler(
    filename = function() paste0("sqm_", input$table_type, "_", Sys.Date(), ".csv"),
    content  = function(file) {
      df <- get_table_data()
      req(df)
      write.csv(df, file, row.names = TRUE)
    }
  )
}

shinyApp(ui = ui, server = server)
