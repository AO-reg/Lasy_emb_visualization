options(shiny.maxRequestSize = 1024^3)
options(shiny.sanitize.errors = FALSE)

library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)

# ---- ここで順序を指定（優先順）----
my_levels <- c(
  "64cells",
  "128cells(early gastrula)",
  "256cells(early gastrula)",
  "early blastula",
  "512cells(mid-blastula)",
  "mid-blastula",
  "late blastula",
  "advanced late blastula",
  "early gastrula",
  "slightly advanced early gastrula",
  "advanced early gastrula",
  "mid-gastrula",
  "advanced mid-gastrula",
  "late gastrula",
  "end of gastrula",
  "early neurula"
)

# -------- ユーティリティ --------
get_expr_mat <- function(obj, assay = "RNA", expr = "data") {
  z <- try(GetAssayData(obj, assay = assay, layer = expr), silent = TRUE)
  if (!inherits(z, "try-error")) return(z)
  GetAssayData(obj, assay = assay, slot = expr)
}

parse_genes <- function(s) {
  if (is.null(s) || !nzchar(s)) return(character(0))
  x <- unlist(strsplit(s, "[,、\\s]+", perl = TRUE))
  unique(x[nzchar(x)])
}

safe_meta <- function(meta, col, n) {
  if (col %in% colnames(meta)) meta[[col]] else rep(NA, n)
}

expand_genes_species <- function(gene_roots, obj_xl, obj_pw, obj_pweu,
                                 assay = "RNA", expr = "data",
                                 xl_suffixes = c(".L", ".S"),
                                 pw_transform = toupper) {
  mx <- get_expr_mat(obj_xl, assay = assay, expr = expr)
  mp <- get_expr_mat(obj_pw, assay = assay, expr = expr)
  me <- if (!is.null(obj_pweu)) get_expr_mat(obj_pweu, assay = assay, expr = expr) else NULL
  
  xl <- list(); pw <- list(); pweu <- list()
  for (gr in gene_roots) {
    xl_cand <- unique(c(gr, paste0(gr, xl_suffixes)))
    xl_hit  <- intersect(xl_cand, rownames(mx))
    if (length(xl_hit)) xl[[gr]] <- xl_hit
    
    pw_cand <- pw_transform(gr)
    pw_hit  <- intersect(pw_cand, rownames(mp))
    if (length(pw_hit)) pw[[gr]] <- pw_hit
    
    if (!is.null(me)) {
      eu_hit <- intersect(pw_cand, rownames(me))
      if (length(eu_hit)) pweu[[gr]] <- eu_hit
    }
  }
  list(xl = xl, pw = pw, pweu = pweu, mx = mx, mp = mp, me = me)
}

create_pseudotime_plot_multi <- function(
    gene_roots,
    obj_xl, obj_pw, obj_pweu = NULL,
    assay = "RNA", expr = "data",
    x_axis = c("Pseudotime", "name", "GMM"),
    pseudotime_col = "Pseudotime", stage_col = "stage",
    name_col = "name", gmm_col = "GMM",
    span = 0.9, point_size = 0.8,
    base_text_size = 14
) {
  x_axis <- match.arg(x_axis)
  expd <- expand_genes_species(gene_roots, obj_xl, obj_pw, obj_pweu, assay, expr)
  if (length(expd$xl) == 0 && length(expd$pw) == 0 && length(expd$pweu) == 0) {
    validate(need(FALSE, sprintf("指定遺伝子が見つかりません: %s", paste(gene_roots, collapse = ", "))))
  }
  
  make_df <- function(object, mat, genes_list, who, allow_axis) {
    if (length(genes_list) == 0) return(NULL)
    if (!(x_axis %in% allow_axis)) return(NULL)
    
    meta  <- object@meta.data
    ncell <- ncol(mat)
    xvec  <- switch(x_axis,
                    "Pseudotime" = safe_meta(meta, pseudotime_col, ncell),
                    "name"       = safe_meta(meta, name_col,       ncell),
                    "GMM"        = safe_meta(meta, gmm_col,        ncell))
    out <- list()
    for (gr in names(genes_list)) {
      for (g in genes_list[[gr]]) {
        vals  <- as.numeric(mat[g, , drop = TRUE])
        gtype <- if (who == "Xl") {
          if (grepl("\\.L$", g)) "Xl_L" else if (grepl("\\.S$", g)) "Xl_S" else "Xl_other"
        } else if (who == "Pw") {
          "Pw"
        } else {
          "Pw_EU"
        }
        out[[paste(who, g)]] <- data.frame(
          gene = vals, x = xvec,
          stage = safe_meta(meta, stage_col, ncell),
          name  = safe_meta(meta, name_col,  ncell),
          GMM   = safe_meta(meta, gmm_col,   ncell),
          gene_id = g, gene_root = gr, species = who, gene_type = gtype,
          stringsAsFactors = FALSE
        )
      }
    }
    dplyr::bind_rows(out)
  }
  
  dfx  <- make_df(obj_xl,   expd$mx,  expd$xl,   "Xl",   allow_axis = c("Pseudotime", "name", "GMM"))
  dfp  <- make_df(obj_pw,   expd$mp,  expd$pw,   "Pw",   allow_axis = c("Pseudotime", "name", "GMM"))
  dfeu <- if (!is.null(obj_pweu)) make_df(obj_pweu, expd$me, expd$pweu, "Pw_EU", allow_axis = c("name")) else NULL
  
  df <- dplyr::bind_rows(dfx, dfp, dfeu)
  validate(need(nrow(df) > 0, "可視化できるデータがありません"))
  
  present_roots <- unique(gene_roots[gene_roots %in% df$gene_root])
  if (length(present_roots) == 0) present_roots <- unique(df$gene_root)
  df$gene_root <- factor(df$gene_root, levels = present_roots)
  
  df$gene_type <- factor(df$gene_type,
                         levels = c("Xl_L", "Xl_S", "Xl_other", "Pw", "Pw_EU"))
  cols <- c("Xl_L"="#3366ff","Xl_S"="#99b3ff","Xl_other"="#669cff","Pw"="#ffa500","Pw_EU"="#ff9ac2")
  
  if (x_axis == "name") {
    # —— name 軸は my_levels を最優先に、未掲載カテゴリは末尾へ追加して固定順に
    levels_use <- c(my_levels, setdiff(unique(df$name), my_levels))
    df$x <- factor(df$name, levels = levels_use)
    
    p <- ggplot(df, aes(x = x, y = gene, fill = gene_type)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.7)) +
      geom_jitter(aes(color = gene_type),
                  size = point_size, alpha = 0.6,
                  position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7)) +
      scale_x_discrete(limits = levels_use, drop = FALSE) +
      scale_fill_manual(values = cols[levels(df$gene_type)], name = "Type") +
      scale_color_manual(values = cols[levels(df$gene_type)], guide = "none") +
      labs(x = "Stage (name)", y = "Expression") +
      theme_minimal(base_size = base_text_size) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.margin = margin(5.5, 5.5, 30, 5.5, "pt")
      ) +
      facet_wrap(~ gene_root, scales = "free_y", ncol = 2)
  } else {
    is_num <- is.numeric(df$x)
    p <- ggplot(df, aes(x = x, y = gene, color = gene_type)) +
      { if (is_num) geom_smooth(se = TRUE, span = span) } +
      geom_point(alpha = 0.6, size = point_size) +
      { if (!is_num && x_axis == "GMM") scale_x_discrete(drop = FALSE) } +
      scale_color_manual(values = cols[levels(df$gene_type)], name = "Type") +
      labs(x = x_axis, y = "Expression") +
      theme_minimal(base_size = base_text_size) +
      theme(
        axis.text.x = if (is_num) element_text() else element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.margin = margin(5.5, 5.5, 30, 5.5, "pt")
      ) +
      facet_wrap(~ gene_root, scales = "free_y", ncol = 2)
  }
  p
}

auto_height <- function(n_genes) {
  rows <- ceiling(max(1, n_genes) / 2)
  320 + 300 * rows
}

# -------- UI --------
ui <- fluidPage(
  titlePanel("Xl / Pw / Pw_EU 擬時間・カテゴリ可視化"),
  sidebarLayout(
    sidebarPanel(width = 3,
                 fileInput("xl_rds",  "Xl Seurat (.rds)",   accept = ".rds"),
                 fileInput("pw_rds",  "Pw Seurat (.rds)",   accept = ".rds"),
                 fileInput("pweu_rds","Pw EU Seurat (.rds)",accept = ".rds"),
                 hr(),
                 textInput("genes", "語根（, ・ 、 ・ 空白 ・ タブで複数可）", value = "fgf8, dkk1, gsc"),
                 selectInput("x_axis", "X軸", choices = c("Pseudotime","name","GMM"), selected = "name"),
                 textInput("assay", "Assay", value = "RNA"),
                 textInput("expr",  "Layer/slot（data / counts / scale.data）", value = "data"),
                 sliderInput("span", "平滑 span", min = 0.2, max = 1.5, step = 0.1, value = 0.9),
                 sliderInput("pt", "点サイズ", min = 0, max = 3, step = 0.1, value = 0.8),
                 sliderInput("basesize", "文字サイズ", min = 10, max = 24, step = 1, value = 14)
    ),
    mainPanel(width = 9,
              uiOutput("dyn_plot_holder")
    )
  )
)

# -------- Server --------
server <- function(input, output, session) {
  obj_xl   <- reactive({ req(input$xl_rds);   readRDS(input$xl_rds$datapath) })
  obj_pw   <- reactive({ req(input$pw_rds);   readRDS(input$pw_rds$datapath) })
  obj_pweu <- reactive({
    if (is.null(input$pweu_rds)) return(NULL)
    readRDS(input$pweu_rds$datapath)
  })
  
  roots <- reactive(parse_genes(input$genes))
  
  output$dyn_plot_holder <- renderUI({
    h <- auto_height(length(roots()))
    plotOutput("p_gene", height = paste0(h, "px"))
  })
  
  output$p_gene <- renderPlot({
    xl   <- obj_xl()
    pw   <- obj_pw()
    pweu <- obj_pweu()
    rs   <- roots()
    validate(need(length(rs) > 0, "語根を入力してください"))
    
    p <- create_pseudotime_plot_multi(
      gene_roots = rs,
      obj_xl = xl, obj_pw = pw, obj_pweu = pweu,
      assay = input$assay, expr = input$expr,
      x_axis = input$x_axis,
      span = input$span, point_size = input$pt,
      base_text_size = input$basesize
    )
    p
  }, res = 96)
}

shinyApp(ui, server)
