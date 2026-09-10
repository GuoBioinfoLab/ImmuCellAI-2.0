# Original per-cohort pie plotting implementation; Figure 2 defaults to three cohorts.
# Override input_files and out_dir before sourcing to render additional metric tables.
options(stringsAsFactors = FALSE)

if (!exists("out_dir", inherits = FALSE)) {
  out_dir <- file.path("benchmarks", "output", "fig2", "plots")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!exists("input_files", inherits = FALSE)) {
  input_files <- data.frame(
    Dataset = c("GSE164522", "GSE176078", "GSE146771"),
    File = file.path("benchmarks", "output", "fig2",
                     c("GSE164522", "GSE176078", "GSE146771"),
                     "all_methods_per_celltype_metrics.txt"),
    stringsAsFactors = FALSE
  )
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package ggplot2 is required.")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  stop("Package RColorBrewer is required.")
}
if (!requireNamespace("magrittr", quietly = TRUE)) {
  stop("Package magrittr is required.")
}
if (!requireNamespace("stringr", quietly = TRUE)) {
  stop("Package stringr is required.")
}
library(ggplot2)
library(magrittr)
if (requireNamespace("scRNAtoolVis", quietly = TRUE)) {
  library(scRNAtoolVis)
}
native_geom_jjpie <- exists("geom_jjpie", mode = "function")

mycol7 <- c("#F66463", "#379DA5", "#6BA5C3", "#FAC74C", "#95A8AC", "#AE997E")
mycol7_scale <- scale_fill_manual(values = mycol7)

method_labels <- c(
  ImmuCellAI2_tcell_VB = "ImmuCellAI 2.0",
  ImmuCellAI2_flat_VB = "ImmuCellAI 2.0",
  ImmuCellAI2_flat_VB_UNKNOWN = "ImmuCellAI 2.0",
  ImmuCellAI2_tcell_VB_UNKNOWN = "ImmuCellAI 2.0",
  ImmuCellAI2_best_ImmuCellAI2_flat_VB = "ImmuCellAI 2.0",
  ImmuCellAI2_best_ImmuCellAI2_flat_VB_UNKNOWN = "ImmuCellAI 2.0",
  BayesPrism_first_state_chain600_burn500 = "BayesPrism",
  DWLS_weighted_lm = "DWLS",
  MuSiC_basic = "MuSiC",
  CITMIC_native = "CITMIC",
  CIBERSORT_default_LM22 = "CIBERSORT",
  CIBERSORT_nuSVR_ref_top1000 = "CIBERSORT",
  ImmuCellAI_native = "ImmuCellAI"
)

# Accept archived result identifiers while displaying the current method name.
legacy_labels <- method_labels
names(legacy_labels) <- gsub("ImmuCellAI2", "ImmuneHier", names(legacy_labels), fixed = TRUE)
names(legacy_labels) <- gsub("_tcell_", "_tcell_only_", names(legacy_labels), fixed = TRUE)
method_labels <- c(method_labels, legacy_labels)

method_order <- c(
  "ImmuCellAI 2.0", "BayesPrism", "CIBERSORT", "CITMIC",
  "DWLS", "ImmuCellAI", "MuSiC"
)

read_one_metrics <- function(dataset, file) {
  if (!file.exists(file)) stop("Input file not found: ", file)
  x <- read.delim(file, check.names = FALSE)
  required <- c("Method", "CellType", "PredictedStates", "Pearson")
  missing_cols <- setdiff(required, colnames(x))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", file, ": ", paste(missing_cols, collapse = ", "))
  }
  x$Dataset <- dataset
  x$Method <- as.character(x$Method)
  x$CellType <- as.character(x$CellType)
  x$PredictedStates <- as.character(x$PredictedStates)
  x$Pearson <- suppressWarnings(as.numeric(x$Pearson))
  x$Spearman <- if ("Spearman" %in% colnames(x)) suppressWarnings(as.numeric(x$Spearman)) else NA_real_
  x$method <- ifelse(x$Method %in% names(method_labels), unname(method_labels[x$Method]), x$Method)
  x$cellType <- x$CellType
  x$correlation <- x$Pearson
  x$status <- "valid"
  x$status[is.na(x$correlation) & trimws(x$PredictedStates) == ""] <- "not_available"
  x$status[is.na(x$correlation) & trimws(x$PredictedStates) != ""] <- "constant_prediction"
  x$correlation[!is.finite(x$correlation)] <- NA_real_
  keep_cols <- c(
    "Dataset", "Method", "CellType", "PredictedStates", "Pearson", "Spearman",
    "method", "cellType", "correlation", "status"
  )
  x[, keep_cols, drop = FALSE]
}

plot_data <- do.call(
  rbind,
  Map(read_one_metrics, input_files$Dataset, input_files$File)
)

drop_celltypes_by_dataset <- data.frame(
  Dataset = c(rep("GSE176078", 3), "GSE164522"),
  cellType = c("NKT", "Th1", "Th17", "NKT"),
  stringsAsFactors = FALSE
)

plot_data <- plot_data[
  !(
    plot_data$Dataset %in% drop_celltypes_by_dataset$Dataset &
      paste(plot_data$Dataset, plot_data$cellType) %in%
        paste(drop_celltypes_by_dataset$Dataset, drop_celltypes_by_dataset$cellType)
  ),
  ,
  drop = FALSE
]

gse107011_label_map <- c(
  "T CD8 Memory" = "CD8Tmemory",
  "T Memory" = "Tmemory",
  "Monocytes C" = "CMonocyte",
  "Monocytes I" = "IMonocyte",
  "T CD4 Naive" = "CD4Tn",
  "Basophils LD" = "Basophil",
  "T CD8 TE" = "CD8Temra",
  "Granulocytes LD" = "Granulocyte",
  "memory B cell" = "Bmemory",
  "T gd" = "gdT",
  "Plasmablasts" = "plasma",
  "T helper" = "T helper",
  "T CD8 Naive" = "CD8Tn"
)

idx_gse107011 <- plot_data$Dataset == "GSE107011" &
  plot_data$cellType %in% names(gse107011_label_map)
plot_data$cellType[idx_gse107011] <- unname(gse107011_label_map[plot_data$cellType[idx_gse107011]])

plot_data$method <- factor(plot_data$method, levels = rev(method_order))

write.csv(
  plot_data,
  file.path(out_dir, "combined_datasets_7tools_pie_plot_data.csv"),
  row.names = FALSE,
  quote = FALSE
)

status_count <- as.data.frame(table(plot_data$Dataset, plot_data$method, plot_data$status))
colnames(status_count) <- c("Dataset", "method", "status", "count")
write.csv(
  status_count,
  file.path(out_dir, "combined_datasets_7tools_status_count.csv"),
  row.names = FALSE,
  quote = FALSE
)

make_circle_poly <- function(cx, cy, r, n = 96) {
  a <- seq(0, 2 * pi, length.out = n + 1)
  data.frame(x = cx + r * cos(a), y = cy + r * sin(a))
}

make_wedge_poly <- function(cx, cy, r, frac, start = pi / 2, n = 96) {
  frac <- min(max(abs(frac), 0), 1)
  if (!is.finite(frac) || frac <= 0) return(NULL)
  a <- seq(start, start - 2 * pi * frac, length.out = max(3, ceiling(n * frac) + 1))
  data.frame(x = c(cx, cx + r * cos(a), cx), y = c(cy, cy + r * sin(a), cy))
}

if (!native_geom_jjpie) {
  geom_jjpie <- function(data, mapping = NULL, width = 1.2, ...) {
    if (is.null(data) || nrow(data) == 0) return(list())
    df <- data
    if (!"x_id" %in% colnames(df)) df$x_id <- as.numeric(df$cellType)
    if (!"y_id" %in% colnames(df)) df$y_id <- as.numeric(df$method)
    if (!"correlation" %in% colnames(df) && "Pearson" %in% colnames(df)) {
      df$correlation <- suppressWarnings(as.numeric(df$Pearson))
    }
    r <- min(width / 2, 0.48)
    bg_list <- vector("list", nrow(df))
    wedge_list <- list()
    for (i in seq_len(nrow(df))) {
      bg <- make_circle_poly(df$x_id[i], df$y_id[i], r = r)
      bg$idx <- i
      bg_list[[i]] <- bg
      if (is.finite(df$correlation[i])) {
        wg <- make_wedge_poly(df$x_id[i], df$y_id[i], r = r, frac = min(abs(df$correlation[i]), 1))
        if (!is.null(wg)) {
          wg$idx <- i
          wg$correlation <- df$correlation[i]
          wedge_list[[length(wedge_list) + 1]] <- wg
        }
      }
    }
    bg_poly <- do.call(rbind, bg_list)
    wedge_poly <- if (length(wedge_list) > 0) do.call(rbind, wedge_list) else data.frame()
    list(
      geom_polygon(
        data = bg_poly,
        aes(x = x, y = y, group = idx),
        inherit.aes = FALSE,
        fill = "white",
        color = "#1f1f1f",
        linewidth = 0.65
      ),
      geom_polygon(
        data = wedge_poly,
        aes(x = x, y = y, group = idx, fill = correlation),
        inherit.aes = FALSE,
        color = "#1f1f1f",
        linewidth = 0.55
      )
    )
  }
}

theme_blue <- theme(
  plot.title = element_text(
    size = 13,
    face = "bold",
    color = "darkred",
    hjust = 0,
    lineheight = 1.2
  ),
  plot.subtitle = element_text(
    size = 13,
    face = "bold",
    color = "grey30",
    lineheight = 1.2,
    hjust = 0
  ),
  panel.background = element_rect(fill = "white"),
  panel.grid.major.y = element_line(
    colour = "gray80",
    linewidth = 1,
    linetype = "dashed"
  ),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(
    vjust = 1,
    face = "bold",
    size = 14,
    color = "darkred"
  ),
  axis.title.y = element_text(
    size = 14,
    face = "bold",
    color = "darkred"
  ),
  axis.text.x = element_text(
    size = 12,
    colour = "black"
  ),
  legend.title = element_text(size = 12, colour = "black"),
  axis.text.y = element_text(size = 12, colour = "black"),
  legend.text = element_text(size = 12, colour = "black"),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
  legend.key = element_blank(),
  strip.background = element_rect(fill = "#FBE4E7", color = "black", linewidth = 1),
  strip.text = element_text(size = 18, colour = "black")
)

draw_one_dataset <- function(df, dataset) {
  df <- df[df$Dataset == dataset, , drop = FALSE]
  df <- df[!is.na(df$method), , drop = FALSE]
  df$method <- factor(as.character(df$method), levels = rev(method_order))
  df <- df[!is.na(df$method), , drop = FALSE]
  df$cellTypeLabel <- gsub("γδT", "gdT", df$cellType, fixed = TRUE)
  pie_width <- 0.65
  cell_width_in <- 0.88
  method_height_in <- 0.84
  outline_size <- 2.2
  outline_stroke <- if (identical(dataset, "GSE107011")) 0.22 else 1.05
  wedge_linewidth <- if (identical(dataset, "GSE107011")) 0.18 else 0.55
  undefined_x_stroke <- if (identical(dataset, "GSE107011")) 0.55 else 1.8
  undefined_x_size <- 1.5
  undefined_slash_size <- 2.8
  undefined_slash_half_width <- 0.13

  cell_order <- unique(as.character(df$cellTypeLabel))
  cell_order <- cell_order[order(tolower(cell_order))]
  df$cellTypeLabel <- factor(df$cellTypeLabel, levels = cell_order)
  df$x_id <- as.numeric(df$cellTypeLabel)
  df$y_id <- as.numeric(df$method)
  df <- df[order(df$x_id, df$y_id), , drop = FALSE]

  constant_df <- df[df$status == "constant_prediction", , drop = FALSE]
  not_available_df <- df[df$status == "not_available", , drop = FALSE]

  slash_df <- data.frame()
  if (nrow(not_available_df) > 0) {
    slash_df <- do.call(rbind, lapply(seq_len(nrow(not_available_df)), function(i) {
      x <- not_available_df$x_id[i]
      y <- not_available_df$y_id[i]
      data.frame(
        x = x - undefined_slash_half_width,
        y = y - undefined_slash_half_width,
        xend = x + undefined_slash_half_width,
        yend = y + undefined_slash_half_width,
        idx = i
      )
    }))
  }

  axis_x <- data.frame(x = seq_along(cell_order), label = cell_order)
  axis_y <- data.frame(y = seq_along(levels(df$method)), label = levels(df$method))
  legend_status <- data.frame(
    x = c(Inf, Inf),
    y = c(Inf, Inf),
    status = factor(c("constant prediction", "not available"),
                    levels = c("constant prediction", "not available"))
  )

  combined_data <- df

  if (native_geom_jjpie) {
    combined_data$cellType <- factor(as.character(combined_data$cellTypeLabel), levels = cell_order)
    combined_data$method <- factor(as.character(combined_data$method), levels = method_order)
    constant_native <- combined_data[combined_data$status == "constant_prediction", , drop = FALSE]
    not_available_native <- combined_data[combined_data$status == "not_available", , drop = FALSE]

    p <- ggplot(combined_data, aes(x = cellType, y = method, fill = correlation)) +
      geom_point(
        data = subset(combined_data, !is.na(correlation)),
        aes(x = cellType, y = method),
        inherit.aes = FALSE,
        shape = 21,
        size = outline_size,
        fill = "white",
        color = "#1f1f1f",
        stroke = outline_stroke
      ) +
      geom_jjpie(
        data = subset(combined_data, !is.na(correlation)),
        aes(piefill = correlation),
        width = pie_width,
        color = "#1f1f1f",
        linewidth = wedge_linewidth
      ) +
      geom_point(
        data = constant_native,
        aes(x = cellType, y = method),
        inherit.aes = FALSE,
        shape = 4,
        size = undefined_x_size,
        color = "black",
        stroke = undefined_x_stroke
      ) +
      geom_text(
        data = not_available_native,
        aes(x = cellType, y = method),
        inherit.aes = FALSE,
        label = "/",
        color = "black",
        size = undefined_slash_size,
        fontface = "bold"
      ) +
      facet_grid(. ~ Dataset, scales = "free", space = "free") +
      scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "Spectral") %>% rev(), limits = c(-1, 1)) +
      scale_x_discrete(expand = expansion(add = c(0.55, 0.55))) +
      scale_y_discrete(
        limits = rev(method_order),
        expand = expansion(add = c(0.55, 0.55))
      ) +
      labs(x = NULL, y = NULL) +
      theme_blue +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black")
      )
  } else {
    p <- ggplot(combined_data, aes(x = x_id, y = y_id, fill = correlation)) +
      geom_jjpie(
        data = subset(combined_data, !is.na(correlation)),
        aes(piefill = correlation),
        width = pie_width
      ) +
      geom_point(data = constant_df, aes(x = x_id, y = y_id),
                 inherit.aes = FALSE,
                 shape = 4, size = undefined_x_size, color = "black", stroke = 1.8) +
      geom_segment(data = slash_df, aes(x = x, y = y, xend = xend, yend = yend, group = idx),
                   inherit.aes = FALSE,
                   color = "black", linewidth = 0.95, lineend = "round") +
      geom_point(data = legend_status, aes(x = x, y = y, shape = status),
                 inherit.aes = FALSE, alpha = 0) +
      facet_grid(. ~ Dataset, scales = "free", space = "free") +
      scale_fill_gradientn(
        colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
        limits = c(-1, 1),
        name = "Pearson"
      ) +
      scale_shape_manual(
        values = c(`constant prediction` = "x", `not available` = "/"),
        name = "Undefined"
      ) +
      guides(shape = guide_legend(override.aes = list(alpha = 1, color = "black", size = 3))) +
      scale_x_continuous(
        breaks = axis_x$x,
        labels = axis_x$label,
        expand = expansion(mult = c(0.01, 0.01))
      ) +
      scale_y_continuous(
        breaks = axis_y$y,
        labels = axis_y$label,
        expand = expansion(mult = c(0.08, 0.08))
      ) +
      coord_fixed(ratio = 1) +
      labs(x = NULL, y = NULL, caption = "Undefined glyphs: x = constant prediction; slash = not available") +
      theme_blue +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black")
      )
  }

  n_cells <- length(cell_order)
  n_methods <- length(levels(df$method))
  if (identical(dataset, "GSE107011")) {
    p <- p + theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
      strip.background = element_rect(fill = "#FBE4E7", color = "black", linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.15)
    )
  }
  width_high <- n_cells * cell_width_in + 5.2
  width_low <- width_high
  height_high <- n_methods * method_height_in + 3.0
  height_low <- height_high
  if (identical(dataset, "GSE107011")) {
    width_high <- 7.2
    width_low <- 7.2
    height_high <- 4.2
    height_low <- 4.2
  }
  prefix <- file.path(out_dir, paste0(dataset, "_7tools_Figure2A_pie"))

  ggsave(file = paste0(prefix, "_high.pdf"), plot = p, width = width_high, height = height_high, limitsize = FALSE)
  ggsave(file = paste0(prefix, "_low.pdf"), plot = p, width = width_low, height = height_low, limitsize = FALSE)
  ggsave(file = paste0(prefix, "_high.png"), plot = p, width = width_high, height = height_high, dpi = 220, limitsize = FALSE)

  write.csv(
    df,
    file.path(out_dir, paste0(dataset, "_7tools_Figure2A_pie_plot_data.csv")),
    row.names = FALSE,
    quote = FALSE
  )

  invisible(p)
}

plots <- lapply(input_files$Dataset, function(ds) draw_one_dataset(plot_data, ds))
names(plots) <- input_files$Dataset

message("Saved output directory: ", out_dir)
message("Saved combined plot data: ", file.path(out_dir, "combined_datasets_7tools_pie_plot_data.csv"))
message("Saved status counts: ", file.path(out_dir, "combined_datasets_7tools_status_count.csv"))
