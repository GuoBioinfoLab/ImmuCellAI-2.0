make_renderer <- function(analysis) {
  combined <- do.call(rbind, analysis[LETTERS[1:4]])
  s <- analysis$E
  method_labels <- c("ImmuCellAI 2.0", "BayesPrism (static ref.)",
                     "CIBERSORT (LM22)", "CITMIC", "DWLS", "ImmuCellAI", "MuSiC (music.basic)")
  targets <- target_labels <- list()
  for (dataset in datasets) {
    d <- combined[combined$dataset == dataset, ]
    types <- unique(d$target)
    labels <- setNames(d$target_label[match(types, d$target)], types)
    targets[[dataset]] <- types[order(tolower(unname(labels[types])), method = "radix")]
    target_labels[[dataset]] <- labels
  }
palette <- c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598",
             "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142")
ramp <- grDevices::colorRamp(palette)
color_r <- function(r) {
  rgb <- ramp(pmax(0, pmin(1, (r + 1) / 2)))
  grDevices::rgb(rgb[, 1], rgb[, 2], rgb[, 3], maxColorValue = 255)
}
u <- function(x) grid::unit(x, "native")
txt <- function(label, x, y, size = 11, color = "#171717", bold = FALSE,
                just = "left", rot = 0) {
  grid::grid.text(label, x = u(x), y = u(y), just = just, rot = rot,
            gp = grid::gpar(fontfamily = "Arial", fontsize = size, col = color,
                      fontface = if (bold) "bold" else "plain"))
}
rect <- function(x, y, w, h, fill = NA, border = NA, lwd = 1) {
  grid::grid.rect(x = u(x), y = u(y), width = u(w), height = u(h),
            gp = grid::gpar(fill = fill, col = border, lwd = lwd))
}
draw_sector <- function(x, y, r, status) {
  if (status != "valid") {
    txt(if (status == "unavailable") "/" else "x", x, y, size = 14,
        bold = status != "unavailable", just = "centre")
    return(invisible(NULL))
  }
  radius <- 0.17
  grid::grid.circle(x = u(x), y = u(y), r = u(radius),
              gp = grid::gpar(fill = "white", col = NA))
  if (abs(r) > 0) {
    a <- seq(pi / 2, pi / 2 - 2 * pi * abs(r),
             length.out = max(3, ceiling(240 * abs(r))))
    grid::grid.polygon(x = u(c(x, x + radius * cos(a), x)),
                 y = u(c(y, y + radius * sin(a), y)),
                 gp = grid::gpar(fill = color_r(r), col = "#252525", lwd = 0.6))
  }
  grid::grid.circle(x = u(x), y = u(y), r = u(radius),
              gp = grid::gpar(fill = NA, col = "#252525", lwd = 1.1))
}

panel_layout <- data.frame(
  dataset = datasets, letter = c("A", "B", "C", "D"),
  left = c(2.50, 2.50, 2.50, 10.50), top = c(16.00, 11.00, 6.15, 6.15),
  letter_x = c(0.35, 0.35, 0.35, 8.45),
  title = c("GSE164522 (CRC) | 26 targets", "GSE146771 (CRC) | 21 targets",
            "GSE176078 (BRCA) | 9 targets", "GSE162025 (NPC) | 11 targets"),
  stringsAsFactors = FALSE
)
draw_panel <- function(layout) {
  dataset <- layout$dataset
  top <- layout$top
  types <- targets[[dataset]]
  left <- layout$left
  dx <- 0.54
  dy <- 0.46
  width <- length(types) * dx
  height <- length(methods) * dy
  bottom <- top - height
  rect(left + width / 2, top + 0.19, width, 0.38, "#F7EAEA", "#202020", 1)
  txt(layout$title, left + width / 2, top + 0.19, 13, "#84252A", TRUE, "centre")
  txt(layout$letter, layout$letter_x, top + 0.19, 21, bold = TRUE)
  rect(left + width / 2, top - dy / 2, width, dy, "#FFF5F4")
  d <- combined[combined$dataset == dataset, ]
  row_labels <- method_labels
  if (dataset %in% c("GSE164522", "GSE146771")) {
    row_labels[methods == "DWLS"] <- "Custom weighted LM*"
  }
  for (i in seq_along(methods)) {
    y <- top - (i - 0.5) * dy
    txt(row_labels[i], left - 0.13, y, 11.5,
        if (i == 1L) "#92252D" else "#171717", i == 1L, "right")
    for (j in seq_along(types)) {
      z <- d[d$method == methods[i] & d$target == types[j], ]
      stopifnot(nrow(z) == 1L)
      draw_sector(left + (j - 0.5) * dx, y, z$Pearson, z$status)
    }
  }
  rect(left + width / 2, top - height / 2, width, height, NA, "#202020", 1.1)
  for (j in seq_along(types)) {
    label <- unname(target_labels[[dataset]][types[j]])
    txt(label, left + (j - 0.5) * dx, bottom - 0.12,
        11.5, just = "right", rot = 45)
  }
}

draw_legend <- function(x = 15.28, y = 9.54) {
  h <- 1.72
  w <- 0.24
  txt("Pearson r", x - 0.36, y + h / 2 + 0.27, 12, bold = TRUE)
  n <- 240L
  for (i in seq_len(n)) {
    r <- -1 + 2 * (i - 0.5) / n
    rect(x, y - h / 2 + (i - 0.5) * h / n, w, h / n + 0.002, color_r(r))
  }
  for (r in c(-1, -0.5, 0, 0.5, 1)) {
    yy <- y + r * h / 2
    grid::grid.lines(x = u(c(x + w / 2, x + w / 2 + 0.05)), y = u(c(yy, yy)),
               gp = grid::gpar(col = "#555555", lwd = 0.6))
    txt(sprintf("%.1f", r), x + w / 2 + 0.10, yy, 10)
  }
  txt("Sector area = |r|", x - 0.48, y - 1.15, 10.5)
  txt("Symbol key", x - 0.48, y - 1.54, 11.5, bold = TRUE)
  txt("x", x - 0.41, y - 1.91, 14, bold = TRUE, just = "centre")
  txt("Undefined r", x - 0.19, y - 1.91, 10.5)
  txt("/", x - 0.41, y - 2.25, 14, just = "centre")
  txt("Unavailable", x - 0.19, y - 2.25, 10.5)
}

draw_summary_panel <- function(dataset, left) {
  width <- 3.34
  limits <- c(-1, 1.25)
  xx <- function(v) left + width * (v - limits[1]) / diff(limits)
  z <- s[s$dataset == dataset, ]
  stopifnot(all(c(z$lower, z$upper) >= limits[1]),
            all(c(z$lower, z$upper) <= limits[2]))
  colors <- c("#A9103D", rep("#367D91", 6))
  labels <- c("ImmuCellAI 2.0", "BayesPrism", "CIBERSORT (LM22)", "CITMIC",
              "DWLS / custom LM*", "ImmuCellAI", "MuSiC (music.basic)")
  cohorts <- c(GSE164522 = "CRC", GSE146771 = "CRC",
               GSE176078 = "BRCA", GSE162025 = "NPC")
  line <- function(x, y, color, width) {
    grid::grid.lines(x = u(x), y = u(y), gp = grid::gpar(col = color, lwd = width))
  }
  txt(paste0(dataset, " (", cohorts[dataset], ")"), left, 6.74, 12.5, bold = TRUE)
  txt(paste0("Available targets: ", unique(z$n_candidate)), left, 6.45, 10, "#555555")
  line(rep(xx(0), 2), c(2.60, 6.13), "#C3C3C3", 0.7)
  line(c(left, left + width), c(2.60, 2.60), "#505050", 0.7)
  for (v in seq(-1, 1, by = 0.5)) {
    line(rep(xx(v), 2), c(2.60, 2.54), "#505050", 0.7)
    txt(sprintf("%.1f", v), xx(v), 2.36, 10.5, just = "centre")
  }
  txt("Pearson r", left + width / 2, 2.08, 12, just = "centre")
  for (i in seq_along(methods)) {
    row <- z[z$method == methods[i], ]
    y <- 6.10 - (i - 1) * 0.52
    line(xx(c(row$lower, row$upper)), c(y, y), colors[i], 1.8)
    grid::grid.points(x = u(xx(row$mean_Pearson)), y = u(y), pch = 18,
                size = grid::unit(0.165, "inches"), gp = grid::gpar(col = colors[i]))
    txt(sprintf("%.3f +/- %.3f  (%d/%d)", row$mean_Pearson, row$SD_Pearson,
                row$n, row$n_candidate), left + 0.045, y + 0.19, 9.5)
    if (dataset == datasets[1]) {
      txt(labels[i], left - 0.18, y, 11,
          if (i == 1L) "#912936" else "#252525", i == 1L, "right")
    }
  }
}

page_width <- 17.2
page_height <- 23.3
pie_offset <- 5.7
draw_figure <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale = c(0, page_width), yscale = c(0, page_height)))
  grid::pushViewport(grid::viewport(x = u(0), y = u(pie_offset), width = u(17.2), height = u(17.6),
                        just = c("left", "bottom"),
                        xscale = c(0, 17.2), yscale = c(0, 17.6)))
  txt("Immune-cell deconvolution: correlations across four cohorts", 0.38, 17.20,
      17, bold = TRUE)
  for (i in seq_len(nrow(panel_layout))) draw_panel(panel_layout[i, ])
  draw_legend()
  grid::popViewport()
  txt("E", 0.35, 7.23, 21, bold = TRUE)
  txt("Pearson correlation: mean +/- SD", 0.82, 7.23, 15, bold = TRUE)
  for (i in seq_along(datasets)) {
    draw_summary_panel(datasets[i], 2.30 + (i - 1) * 3.62)
  }
  grid::grid.lines(x = u(c(0.38, 16.80)), y = u(c(1.85, 1.85)),
             gp = grid::gpar(col = "#D0D0D0", lwd = 0.65))
  grid::popViewport()
}


  bounds <- function(panel) {
    if (panel == "complete") return(c(0, page_width, 0, page_height))
    if (panel == "E") return(c(0, page_width, 1.85, 7.65))
    p <- panel_layout[panel_layout$letter == panel, ]
    right <- p$left + length(targets[[p$dataset]]) * .54 + 2.2
    c(p$letter_x - .25, right, p$top - 7 * .46 - 1.6, p$top + .9)
  }
  draw <- function(panel) {
    if (panel == "complete") return(draw_figure())
    b <- bounds(panel)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(xscale = b[1:2], yscale = b[3:4]))
    on.exit(grid::popViewport(), add = TRUE)
    if (panel == "E") {
      txt("E", .35, 7.23, 21, bold = TRUE)
      txt("Pearson correlation: mean +/- SD", .82, 7.23, 15, bold = TRUE)
      for (i in seq_along(datasets)) draw_summary_panel(datasets[i], 2.30 + (i - 1) * 3.62)
    } else {
      p <- panel_layout[panel_layout$letter == panel, ]
      draw_panel(p)
      draw_legend(p$left + length(targets[[p$dataset]]) * .54 + .6, p$top - .5)
    }
  }
  list(draw = draw, size = function(panel) {
    b <- bounds(panel)
    c(diff(b[1:2]), diff(b[3:4]))
  })
}

save_panel <- function(renderer, panel, output_dir, output_format, dpi) {
  size <- renderer$size(panel)
  formats <- if (output_format == "both") c("pdf", "png") else output_format
  files <- file.path(output_dir, paste0("Fig2_", panel, ".", formats))
  for (i in seq_along(formats)) {
    if (formats[i] == "pdf") {
      grDevices::cairo_pdf(files[i], width = size[1], height = size[2],
                          family = "Arial", bg = "white", onefile = TRUE)
    } else {
      grDevices::png(files[i], width = size[1], height = size[2], units = "in",
                     res = dpi, type = "cairo", family = "Arial", bg = "white")
    }
    device <- grDevices::dev.cur()
    tryCatch(renderer$draw(panel), finally = grDevices::dev.off(device))
  }
  files
}
