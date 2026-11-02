make_renderer <- function(analysis) {
methods <- c("ImmuCellAI 2.0", "BayesPrism", "CIBERSORT", "CITMIC", "DWLS", "ImmuCellAI", "MuSiC")
groups <- c("BostonGene blood", "Newman PBMC", "Monaco PBMC")
counts <- setNames(c(26L, 6L, 5L), groups)
panels <- setNames(c("A", "B", "C"), groups)
display <- c(Nonplasma_B = "Nonplasma B", Plasma_lineage = "Plasma lineage",
             T_total = "T total", DC_total = "DC total")
colors <- setNames(c("#DF5A5A", "#39959D", "#AD9C85", "#96A5A5", "#6F9FB8", "#8F6EA5", "#E7BB53"), methods)
spectral <- c("#5e4fa2", "#3288bd", "#66c2a5", "#abdda4", "#e6f598", "#ffffbf",
              "#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#9e0142")
ink <- "#202020"
red <- "#8C252D"
width <- 1200

  records <- do.call(rbind, analysis[LETTERS[1:3]])
  summary <- analysis$D
  runtimes <- analysis$E
u <- function(x) grid::unit(x, "native")
qa_boxes <- list()
render_height <- 0
ramp <- grDevices::colorRamp(spectral)
corr_color <- function(r) {
  rgb <- ramp(pmax(0, pmin(1, (r + 1) / 2)))
  grDevices::rgb(rgb[, 1], rgb[, 2], rgb[, 3], maxColorValue = 255)
}
text <- function(label, x, y, size = 12, bold = FALSE, color = ink,
                 align = "left", angle = 0) {
  label <- as.character(label)
  gp <- grid::gpar(fontfamily = "Arial", fontsize = size, col = color,
                   fontface = if (bold) "bold" else "plain")
  tw <- grid::convertWidth(grid::grobWidth(grid::textGrob(label, gp = gp)), "inches", valueOnly = TRUE) * 72
  dx <- switch(align, center = -tw / 2, right = -tw, 0)
  bounds <- expand.grid(bx = c(dx, dx + tw), by = c(-size * .25, size))
  theta <- angle * pi / 180
  bx <- x + bounds$bx * cos(theta) - bounds$by * sin(theta)
  by <- y + bounds$bx * sin(theta) + bounds$by * cos(theta)
  qa_boxes[[length(qa_boxes) + 1L]] <<- c(min(bx), min(by), max(bx), max(by))
  grid::grid.text(label, x = u(x), y = u(y),
                  just = c(if (align == "center") "centre" else align, "bottom"),
                  rot = angle, gp = gp)
}
line <- function(x1, y1, x2, y2, color = ink, width = 1, dash = FALSE) {
  grid::grid.lines(u(c(x1, x2)), u(c(y1, y2)),
                   gp = grid::gpar(col = color, lwd = width, lty = if (dash) "44" else "solid"))
}
rect <- function(x, y, w, h, fill = NA, stroke = ink, lw = 1) {
  grid::grid.rect(u(x), u(y), width = u(w), height = u(h), just = c("left", "bottom"),
                  gp = grid::gpar(fill = fill, col = stroke, lwd = lw))
}
circle_sector <- function(x, y, radius, r) {
  grid::grid.circle(u(x), u(y), r = u(radius),
                    gp = grid::gpar(fill = "white", col = ink, lwd = 1.1))
  if (r == 0) return(invisible(NULL))
  n <- max(3L, ceiling(100 * abs(r)))
  theta <- pi / 2 - 2 * pi * abs(r) * (0:n) / n
  grid::grid.polygon(u(c(x, x + radius * cos(theta))), u(c(y, y + radius * sin(theta))),
                     gp = grid::gpar(fill = corr_color(r), col = ink, lwd = 1.1))
}
title_strip <- function(x, top, w, title, size = 17) {
  rect(x, top, w, 34, "#FAEBEF", lw = 1.3)
  text(title, x + w / 2, top + 11, size, TRUE, red, "center")
}
method_label <- function(method) if (method == "MuSiC") "MuSiC (basic)" else method
target_label <- function(targets) {
  labels <- targets
  idx <- match(targets, names(display))
  labels[!is.na(idx)] <- unname(display[idx[!is.na(idx)]])
  labels
}
pie_panel <- function(dataset, x, y, w, h, title, labels = TRUE, font = 14, title_size = 16) {
  rows <- records[records$Dataset == dataset, ]
  targets <- unique(rows$cellType)
  targets <- targets[order(tolower(target_label(targets)), method = "radix")]
  stopifnot(length(targets) == counts[[dataset]])
  rect(x, y, w, h, "#FFFFFF", lw = 1.3)
  title_strip(x, y + h, w, title, title_size)
  dx <- w / length(targets)
  dy <- h / length(methods)
  for (j in seq_along(methods)) {
    method <- methods[j]
    cy <- y + h - (j - .5) * dy
    if (labels) text(method_label(method), x - 13, cy - 5, font,
                     j == 1L, if (j == 1L) red else ink, "right")
    for (k in seq_along(targets)) {
      row <- rows[rows$method == method & rows$cellType == targets[k], ]
      stopifnot(nrow(row) == 1L)
      cx <- x + (k - .5) * dx
      if (is.finite(row$r)) {
        circle_sector(cx, cy, min(12.5, dx * .33), row$r)
      } else {
        symbol <- if (row$status %in% c("constant_prediction", "constant_truth", "insufficient_pairs")) "x" else "/"
        text(symbol, cx, cy - 5, 16, TRUE, align = "center")
      }
    }
  }
  for (k in seq_along(targets)) {
    label <- sub("Neutrophils LD", "Neutrophils", target_label(targets[k]), fixed = TRUE)
    text(label, x + (k - .5) * dx, y - 16, font, align = "right", angle = 48)
  }
}
key <- function() {
  x <- 1120; y <- 1060; w <- 15; h <- 125
  text("Pearson r", x - 8, y + h + 15, 12, TRUE)
  for (i in 0:149) rect(x, y + h * i / 150, w, h / 150 + .3, corr_color(-1 + 2 * i / 149), NA)
  for (value in c(-1, -.5, 0, .5, 1)) {
    text(format(value, trim = TRUE), x + 23, y + (value + 1) * h / 2 - 3, 10)
  }
  text("Symbol key", 1103, 1024, 12, TRUE, red)
  text("x", 1107, 1001, 13, TRUE)
  text("Undefined r", 1120, 1001, 9)
  text("/", 1107, 977, 13, TRUE)
  text("Unavailable", 1120, 981, 9)
  text("or no reference", 1120, 967, 9)
  text("Sector area: |r|", 1103, 941, 9)
}
runtime_panel <- function() {
  x <- 910; y <- 549; w <- 245; h <- 251
  title_strip(x, y + h, w, "Runtime efficiency", 16)
  rect(x, y, w, h, "#FFFFFF", lw = 1.3)
  for (tick in c(0, .5, 1, 1.5, 2)) {
    xx <- x + w * tick / 2.4
    line(xx, y, xx, y + h, "#D9D9D9", .8, TRUE)
    text(sprintf("%.1f", tick), xx, y - 17, 10, align = "center")
  }
  dy <- h / length(methods)
  scenarios <- c("single_sample_1core", "fifty_samples_8cores")
  for (j in seq_along(methods)) {
    method <- methods[j]
    cy <- y + h - (j - .5) * dy
    label <- if (method %in% c("CIBERSORT", "DWLS")) paste0(method, "*") else method
    text(label, x - 12, cy - 4, 12.5, j == 1L, if (j == 1L) red else ink, "right")
    for (k in seq_along(scenarios)) {
      row <- runtimes[runtimes$MethodShort == method & runtimes$Scenario == scenarios[k], ]
      stopifnot(nrow(row) == 1L)
      rect(x, cy - if (k == 1L) 12 else 0, row$LogTime / 2.4 * w, 12,
           c("#379DA5", "#E96061")[k], lw = .8)
    }
  }
  text("log10(seconds per sample + 1)", x + w / 2, y - 36, 11, align = "center")
  for (i in 1:2) {
    yy <- y - 60 - 18 * (i - 1)
    rect(x + 15, yy, 10, 10, c("#379DA5", "#E96061")[i], lw = .5)
    text(c("1 sample / 1 core", "50 samples / 8 cores")[i], x + 33, yy, 11)
  }
}
summary_row <- function() {
  text("D", 29, 414, 23, TRUE)
  text("Pearson Correlation Across Immune Cell Targets (Mean \u00b1 SD)", 74, 416, 18, TRUE)
  labels <- c("BostonGene Blood Cohort", "Newman's Cohort (Whole Blood)", "Monaco's Cohort (PBMC)")
  ymin <- min(0, floor(min(summary$lower) * 4) / 4)
  ymax <- ceiling((max(summary$upper) + .12) * 4) / 4
  for (i in seq_along(groups)) {
    dataset <- groups[i]
    x <- 66 + (i - 1) * 378; y <- 170; w <- 344; h <- 194
    fy <- function(value) y + h * (value - ymin) / (ymax - ymin)
    rect(x, y, w, h, "#FFFFFF", lw = .85)
    title_strip(x, y + h, w, labels[i], 13)
    for (tick in seq(ymin, ymax, by = .25)) {
      yy <- fy(tick)
      line(x, yy, x + w, yy, "#D9D9D9", .65, TRUE)
      text(sprintf("%.2f", tick), x - 9, yy - 3, 9, align = "right")
    }
    line(x, fy(0), x + w, fy(0), "#8B8B8B", .8)
    rect(x, y, w, h, lw = .85)
    if (i == 1L) text("Mean Pearson r", 24, y + h / 2, 12.5, TRUE, align = "center", angle = 90)
    for (j in seq_along(methods)) {
      method <- methods[j]
      row <- summary[summary$Dataset == dataset & summary$method == method, ]
      cx <- x + w * (j - .5) / length(methods)
      bottom <- min(fy(0), fy(row$mean))
      rect(cx - 15, bottom, 30, abs(fy(row$mean) - fy(0)), colors[[method]], lw = .8)
      low <- fy(row$lower); high <- fy(row$upper)
      line(cx, low, cx, high, width = 1.25)
      line(cx - 5, low, cx + 5, low, width = 1.25)
      line(cx - 5, high, cx + 5, high, width = 1.25)
      text(paste0(row$n_valid, "/", row$n_candidates), cx, high + 8, 12.5, align = "center")
      text(method_label(method), cx + 3, y - 15, 10.8, align = "right", angle = 48)
    }
  }
}
draw_figure <- function(height, full) {
  qa_boxes <<- list()
  render_height <<- height
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale = c(0, width), yscale = c(0, height)))
  on.exit(grid::popViewport(), add = TRUE)
  if (full) {
    pie_panel(groups[1], 148, 930, 937, 280,
              "BostonGene Blood Cohort (n = 45)", font = 13, title_size = 18)
    text("A", 29, 1220, 25, TRUE)
    key()
    pie_panel(groups[2], 149, 540, 282, 260,
              "Newman's Cohort (Whole Blood, n = 12)", font = 14, title_size = 12.5)
    text("B", 29, 816, 25, TRUE)
    pie_panel(groups[3], 491, 540, 242, 260,
              "Monaco's Cohort (PBMC, n = 6)", labels = FALSE, font = 14, title_size = 13)
    text("C", 455, 816, 25, TRUE)
    runtime_panel()
    text("E", 754, 816, 25, TRUE)
  }
  summary_row()
  boxes <- do.call(rbind, qa_boxes)
  stopifnot(all(boxes[, 1] >= 0), all(boxes[, 2] >= 0),
            all(boxes[, 3] <= width), all(boxes[, 4] <= render_height))
}

  bounds <- function(panel) {
    switch(panel, complete = c(0, width, 0, 1280), A = c(0, width, 780, 1280),
           B = c(0, 510, 10, 470), C = c(0, 510, 10, 470),
           D = c(0, width, 30, 452), E = c(750, 1190, 440, 865))
  }
  draw <- function(panel) {
    if (panel == "complete") return(draw_figure(1280, TRUE))
    b <- bounds(panel)
    qa_boxes <<- list()
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(xscale = b[1:2], yscale = b[3:4]))
    on.exit(grid::popViewport(), add = TRUE)
    if (panel == "A") {
      pie_panel(groups[1], 148, 930, 937, 280,
                "BostonGene Blood Cohort (n = 45)", font = 13, title_size = 18)
      text("A", 29, 1220, 25, TRUE)
      key()
    } else if (panel %in% c("B", "C")) {
      index <- match(panel, c("B", "C")) + 1L
      title <- c("Newman's Cohort (Whole Blood, n = 12)",
                 "Monaco's Cohort (PBMC, n = 6)")[index - 1L]
      pie_panel(groups[index], 149, 140, 330, 260, title, font = 14, title_size = 13)
      text(panel, 29, 416, 25, TRUE)
    } else if (panel == "D") {
      summary_row()
    } else {
      runtime_panel()
      text("E", 754, 816, 25, TRUE)
    }
    boxes <- do.call(rbind, qa_boxes)
    stopifnot(all(boxes[, 1] >= b[1]), all(boxes[, 3] <= b[2]),
              all(boxes[, 2] >= b[3]), all(boxes[, 4] <= b[4]))
  }
  list(draw = draw, size = function(panel) {
    b <- bounds(panel)
    c(diff(b[1:2]), diff(b[3:4])) / 72
  })
}

save_panel <- function(renderer, panel, output_dir, output_format, dpi) {
  size <- renderer$size(panel)
  formats <- if (output_format == "both") c("pdf", "png") else output_format
  files <- file.path(output_dir, paste0("Fig3_", panel, ".", formats))
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
