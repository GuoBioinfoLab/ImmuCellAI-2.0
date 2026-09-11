source("analysis/common/age_groups.R")
source("analysis/prepare_source_archive.R")

ages <- c(-1, 0, 1, 1.9, 2, 10, 10.1, 20, 20.1, 30, 30.1, 49.9, 50, 69.9, 70, 95, NA, Inf)
expected <- c(NA, "0-1", "0-1", "0-1", NA, NA, "10-20", "10-20", "20-30",
              "20-30", "30-50", "30-50", "50-70", "50-70", "70+", "70+", NA, NA)
stopifnot(identical(as.character(figure5_age_group(ages)), expected))
cells <- unlist(figure5_major_cells, use.names = FALSE)
stopifnot(length(cells) == 53L, !anyDuplicated(cells), length(figure5_major_cells) == 8L)
x <- matrix(seq_len(106), nrow = 2, dimnames = list(c("a", "b"), cells))
aggregated <- vapply(figure5_major_cells, function(z) rowSums(x[, z, drop = FALSE]), numeric(2))
stopifnot(identical(unname(rowSums(x)), unname(rowSums(aggregated))))

assert_error <- function(expr, pattern) {
  msg <- tryCatch({force(expr); NA_character_}, error = conditionMessage)
  stopifnot(!is.na(msg), grepl(pattern, msg, fixed = TRUE))
}
fixture <- tempfile("archive-input-")
dir.create(file.path(fixture, "module", "source_archive"), recursive = TRUE)
source_file <- file.path(fixture, "module", "source_archive", "example.R")
writeLines(c('x <- "<LOCAL_R_ROOT>/Fig5/input.tsv"',
             'stop("This must never execute during preparation")'), source_file)
dest <- tempfile("archive-output-")
root <- normalizePath(tempdir(), winslash = "/")
prepared <- prepare_source_archive(fixture, dest, c(LOCAL_R_ROOT = root))
stopifnot(length(prepared) == 1L, file.exists(prepared))
stopifnot(any(grepl(root, readLines(prepared), fixed = TRUE)))
stopifnot(any(grepl("<LOCAL_R_ROOT>", readLines(source_file), fixed = TRUE)))
assert_error(prepare_source_archive(fixture, dest, c(LOCAL_R_ROOT = root)), "overwrite")
missing_dest <- tempfile("archive-missing-")
assert_error(prepare_source_archive(fixture, missing_dest, c(LOCAL_OTHER_ROOT = root)), "Unconfigured")
stopifnot(!dir.exists(missing_dest))
assert_error(prepare_source_archive(fixture, tempfile(), c(LOCAL_R_ROOT = 'C:/bad"path')), "absolute paths")

r_files <- list.files("analysis", pattern = "[.][Rr]$", recursive = TRUE, full.names = TRUE)
r_files <- r_files[!grepl("/local_sources/", r_files, fixed = TRUE)]
for (f in r_files) parse(f, encoding = "UTF-8")
roots <- setNames(rep(root, 5), c("LOCAL_R_ROOT", "LOCAL_CLUSTER_ROOT", "LOCAL_PROJECT_ROOT",
                                "LOCAL_LEGACY_PROJECT_ROOT", "LOCAL_LEGACY_AUX_ROOT"))
all_prepared <- prepare_source_archive("analysis", tempfile("all-source-output-"), roots)
stopifnot(length(all_prepared) == 59L)
for (f in all_prepared) {
  stopifnot(!any(grepl("<LOCAL_[A-Z_]+>", readLines(f, warn = FALSE))))
}
cat("PASS: six age groups, 53-state aggregation, archive preparation, overwrite guards.\n")
cat("PASS: parsed", length(r_files), "R scripts without executing analyses.\n")
