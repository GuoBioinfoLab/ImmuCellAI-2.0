read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                                     fileEncoding = "UTF-8-BOM", quote = "\"", comment.char = "")
table_key <- function(x, cols) do.call(paste, c(x[cols], sep = "|"))
check_num <- function(a, b, label = "Numeric comparison", tolerance = 1e-12) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b) || !identical(is.na(a), is.na(b))) stop(label, ": missingness or length differs")
  finite <- is.finite(a) & is.finite(b)
  if (any(is.infinite(a) != is.infinite(b)) || any(a[!finite & !is.na(a)] != b[!finite & !is.na(b)]))
    stop(label, ": nonfinite values differ")
  if (any(abs(a[finite] - b[finite]) > tolerance)) stop(label, ": max difference ", max(abs(a[finite] - b[finite])))
  invisible(TRUE)
}
check_table <- function(a, b, keys, cols = intersect(names(a), names(b)), tolerance = 1e-12) {
  ka <- table_key(a, keys); kb <- table_key(b, keys)
  stopifnot(!anyDuplicated(ka), !anyDuplicated(kb), setequal(ka, kb))
  b <- b[match(ka, kb), , drop = FALSE]
  for (column in cols) {
    if (is.numeric(a[[column]]) && is.numeric(b[[column]])) check_num(a[[column]], b[[column]], column, tolerance)
    else stopifnot(identical(as.character(a[[column]]), as.character(b[[column]])))
  }
  invisible(TRUE)
}
auc_rank <- function(y, score) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  stopifnot(n1 > 0, n0 > 0, all(is.finite(score)), all(y %in% c(0, 1)))
  (sum(rank(score, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

prepare_analysis <- function(data_dir) {
  get <- function(path) read_tsv(file.path(data_dir, path))
  patients <- get("shared/TCGA_patients.tsv")
  complete <- with(patients, !is.na(age) & age > 0 & age < 120 & sex %in% c("female", "male") & stage %in% c("I", "II", "III", "IV"))
  stopifnot(nrow(patients) == 4746L, !anyDuplicated(patients$patient), length(unique(patients$cancer)) == 14L,
            identical(complete, tolower(as.character(patients$clinical_complete)) == "true"))
  check_table(patients[complete, ], get("E/clinical_complete_patients.tsv"), "patient")
  list(patients = patients, complete = patients[complete, ], verification = list(patients_D = 4746L, patients_E = 4593L, clustering_rerun = FALSE))
}
