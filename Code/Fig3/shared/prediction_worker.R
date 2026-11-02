# Isolated child process; never source this worker from an interactive session.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
job <- readRDS(args[1]); task <- args[2]
stopifnot(file.exists(job$engine), all(file.exists(job$source_files)),
          identical(job$source_hashes, unname(tools::md5sum(job$source_files))))
root <- normalizePath(job$output, winslash = "/", mustWork = TRUE)
allowed <- paste0(normalizePath(file.path(job$code_dir, "output_analysis"), winslash = "/"), "/")
stopifnot(startsWith(tolower(paste0(root, "/")), tolower(allowed)))
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
.libPaths(c(file.path(job$library_root, c("Rlib_GSE176078_comparison", "Rlib_tools", "Rlib_test")), .libPaths()))
Sys.setenv(R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep),
           LANG = "en_US.UTF-8", LC_ALL = "English_United States.utf8")
source(job$engine, local = .GlobalEnv, encoding = "UTF-8")
