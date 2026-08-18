args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args)) args[1] else "analysis/tumor_immune_phenotypes_immuicbscore/config.R"
source(config_file)
source(file.path(repo_root, "analysis", "common", "analysis_utils.R"))
require_packages(c("ConsensusClusterPlus", "data.table"))

fraction <- read_fraction_matrix(tcga_fraction_file)
if (ncol(fraction) != 53L) warning("Expected 53 immune states; found ", ncol(fraction), ".")

# ConsensusClusterPlus expects features in rows and samples in columns.
d <- t(fraction)
out <- file.path(output_dir, "consensus_clustering")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(consensus_seed)
results <- ConsensusClusterPlus::ConsensusClusterPlus(
  d,
  distance = "euclidean",
  maxK = 6,
  reps = consensus_reps,
  pItem = 0.8,
  pFeature = 1,
  title = out,
  innerLinkage = "complete",
  finalLinkage = "complete",
  clusterAlg = "pam",
  corUse = "complete.obs",
  seed = consensus_seed,
  plot = "pdf",
  writeTable = TRUE
)
saveRDS(results, file.path(out, "consensus_results.rds"))
class_k4 <- data.frame(
  sample = names(results[[4L]]$consensusClass),
  cluster = paste0("IP", as.integer(results[[4L]]$consensusClass)),
  stringsAsFactors = FALSE
)
data.table::fwrite(class_k4, file.path(out, "consensus_class_k4.tsv"), sep = "\t")
