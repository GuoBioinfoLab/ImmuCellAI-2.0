# R API

## Recommended entry point

~~~r
run_immucellai2(
  bulk,
  reference = NULL,
  marker.genes = NULL,
  hierarchy.mode = "tcell",
  add.unknown = FALSE,
  n.cores = 1,
  pseudo.depth = 1e5,
  n.iter = 50,
  vb.tol = 1e-6,
  alpha.major = 10,
  alpha.sub = 5,
  alpha.state = 1,
  seed = 123,
  verbose = FALSE
)
~~~

### Required argument

- bulk: a genes-by-samples numeric matrix or a path to a tab-delimited file.
  Matrix row names must be gene symbols and column names must be unique sample
  identifiers.

### Reference and marker arguments

- reference: optional genes-by-cell-states matrix or file path. When omitted,
  the packaged 53-state reference atlas is used.
- marker.genes: optional character vector or one-gene-per-line text file.
  When omitted, the packaged 5,510-gene panel is used.

### Inference arguments

- hierarchy.mode: tcell, flat, hierarchical, or hybrid.
  tcell is the principal manuscript setting.
  The legacy value tcell_only remains accepted for backward compatibility.
- add.unknown: add a residual latent state. The principal analyses use FALSE.
- pseudo.depth: fixed pseudo-count mass used after within-sample scaling.
- n.iter: maximum number of deterministic posterior-mean updates.
- vb.tol: stop when the largest absolute state-fraction change is smaller
  than this value.
- alpha.major, alpha.sub, alpha.state: Dirichlet prior concentration
  parameters at the three hierarchy levels.
- n.cores: number of samples processed concurrently.
- seed: random seed used for reproducible worker initialization. The VB
  updates themselves are deterministic.

## Packaged resources

~~~r
reference <- load_immucellai2_reference()
markers <- load_immucellai2_markers()
~~~

reference is a genes-by-53-states matrix. markers is a vector of 5,510
symbols. The genes actually used for a dataset are the intersection of bulk
genes, reference genes, and marker genes.

## Reading input

~~~r
bulk <- read_expression_matrix("bulk_TPM.txt")
~~~

The reader:

1. treats the first column as gene identifiers;
2. requires numeric sample columns;
3. rejects negative values;
4. replaces missing or non-finite values by zero with a warning;
5. collapses duplicate genes by arithmetic mean;
6. preserves sample names exactly.

## Result object

~~~r
fit <- run_immucellai2(bulk)
names(fit)
~~~

Important fields are:

- state.fraction: samples by 53 terminal immune states;
- subtype.fraction: samples by hierarchy subtype;
- major.fraction: samples by major branch;
- hierarchy: state-to-subtype-to-major mapping;
- genes.used: exact gene vector supplied to inference;
- settings: all wrapper settings and marker counts;
- phi.final.list: reference probability matrix returned for each sample.

The primary result is state.fraction. With add.unknown = FALSE, each row
sums to one across the 53 modeled immune states.

## Writing results

~~~r
write_deconvolution_outputs(fit, "ImmuCellAI2_results")
~~~

Save the additional reproducibility fields explicitly:

~~~r
saveRDS(fit$settings, "ImmuCellAI2_results/settings.rds")
writeLines(fit$genes.used, "ImmuCellAI2_results/genes_used.txt")
capture.output(sessionInfo(), file = "ImmuCellAI2_results/sessionInfo.txt")
~~~

## Custom reference atlas

~~~r
fit <- run_immucellai2(
  bulk = bulk,
  reference = "custom_reference_TPM.txt",
  marker.genes = "custom_markers.txt",
  hierarchy.mode = "flat"
)
~~~

Custom references require independent validation. If custom state names are
not among the packaged 53 labels, use lower-level functions and provide an
explicit hierarchy table with columns state_name, subtype, and major_lineage.
