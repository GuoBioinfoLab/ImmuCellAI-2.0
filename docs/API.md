# R API

The examples use ImmuCellAI2.0 v0.2.0 or later. Namespace-qualified calls do not
require attaching the package with `library()`. After upgrading, restart R so
that an older loaded namespace is not reused.

## Recommended entry point

~~~r
ImmuCellAI2.0::run_immucellai2(
  bulk,
  reference = NULL,
  marker.genes = NULL,
  hierarchy.mode = "tcell",
  add.unknown = FALSE,
  n.cores = 1,
  pseudo.depth = 1e5,
  n.iter = 50,
  vb.tol = 1e-4,
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

- reference: optional genes-by-cell-subtypes matrix or file path. When omitted,
  the packaged 53-subtype reference atlas is used.
- marker.genes: optional character vector or one-gene-per-line text file.
  When omitted, the packaged 5,510-gene panel is used.

### Inference arguments

- hierarchy.mode: tcell, flat, hierarchical, or hybrid.
  tcell is the principal manuscript setting.
  The legacy value tcell_only remains accepted for backward compatibility.
- add.unknown: add a residual component. The principal analyses use FALSE.
- pseudo.depth: fixed pseudo-count mass used after within-sample scaling.
- n.iter: maximum number of deterministic posterior-mean updates.
- vb.tol: stop when the largest absolute terminal-subtype fraction change is smaller
  than this value.
- alpha.major, alpha.sub, alpha.state: Dirichlet prior concentration
  parameters at the three hierarchy levels.
- n.cores: number of samples processed concurrently.
- seed: random seed used for reproducible worker initialization. The VB
  updates themselves are deterministic.

## Packaged resources

~~~r
reference <- ImmuCellAI2.0::load_immucellai2_reference()
markers <- ImmuCellAI2.0::load_immucellai2_markers()
~~~

reference is a genes-by-53-subtypes matrix. markers is a vector of 5,510
symbols. The genes actually used for a dataset are the intersection of bulk
genes, reference genes, and marker genes.

## Reading input

~~~r
bulk <- ImmuCellAI2.0::read_expression_matrix("bulk_TPM.txt")
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
fit <- ImmuCellAI2.0::run_immucellai2(bulk)
names(fit)
~~~

Important fields are:

- subtype.fraction: samples by 53 terminal immune subtypes;
- type.fraction: samples by middle-level cell type;
- lineage.fraction: samples by top-level lineage;
- hierarchy: mapping for lineage -> cell type -> subtype;
- genes.used: exact gene vector supplied to inference;
- settings: all wrapper settings and marker counts;
- phi.final.list: reference probability matrix returned for each sample.

The primary result is subtype.fraction. With add.unknown = FALSE, each row
sums to one across the 53 modeled immune subtypes. Type and lineage matrices
contain aggregated weights, not conditional within-parent proportions. Their
columns follow the selected hierarchy, which is unchanged in v0.2.0.
For legacy objects and the unchanged hierarchy-table column names, see
[input/output migration notes](INPUT_OUTPUT.md#upgrading-from-v019-or-earlier).

## Writing results

~~~r
ImmuCellAI2.0::write_deconvolution_outputs(fit, "ImmuCellAI2_results")
~~~

This writes `subtype_fraction.txt`, `type_fraction.txt`,
`lineage_fraction.txt`, and `hierarchy_used.txt`. Do not reuse a directory
containing old `state_fraction.txt` or `major_lineage_fraction.txt` files.

Save the additional reproducibility fields explicitly:

~~~r
saveRDS(fit$settings, "ImmuCellAI2_results/settings.rds")
writeLines(fit$genes.used, "ImmuCellAI2_results/genes_used.txt")
capture.output(sessionInfo(), file = "ImmuCellAI2_results/sessionInfo.txt")
~~~

## Custom reference atlas

~~~r
fit <- ImmuCellAI2.0::run_immucellai2(
  bulk = bulk,
  reference = "custom_reference_TPM.txt",
  marker.genes = "custom_markers.txt",
  hierarchy.mode = "flat"
)
~~~

Custom references require independent validation. If custom subtype names are
not among the packaged 53 labels, use lower-level functions and provide an
explicit hierarchy table with columns state_name, subtype, and major_lineage.
