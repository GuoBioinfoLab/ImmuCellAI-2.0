# R API

After installing the package, use the namespace-qualified calls below without
attaching the package with `library()`.

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
  the packaged 53-subtype reference atlas is used. Supplying a reference selects
  the simplified custom-reference workflow described below.
- marker.genes: optional character vector or one-gene-per-line text file.
  For the bundled reference, the packaged 5,510-gene panel is used when omitted.
  This option is not used with a custom reference.

### Inference arguments

- hierarchy.mode: `tcell` or `flat`.
  tcell is the principal manuscript setting.
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

For the bundled reference, important fields are:

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
columns follow the selected hierarchy. See the
[input/output specification](INPUT_OUTPUT.md) for details.

## Writing results

~~~r
ImmuCellAI2.0::write_deconvolution_outputs(fit, "ImmuCellAI2_results")
~~~

This writes `subtype_fraction.txt`, `type_fraction.txt`,
`lineage_fraction.txt`, and `hierarchy_used.txt`.

Save the additional reproducibility fields explicitly:

~~~r
saveRDS(fit$settings, "ImmuCellAI2_results/settings.rds")
writeLines(fit$genes.used, "ImmuCellAI2_results/genes_used.txt")
capture.output(sessionInfo(), file = "ImmuCellAI2_results/sessionInfo.txt")
~~~

## Custom reference atlas

For users wishing to "retrain" the method for a new dataset or species, such as
livestock immunity studies, prepare a species-matched reference expression
matrix with one column per cell population and a bulk expression matrix with
matching gene identifiers. Use non-negative, non-log-transformed TPM values and
supply both files to `run_immucellai2_custom()` as shown below.

The repository includes a [bulk TPM matrix](../example_data/custom_reference/bulk_TPM.txt)
and a [reference TPM matrix](../example_data/custom_reference/custom_reference_TPM.txt).
This is a synthetic input-format demonstration with 1,500 genes, three reference
populations, and four bulk samples, not a biological validation dataset. Each
expression column sums to one million. Run from the repository root:

~~~r
fit <- ImmuCellAI2.0::run_immucellai2_custom(
  bulk = "example_data/custom_reference/bulk_TPM.txt",
  reference = "example_data/custom_reference/custom_reference_TPM.txt"
)
head(fit$subtype.fraction)
ImmuCellAI2.0::write_deconvolution_outputs(fit, "custom_reference_results")
~~~

The files can be regenerated with `Rscript examples/00_make_custom_reference_demo.R`.
The supplied gene labels are synthetic; for real data, use matching gene
identifiers in the bulk and reference matrices.

Both inputs can be numeric matrices or tab-delimited files with gene identifiers in
the first column. Rows are genes; bulk columns are samples and reference columns
are the populations to estimate. Column names must be unique and non-empty.
Reference population names are not restricted to the packaged 53 labels. Use
non-negative expression on a compatible scale, preferably TPM, in both inputs.

All shared genes are used without a predefined marker panel. Genes with zero
total expression across all bulk samples or across all reference columns are
removed. Samples or reference profiles with no expression on the retained genes
are rejected. Duplicate gene symbols are collapsed by their mean.

The workflow uses the existing deterministic flat Bayesian posterior-mean
updates with a symmetric Dirichlet prior. It does not impose a biological
hierarchy, add an unknown component, or adapt reference profiles. No marker,
hierarchy, or UNKNOWN arguments are exposed by this entry point. Optional
numerical controls are `n.cores`, `pseudo.depth`, `n.iter`, `vb.tol`, `alpha`
(total prior concentration), `seed`, and `verbose`.

`fit$subtype.fraction` contains samples in rows and the supplied reference
populations in columns. The number and order of columns match the reference;
weights sum to one within each sample. The result includes `genes.used`,
`inference.method`, `phi.final.list`, and `settings`, but no `type.fraction`,
`lineage.fraction`, or hierarchy table. The output writer saves only
`subtype_fraction.txt` for this workflow.

`run_immucellai2(bulk, reference = reference)` dispatches to the same custom
workflow. Do not pass marker, hierarchy, or UNKNOWN options with a supplied
reference. Omitting `reference` retains the standard bundled-reference workflow.
Alternative references require validation for the intended application.
