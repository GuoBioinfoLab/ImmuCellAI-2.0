\name{run_immucellai2_custom}
\alias{run_immucellai2_custom}
\title{Deconvolve with a user-supplied reference}
\usage{
run_immucellai2_custom(
  bulk,
  reference,
  n.cores = 1,
  pseudo.depth = 1e+05,
  n.iter = 50,
  vb.tol = 1e-04,
  alpha = 1,
  seed = 123,
  verbose = FALSE
)
}
\arguments{
\item{bulk}{Genes-by-samples non-negative numeric matrix or tab-delimited file.}
\item{reference}{Genes-by-populations non-negative numeric matrix or
tab-delimited file. Gene identifiers must match those in bulk. TPM on a
compatible scale is recommended for both inputs.}
\item{n.cores}{Number of samples processed concurrently.}
\item{pseudo.depth}{Pseudo-count depth after within-sample scaling.}
\item{n.iter}{Maximum number of deterministic posterior-mean updates.}
\item{vb.tol}{Convergence tolerance for the maximum fraction change.}
\item{alpha}{Total concentration of the symmetric Dirichlet prior.}
\item{seed}{Seed for worker initialization; the inference is deterministic.}
\item{verbose}{Print inference progress.}
}
\value{
A list with subtype.fraction (samples by supplied reference columns),
genes.used, inference.method, phi.final.list, and settings. No type or
lineage aggregations are produced. Zero-expression samples and reference
columns on retained genes are rejected rather than silently omitted.
}
\description{
Uses all shared, informative genes and the deterministic flat Bayesian
posterior-mean engine. Reference columns are treated as separate populations
under a symmetric Dirichlet prior. No marker panel, biological hierarchy, or unknown
component is used. Arbitrary reference column names are supported.
}
