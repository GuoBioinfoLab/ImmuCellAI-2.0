# TCGA immune-phenotype analysis

> The complete executable workflow has moved to
> [`analysis/tumor_immune_phenotypes_immuicbscore/`](../tumor_immune_phenotypes_immuicbscore/README.md).
> This page is retained as a concise record of the clustering rationale.

The TCGA application uses the samples from cancer types selected for established
relevance to immune checkpoint blockade. ImmuCellAI 2.0 state fractions are
clustered across samples to identify four descriptive immune phenotypes.

## Required analysis record

Retain the following files together:

- sample-by-53-state fraction matrix;
- included TCGA cancer-type list and sample identifiers;
- clinical sample mapping;
- ConsensusClusterPlus settings and random seed;
- consensus class file for K = 4;
- survival endpoint definition;
- scripts used for the heatmap, Kaplan-Meier plot, and score association.

## Consensus clustering settings

The analysis used:

~~~r
ConsensusClusterPlus(
  d,
  distance = "euclidean",
  maxK = 6,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  innerLinkage = "complete",
  clusterAlg = "pam",
  corUse = "complete.obs",
  seed = 123,
  plot = "pdf"
)
~~~

Before clustering, each immune-state row is standardized across samples for
heatmap display. The unscaled compositional fractions are retained for
statistical comparisons and score calculation.

## Interpretation

The four immune phenotypes are data-driven sample clusters, not cell lineages
or temporal states. Their labels are assigned after inspection of dominant
immune composition patterns. The K = 4 choice should be supported by consensus
CDF/delta-area behavior, cluster stability, sample size, and interpretability,
rather than survival separation alone.
