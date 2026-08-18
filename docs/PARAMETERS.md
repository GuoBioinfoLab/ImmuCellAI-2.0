# Parameter reference

| Parameter | Manuscript value | Meaning |
|---|---:|---|
| `hierarchy.mode` | `tcell` | Hierarchical refinement within CD4 and CD8 branches only |
| `inference.method` | `vb` | Deterministic posterior-mean iteration |
| `add.unknown` | `FALSE` | Optional residual state disabled |
| `pseudo.depth` | `100000` | Total pseudo-count mass assigned to each sample |
| `n.iter` | `50` | Maximum VB iterations |
| `vb.tol` | `1e-6` | Maximum absolute fraction change required for convergence |
| `alpha.major` | `10` | Major-branch Dirichlet prior concentration |
| `alpha.sub` | `5` | Subtype Dirichlet prior concentration |
| `alpha.state` | `1` | State-level Dirichlet prior concentration |
| `seed` | `123` | Reproducibility seed |
| `n.cores` | user selected | Number of samples processed concurrently |

`n.cores` changes throughput across samples but not the mathematical update for
one sample. The VB updates within an individual sample are sequential.

## Sensitivity modes

The four-mode comparison used during development is:

1. `tcell`, UNKNOWN off;
2. `tcell`, UNKNOWN on;
3. `flat`, UNKNOWN off;
4. `flat`, UNKNOWN on.

These alternatives should be described as sensitivity analyses. Selecting a
different mode separately for every test cohort can introduce evaluation bias.
