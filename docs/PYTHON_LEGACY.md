# Legacy Python implementation

The `immucellai2/` directory contains the original Python Gibbs/MCMC code. It is
retained to avoid silently breaking earlier users, but it is not the
deterministic variational Bayesian implementation described in the current
ImmuCellAI 2.0 manuscript.

For manuscript reproduction and new analyses, use the R package under
`R-package/ImmuCellAI2.0/` and the `run_immucellai2()` entry point.

The legacy Python package also expects reference resources from its historical
distribution. Results from it should not be presented as output of the current
VB method unless the implementation and parameterization are stated explicitly.
