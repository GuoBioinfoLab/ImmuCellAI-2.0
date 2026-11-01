# ImmuCellAI2.0 0.1.9

- Set the default convergence tolerance to `1e-4` consistently in the standard
  entry point, matrix interface, and flat and hierarchical VB inference functions.
- Update the standard example, configurable analysis wrapper, and documentation
  to use the same tolerance.
- Retain the standard workflow's maximum of 50 iterations and support for explicit
  tolerance overrides. Archived analyses should use their recorded settings for
  exact reproduction.
- Add regression tests for default consistency, tolerance propagation, stopping
  behavior, explicit overrides, and iteration limits.
