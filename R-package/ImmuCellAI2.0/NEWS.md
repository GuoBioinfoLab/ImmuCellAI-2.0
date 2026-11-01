# ImmuCellAI2.0 0.2.0

- Rename the three fraction outputs to match lineage -> cell type -> subtype:
  `state.fraction` becomes `subtype.fraction`, the old `subtype.fraction`
  becomes `type.fraction`, and `major.fraction` becomes `lineage.fraction`.
- Export `subtype_fraction.txt`, `type_fraction.txt`, and `lineage_fraction.txt`.
  All fraction values, hierarchy assignments, references, and inference settings
  are unchanged. Hierarchy-table columns keep their legacy names for custom
  input compatibility: state_name (terminal subtype), subtype (parent cell type),
  major_lineage (lineage).
- This is a breaking output-schema change. No old fraction aliases are returned.
  The writer rejects legacy result fields and legacy output directories to avoid
  silently mixing the two meanings of subtype_fraction.txt.
- To migrate a saved legacy fit without recomputing, rename all three fields
  simultaneously in R:

  ```r
  old_names <- c("state.fraction", "subtype.fraction", "major.fraction")
  stopifnot(all(old_names %in% names(fit)))
  names(fit)[match(old_names, names(fit))] <-
    c("subtype.fraction", "type.fraction", "lineage.fraction")
  ```

  Export into a new directory after migration. Updated examples require v0.2.0.
- Update current consumers and add regression tests for all output levels,
  aggregation, UNKNOWN, file round-trips, and legacy export protection.

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
