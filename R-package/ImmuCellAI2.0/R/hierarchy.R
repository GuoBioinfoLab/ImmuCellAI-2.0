#' Assign broad immune lineage labels for 53-state reference names
#'
#' @param state.labels Character vector of state names.
#' @return data.frame with state_name, subtype, major_lineage.
base_53_hierarchy <- function(state.labels) {
  state.labels <- as.character(state.labels)
  hierarchy <- data.frame(
    state_name = state.labels,
    subtype = NA_character_,
    major_lineage = NA_character_,
    stringsAsFactors = FALSE
  )

  b.cells <- c("BGC", "Bex", "Bnaive", "Breg", "FOB", "MBC", "MZB", "PB", "PC")
  t.cells <- c(
    "CD4Tcm", "CD4Tem", "CD4Temra", "CD4Tn", "CD4Tnaive", "CD4Trm",
    "CD8Tcm", "CD8Tem", "CD8Temra", "CD8Tn", "CD8Tnaive", "CD8Trm",
    "MAIT", "NKT", "Tc", "Tex", "exhausted_T", "Tfh",
    "Th1", "Th1/17", "Th1/Th17", "Th1.Th17", "Th1_Th17", "Th17", "Th2",
    "Tr1", "Treg", "gdT"
  )
  nk.cells <- c("NKreg", "cNK")
  ilc.cells <- c("ILC1", "ILC2", "ILC3")
  granulocytes <- c("Basophil", "Eosinophil", "Neutrophil", "Mast cell")
  macrophages <- c("M0", "M1", "M2", "TAM")
  dc.cells <- c("Langerhans", "cDC1", "cDC2", "monoDC", "pDC")
  monocytes <- c("cMo", "CMonocyte", "intMo", "iMo", "IMonocyte", "ncMo", "NMonocyte")
  mdsc.cells <- c("MDSC")
  unknown.cells <- c("UNKNOWN", "Unknown", "unknown")

  hierarchy$subtype[hierarchy$state_name %in% b.cells] <- "B_cell"
  hierarchy$major_lineage[hierarchy$state_name %in% b.cells] <- "Lymphoid"

  hierarchy$subtype[hierarchy$state_name %in% t.cells] <- "T_cell"
  hierarchy$major_lineage[hierarchy$state_name %in% t.cells] <- "Lymphoid"

  hierarchy$subtype[hierarchy$state_name %in% nk.cells] <- "NK"
  hierarchy$major_lineage[hierarchy$state_name %in% nk.cells] <- "Lymphoid"

  hierarchy$subtype[hierarchy$state_name %in% ilc.cells] <- "ILC"
  hierarchy$major_lineage[hierarchy$state_name %in% ilc.cells] <- "Lymphoid"

  hierarchy$subtype[hierarchy$state_name %in% granulocytes] <- "Granulocyte"
  hierarchy$major_lineage[hierarchy$state_name %in% granulocytes] <- "Myeloid"

  hierarchy$subtype[hierarchy$state_name %in% macrophages] <- "Macrophage"
  hierarchy$major_lineage[hierarchy$state_name %in% macrophages] <- "Myeloid"

  hierarchy$subtype[hierarchy$state_name %in% dc.cells] <- "DC"
  hierarchy$major_lineage[hierarchy$state_name %in% dc.cells] <- "Myeloid"

  hierarchy$subtype[hierarchy$state_name %in% monocytes] <- "Monocyte"
  hierarchy$major_lineage[hierarchy$state_name %in% monocytes] <- "Myeloid"

  hierarchy$subtype[hierarchy$state_name %in% mdsc.cells] <- "MDSC"
  hierarchy$major_lineage[hierarchy$state_name %in% mdsc.cells] <- "Myeloid"

  hierarchy$subtype[hierarchy$state_name %in% unknown.cells] <- "UNKNOWN"
  hierarchy$major_lineage[hierarchy$state_name %in% unknown.cells] <- "UNKNOWN"

  unassigned <- is.na(hierarchy$subtype) | is.na(hierarchy$major_lineage)
  if (any(unassigned)) {
    warning(
      "Some states were not assigned by the default hierarchy: ",
      paste(hierarchy$state_name[unassigned], collapse = ", "),
      ". They will be assigned to subtype = state_name and major_lineage = Other."
    )
    hierarchy$subtype[unassigned] <- hierarchy$state_name[unassigned]
    hierarchy$major_lineage[unassigned] <- "Other"
  }

  rownames(hierarchy) <- NULL
  hierarchy
}

#' Create a flat hierarchy
#'
#' In this mode every state is its own subtype. The major lineage labels are
#' still retained. This is the recommended control/default mode for Experiment 1
#' when testing whether hierarchy itself improves performance.
#'
#' @param state.labels Character vector of state names.
#' @return data.frame with state_name, subtype, major_lineage.
#' @export
create_flat_53_hierarchy <- function(state.labels) {
  hierarchy <- base_53_hierarchy(state.labels)
  hierarchy$subtype <- hierarchy$state_name

  unknown.idx <- hierarchy$state_name %in% c("UNKNOWN", "Unknown", "unknown")
  hierarchy$subtype[unknown.idx] <- "UNKNOWN"

  ## True flat control: all known immune states compete in one top-level pool.
  ## UNKNOWN remains an independent major lineage when present.
  hierarchy$major_lineage[!unknown.idx] <- "Immune"
  hierarchy$major_lineage[unknown.idx] <- "UNKNOWN"

  rownames(hierarchy) <- NULL
  hierarchy
}

#' Refine T-cell subtype labels in a hierarchy table
#'
#' This creates the refined T-cell hierarchy used for the strict hierarchical
#' benchmark: CD4/CD8 memory and naive groups are separate; helper, Treg, Tr1,
#' Tfh, Tc, Tex, MAIT, NKT and gdT are separate subtypes.
#'
#' @param hierarchy data.frame with state_name, subtype, major_lineage.
#' @return Modified hierarchy.
#' @export
refine_tcell_subtypes <- function(hierarchy) {
  hierarchy <- as.data.frame(hierarchy, stringsAsFactors = FALSE)
  required <- c("state_name", "subtype", "major_lineage")
  miss <- setdiff(required, colnames(hierarchy))
  if (length(miss) > 0) stop("hierarchy is missing columns: ", paste(miss, collapse = ", "))

  hierarchy$subtype[hierarchy$state_name %in% c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Trm")] <- "CD4 Memory T cell"
  hierarchy$subtype[hierarchy$state_name %in% c("CD8Tcm", "CD8Tem", "CD8Temra", "CD8Trm")] <- "CD8 Memory T cell"
  hierarchy$subtype[hierarchy$state_name %in% c("CD4Tn", "CD4Tnaive")] <- "CD4 T naive"
  hierarchy$subtype[hierarchy$state_name %in% c("CD8Tn", "CD8Tnaive")] <- "CD8 T naive"
  hierarchy$subtype[hierarchy$state_name %in% c("Treg")] <- "Treg"
  hierarchy$subtype[hierarchy$state_name %in% c("Tr1")] <- "Tr1"
  hierarchy$subtype[hierarchy$state_name %in% c("Th1", "Th2", "Th17", "Th1/17", "Th1/Th17", "Th1.Th17", "Th1_Th17")] <- "Helper T"
  hierarchy$subtype[hierarchy$state_name %in% c("Tfh")] <- "Tfh"
  hierarchy$subtype[hierarchy$state_name %in% c("Tc")] <- "Tc"
  hierarchy$subtype[hierarchy$state_name %in% c("Tex", "exhausted_T")] <- "Tex"
  hierarchy$subtype[hierarchy$state_name %in% c("MAIT")] <- "MAIT"
  hierarchy$subtype[hierarchy$state_name %in% c("NKT")] <- "NKT"
  hierarchy$subtype[hierarchy$state_name %in% c("gdT")] <- "gdT"

  refined.t.subtypes <- c(
    "CD4 Memory T cell", "CD8 Memory T cell", "CD4 T naive", "CD8 T naive",
    "Treg", "Tr1", "Helper T", "Tfh", "Tc", "Tex", "MAIT", "NKT", "gdT"
  )
  hierarchy$major_lineage[hierarchy$subtype %in% refined.t.subtypes] <- "Lymphoid"
  rownames(hierarchy) <- NULL
  hierarchy
}

#' Refine monocyte subtype labels in a hierarchy table
#'
#' This splits cMo/intMo/ncMo into independent subtypes under Myeloid.
#'
#' @param hierarchy data.frame with state_name, subtype, major_lineage.
#' @return Modified hierarchy.
#' @export
refine_monocyte_subtypes <- function(hierarchy) {
  hierarchy <- as.data.frame(hierarchy, stringsAsFactors = FALSE)
  required <- c("state_name", "subtype", "major_lineage")
  miss <- setdiff(required, colnames(hierarchy))
  if (length(miss) > 0) stop("hierarchy is missing columns: ", paste(miss, collapse = ", "))

  hierarchy$subtype[hierarchy$state_name %in% c("cMo", "CMonocyte")] <- "cMo"
  hierarchy$subtype[hierarchy$state_name %in% c("intMo", "iMo", "IMonocyte")] <- "intMo"
  hierarchy$subtype[hierarchy$state_name %in% c("ncMo", "NMonocyte")] <- "ncMo"
  hierarchy$major_lineage[hierarchy$subtype %in% c("cMo", "intMo", "ncMo")] <- "Myeloid"
  rownames(hierarchy) <- NULL
  hierarchy
}

#' Create a hybrid 53-state hierarchy
#'
#' Hybrid mode reduces strong hierarchy where it previously hurt performance:
#' CD8 memory/naive remains grouped, CD4 memory is grouped with helper states,
#' and monocytes are split into independent cMo/intMo/ncMo subtypes.
#'
#' @param state.labels Character vector of state names.
#' @return data.frame with state_name, subtype, major_lineage.
#' @export
create_hybrid_53_hierarchy <- function(state.labels) {
  hierarchy <- base_53_hierarchy(state.labels)

  hierarchy$subtype[hierarchy$state_name %in% c("CD4Tn", "CD4Tnaive")] <- "CD4 T naive"
  hierarchy$subtype[hierarchy$state_name %in% c(
    "CD4Tcm", "CD4Tem", "CD4Temra", "CD4Trm",
    "Th1", "Th2", "Th17", "Th1/17", "Th1/Th17", "Th1.Th17", "Th1_Th17"
  )] <- "CD4 Memory/Helper T"
  hierarchy$subtype[hierarchy$state_name %in% c("CD8Tn", "CD8Tnaive")] <- "CD8 T naive"
  hierarchy$subtype[hierarchy$state_name %in% c("CD8Tcm", "CD8Tem", "CD8Temra", "CD8Trm")] <- "CD8 Memory T cell"
  hierarchy$subtype[hierarchy$state_name %in% c("Tfh")] <- "Tfh"
  hierarchy$subtype[hierarchy$state_name %in% c("Treg")] <- "Treg"
  hierarchy$subtype[hierarchy$state_name %in% c("Tr1")] <- "Tr1"
  hierarchy$subtype[hierarchy$state_name %in% c("Tc")] <- "Tc"
  hierarchy$subtype[hierarchy$state_name %in% c("Tex", "exhausted_T")] <- "Tex"
  hierarchy$subtype[hierarchy$state_name %in% c("MAIT")] <- "MAIT"
  hierarchy$subtype[hierarchy$state_name %in% c("NKT")] <- "NKT"
  hierarchy$subtype[hierarchy$state_name %in% c("gdT")] <- "gdT"
  hierarchy <- refine_monocyte_subtypes(hierarchy)

  lymphoid.subtypes <- c(
    "CD4 T naive", "CD4 Memory/Helper T", "CD8 T naive", "CD8 Memory T cell",
    "Tfh", "Treg", "Tr1", "Tc", "Tex", "MAIT", "NKT", "gdT"
  )
  hierarchy$major_lineage[hierarchy$subtype %in% lymphoid.subtypes] <- "Lymphoid"
  rownames(hierarchy) <- NULL
  hierarchy
}

#' Create a T-cell-only refined hierarchy
#'
#' Only CD4/CD8 T-cell related states are modeled with an internal hierarchy.
#' All other immune states are left as independent top-level branches, avoiding
#' shrinkage across unrelated immune lineages.
#'
#' @param state.labels Cell-type names.
#' @return data.frame with state_name, subtype, and major_lineage.
#' @export
create_tcell_only_53_hierarchy <- function(state.labels) {
  hierarchy <- data.frame(
    state_name = as.character(state.labels),
    subtype = as.character(state.labels),
    major_lineage = as.character(state.labels),
    stringsAsFactors = FALSE
  )

  cd4.memory <- c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Trm")
  cd4.naive <- c("CD4Tn", "CD4Tnaive")
  cd4.helper <- c("Th1", "Th1/17", "Th1/Th17", "Th1.Th17", "Th1_Th17", "Th17", "Th2")
  cd4.special <- c("Tfh", "CD4Tfh", "Treg", "Tr1")

  cd8.memory <- c("CD8Tcm", "CD8Tem", "CD8Temra", "CD8Trm")
  cd8.naive <- c("CD8Tn", "CD8Tnaive")
  cd8.effector <- c("Tc", "Tex", "exhausted_T")

  hierarchy$major_lineage[hierarchy$state_name %in% c(cd4.memory, cd4.naive, cd4.helper, cd4.special)] <- "CD4T"
  hierarchy$subtype[hierarchy$state_name %in% cd4.memory] <- "CD4 Memory T"
  hierarchy$subtype[hierarchy$state_name %in% cd4.naive] <- "CD4 Naive T"
  hierarchy$subtype[hierarchy$state_name %in% cd4.helper] <- "CD4 Helper T"
  hierarchy$subtype[hierarchy$state_name %in% c("Tfh", "CD4Tfh")] <- "Tfh"
  hierarchy$subtype[hierarchy$state_name %in% "Treg"] <- "Treg"
  hierarchy$subtype[hierarchy$state_name %in% "Tr1"] <- "Tr1"

  hierarchy$major_lineage[hierarchy$state_name %in% c(cd8.memory, cd8.naive, cd8.effector)] <- "CD8T"
  hierarchy$subtype[hierarchy$state_name %in% cd8.memory] <- "CD8 Memory T"
  hierarchy$subtype[hierarchy$state_name %in% cd8.naive] <- "CD8 Naive T"
  hierarchy$subtype[hierarchy$state_name %in% "Tc"] <- "Tc"
  hierarchy$subtype[hierarchy$state_name %in% c("Tex", "exhausted_T")] <- "Tex"

  rownames(hierarchy) <- NULL
  hierarchy
}

#' Create default 53-state hierarchy
#'
#' @param state.labels Cell-type names, usually colnames(reference).
#' @param hierarchy.mode One of flat, hierarchical, or hybrid.
#' @param refine.tcell.subtypes Backward-compatible option. Used only when hierarchy.mode = "hierarchical".
#' @return data.frame with state_name, subtype, and major_lineage.
#' @export
create_default_53_hierarchy <- function(
  state.labels,
  hierarchy.mode = c("flat", "hierarchical", "hybrid", "tcell_only"),
  refine.tcell.subtypes = TRUE
) {
  hierarchy.mode <- match.arg(hierarchy.mode)

  if (hierarchy.mode == "flat") {
    return(create_flat_53_hierarchy(state.labels))
  }

  if (hierarchy.mode == "hybrid") {
    return(create_hybrid_53_hierarchy(state.labels))
  }

  if (hierarchy.mode == "tcell_only") {
    return(create_tcell_only_53_hierarchy(state.labels))
  }

  hierarchy <- base_53_hierarchy(state.labels)
  if (isTRUE(refine.tcell.subtypes)) {
    hierarchy <- refine_tcell_subtypes(hierarchy)
  }
  rownames(hierarchy) <- NULL
  hierarchy
}

validate_hierarchy <- function(hierarchy, reference) {
  required <- c("state_name", "subtype", "major_lineage")
  miss <- setdiff(required, colnames(hierarchy))
  if (length(miss) > 0) stop("hierarchy is missing columns: ", paste(miss, collapse = ", "))
  if (!all(hierarchy$state_name %in% colnames(reference))) {
    stop("Some hierarchy$state_name values are not in colnames(reference): ",
         paste(setdiff(hierarchy$state_name, colnames(reference)), collapse = ", "))
  }
  hierarchy <- hierarchy[match(colnames(reference), hierarchy$state_name), , drop = FALSE]
  rownames(hierarchy) <- NULL
  hierarchy
}

#' Aggregate state fractions by hierarchy
#'
#' @param state.mat Samples x states matrix.
#' @param hierarchy Hierarchy table.
#' @return list with subtype and major matrices.
#' @export
aggregate_by_hierarchy <- function(state.mat, hierarchy) {
  state.mat <- as.matrix(state.mat)
  common <- intersect(colnames(state.mat), hierarchy$state_name)
  h <- hierarchy[match(common, hierarchy$state_name), , drop = FALSE]
  x <- state.mat[, common, drop = FALSE]

  aggregate_cols <- function(mat, group) {
    group <- as.character(group)
    group.names <- sort(unique(group))
    out <- matrix(0, nrow = nrow(mat), ncol = length(group.names))
    rownames(out) <- rownames(mat)
    colnames(out) <- group.names
    for (g in group.names) {
      out[, g] <- rowSums(mat[, group == g, drop = FALSE])
    }
    out
  }

  subtype <- aggregate_cols(x, h$subtype)
  major <- aggregate_cols(x, h$major_lineage)
  list(subtype = subtype, major = major)
}
