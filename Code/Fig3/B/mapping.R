# Exact broad-group mappings from the historical external-cohort runner.
external_state_map <- list(
  Neutrophils = "Neutrophil",
  Monocytes = c("cMo", "intMo", "ncMo"),
  CD8T = c("CD8Tcm", "CD8Tem", "CD8Temra", "CD8Tn", "CD8Trm", "Tc", "Tex"),
  CD4T = c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Tn", "CD4Trm",
           "Tfh", "Th1", "Th1/17", "Th17", "Th2", "Tr1", "Treg"),
  B = c("BGC", "Bex", "Bnaive", "Breg", "FOB", "MBC", "MZB", "PB", "PC"),
  NK = c("cNK", "NKreg")
)
external_citmic_map <- list(
  Neutrophils = "Neutrophils",
  Monocytes = "Monocytes",
  CD8T = "CD8+ T cells",
  CD4T = "CD4+ T cells",
  B = "B cells",
  NK = "NK cells"
)
external_immucellai_map <- list(
  Neutrophils = "Neutrophil",
  Monocytes = "Monocyte",
  CD8T = "CD8_T",
  CD4T = "CD4_T",
  B = "Bcell",
  NK = "NK"
)
external_lm22_map <- list(
  Neutrophils = "Neutrophils",
  Monocytes = "Monocytes",
  CD8T = "T cells CD8",
  CD4T = c(
    "T cells CD4 naive",
    "T cells CD4 memory resting",
    "T cells CD4 memory activated",
    "T cells follicular helper",
    "T cells regulatory (Tregs)"
  ),
  B = c("B cells naive", "B cells memory", "Plasma cells"),
  NK = c("NK cells resting", "NK cells activated")
)

state_target_mapper <- function(target) {
  out <- lapply(as.character(target), function(x) external_state_map[[x]])
  names(out) <- target
  out
}
citmic_target_mapper <- function(target) {
  out <- lapply(as.character(target), function(x) external_citmic_map[[x]])
  names(out) <- target
  out
}
immucellai_target_mapper <- function(target) {
  out <- lapply(as.character(target), function(x) external_immucellai_map[[x]])
  names(out) <- target
  out
}
lm22_target_mapper <- function(target) {
  out <- lapply(as.character(target), function(x) external_lm22_map[[x]])
  names(out) <- target
  out
}
