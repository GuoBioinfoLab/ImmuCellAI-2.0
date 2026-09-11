# Exact interval boundaries retained from the original Figure 5 analysis.
figure5_age_levels <- c("0-1", "10-20", "20-30", "30-50", "50-70", "70+")
figure5_age_group <- function(age) {
  age <- as.numeric(age)
  out <- rep(NA_character_, length(age))
  valid <- is.finite(age)
  out[valid & age >= 0 & age < 2] <- "0-1"
  out[valid & age > 10 & age <= 20] <- "10-20"
  out[valid & age > 20 & age <= 30] <- "20-30"
  out[valid & age > 30 & age < 50] <- "30-50"
  out[valid & age >= 50 & age < 70] <- "50-70"
  out[valid & age >= 70] <- "70+"
  factor(out, levels = figure5_age_levels)
}

# These eight display groups are not the computational T-cell hierarchy.
figure5_major_cells <- list(
  Tcell = c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Tn", "CD4Trm",
            "CD8Tcm", "CD8Tem", "CD8Temra", "CD8Tn", "CD8Trm",
            "MAIT", "NKT", "Tc", "Tex", "Tfh", "Th1", "Th1/17",
            "Th17", "Th2", "Tr1", "Treg", "gdT"),
  Bcell = c("Bnaive", "Breg", "FOB", "BGC", "MZB", "MBC", "PC", "PB", "Bex"),
  ILC = c("ILC1", "ILC2", "ILC3"),
  NK = c("NKreg", "cNK"),
  Granulocyte = c("Neutrophil", "Basophil", "Eosinophil", "Mast cell"),
  Monocyte = c("cMo", "intMo", "ncMo", "MDSC"),
  Macrophage = c("M0", "M1", "M2", "TAM"),
  DC = c("pDC", "cDC1", "cDC2", "monoDC", "Langerhans")
)
