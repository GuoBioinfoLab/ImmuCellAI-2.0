# Exact operational maps extracted from the historical runner.
state_target_mapper <- function(target) {
  mapper <- list(
    Bcells = c("BGC", "Bex", "Bnaive", "Breg", "FOB", "MBC", "MZB", "PB", "PC"),
    CD4Tcm = "CD4Tcm",
    CD4Tem = "CD4Tem",
    CD4Temra = "CD4Temra",
    CD4Tfh = "Tfh",
    CD4Tnaive = "CD4Tn",
    CD4Trm = "CD4Trm",
    CD8Tcm = "CD8Tcm",
    CD8Tem = "CD8Tem",
    CD8Temra = "CD8Temra",
    CD8Tnaive = "CD8Tn",
    CD8Trm = "CD8Trm",
    MAIT = "MAIT",
    Monocyte = c("cMo", "intMo", "ncMo"),
    NK = c("cNK", "NKreg"),
    NKT = "NKT",
    Th1 = "Th1",
    Th17 = "Th17",
    monoDC = "monoDC",
    `γδT` = "gdT",
    CD4Tmemory = c("CD4Tcm", "CD4Tem", "CD4Temra", "CD4Trm"),
    CD4Tn = "CD4Tn",
    CD8Tmemory = c("CD8Tcm", "CD8Tem", "CD8Temra", "CD8Trm"),
    CD8Tn = "CD8Tn",
    CMonocyte = "cMo",
    FOB = "FOB",
    GC_B = "BGC",
    ILC3 = "ILC3",
    NMonocyte = "ncMo",
    TAM = "TAM",
    Thelper = c("Tfh", "Th1", "Th1/17", "Th17", "Th2"),
    Treg = "Treg",
    cDC1 = "cDC1",
    cDC2 = "cDC2",
    cytotoxicNK = "cNK",
    exhausted_T = "Tex",
    macrophage = c("M0", "M1", "M2"),
    mast_cell = "Mast cell",
    pDC = "pDC",
    plasma = c("PB", "PC"),
    regulatoryNK = "NKreg"
  )
  out <- lapply(as.character(target), function(x) if (x %in% names(mapper)) mapper[[x]] else x)
  names(out) <- target
  out
}

citmic_target_mapper <- function(target) {
  mapper <- list(
    Bcells = "B cells",
    CD4Tcm = "CD4+ Tcm",
    CD4Tem = "CD4+ Tem",
    CD4Temra = character(0),
    CD4Tfh = "Tfh",
    CD4Tnaive = "Naive CD4+ T cell",
    CD4Trm = character(0),
    CD8Tcm = "CD8+ Tcm",
    CD8Tem = "CD8+ Tem",
    CD8Temra = character(0),
    CD8Tnaive = "Naive CD8+ T cell",
    CD8Trm = character(0),
    MAIT = "MAIT cells",
    Monocyte = "Monocytes",
    NK = "NK cells",
    NKT = "NKT cells",
    Th1 = "Th1 cells",
    Th17 = "Th17 cells",
    monoDC = character(0),
    `γδT` = "Tgd cells",
    CD4Tmemory = c("CD4+ memory T cells", "CD4+ Tcm", "CD4+ Tem"),
    CD4Tn = "Naive CD4+ T cell",
    CD8Tmemory = c("CD8+ Tcm", "CD8+ Tem"),
    CD8Tn = "Naive CD8+ T cell",
    CMonocyte = character(0),
    FOB = character(0),
    GC_B = character(0),
    ILC3 = character(0),
    NMonocyte = character(0),
    TAM = character(0),
    Thelper = c("T helper cells", "Tfh", "Th1 cells", "Th17 cells", "Th2 cells"),
    Treg = c("Tregs", "nTreg", "iTregs"),
    cDC1 = character(0),
    cDC2 = character(0),
    cytotoxicNK = "NK CD56dim cells",
    exhausted_T = "Exhausted T cells",
    macrophage = c("Macrophages", "M1 Macrophages", "M2 Macrophages"),
    mast_cell = "Mast cell",
    pDC = "pDCs",
    plasma = "Plasma cells",
    regulatoryNK = "NK CD56bright cells"
  )
  out <- lapply(as.character(target), function(x) if (x %in% names(mapper)) mapper[[x]] else x)
  names(out) <- target
  out
}

immucellai_target_mapper <- function(target) {
  mapper <- list(
    Bcells = "Bcell",
    CD4Tcm = character(0),
    CD4Tem = character(0),
    CD4Temra = character(0),
    CD4Tfh = "Tfh",
    CD4Tnaive = "CD4_naive",
    CD4Trm = character(0),
    CD8Tcm = character(0),
    CD8Tem = character(0),
    CD8Temra = character(0),
    CD8Tnaive = "CD8_naive",
    CD8Trm = character(0),
    MAIT = "MAIT",
    Monocyte = "Monocyte",
    NK = "NK",
    NKT = "NKT",
    Th1 = "Th1",
    Th17 = "Th17",
    monoDC = character(0),
    `γδT` = "Gamma_delta",
    CD4Tmemory = character(0),
    CD4Tn = "CD4_naive",
    CD8Tmemory = character(0),
    CD8Tn = "CD8_naive",
    CMonocyte = character(0),
    FOB = character(0),
    GC_B = character(0),
    ILC3 = character(0),
    NMonocyte = character(0),
    TAM = character(0),
    Thelper = c("Th1", "Th2", "Th17", "Tfh"),
    Treg = c("nTreg", "iTreg", "Tr1"),
    cDC1 = character(0),
    cDC2 = character(0),
    cytotoxicNK = character(0),
    exhausted_T = "Exhausted",
    macrophage = "Macrophage",
    mast_cell = character(0),
    pDC = character(0),
    plasma = character(0),
    regulatoryNK = character(0)
  )
  out <- lapply(as.character(target), function(x) if (x %in% names(mapper)) mapper[[x]] else x)
  names(out) <- target
  out
}

lm22_mapper <- function(dataset) {
  if(dataset=='GSE176078') return(list(Bmemory='B cells memory',Bnaive='B cells naive',
    CD4Tnaive='T cells CD4 naive',M1='Macrophages M1',Monocyte='Monocytes',
    NK=c('NK cells resting','NK cells activated'),Tfh='T cells follicular helper',Treg='T cells regulatory (Tregs)'))
  if(dataset=='GSE146771') return(list(CD4Tmemory=c('T cells CD4 memory resting','T cells CD4 memory activated'),
    CD4Tn='T cells CD4 naive',macrophage=c('Macrophages M0','Macrophages M1','Macrophages M2'),
    mast_cell=c('Mast cells resting','Mast cells activated'),plasma='Plasma cells',Treg='T cells regulatory (Tregs)'))
  stopifnot(dataset=='GSE164522')
  m<-list(Bcells=c('B cells naive','B cells memory'),CD4Tfh='T cells follicular helper',
    CD4Tnaive='T cells CD4 naive',Monocyte='Monocytes',NK=c('NK cells resting','NK cells activated'),
    Treg='T cells regulatory (Tregs)',mast_cell=c('Mast cells resting','Mast cells activated'),plasma='Plasma cells')
  m[[paste0(intToUtf8(c(947,948)),'T')]]<-'T cells gamma delta'
  m
}

lm22_target_mapper <- function(target) {
  out <- lm22_mapper("GSE164522")[target]
  names(out) <- target
  out
}
