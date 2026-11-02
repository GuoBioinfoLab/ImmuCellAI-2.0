patch_immucellai_gsva <- function() {
    if (!isTRUE(run_immucellai))
        return(invisible(FALSE))
    ns <- asNamespace("ImmuCellAI")
    patched_sample_abundance <- function(sample, paper_marker, marker_exp, data_type, customer) {
        sam <- apply(sample, 2, as.numeric)
        row.names(sam) <- row.names(sample)
        tt <- intersect(row.names(sam), as.vector(unlist(paper_marker)))
        genes <- intersect(tt, row.names(marker_exp))
        sam_exp <- as.matrix(sam[genes, , drop = FALSE])
        colnames(sam_exp) <- colnames(sample)
        marker_exp <- marker_exp[genes, , drop = FALSE]
        marker_tag_mat <- c()
        for (cell in names(paper_marker)) {
            tag <- marker_tag(genes, as.vector(unlist(paper_marker[cell])))
            marker_tag_mat <- cbind(marker_tag_mat, tag)
        }
        row.names(marker_tag_mat) <- row.names(marker_exp)
        colnames(marker_tag_mat) <- names(paper_marker)
        exp_new <- apply(sam_exp, 2, function(x) sample_ratio(x, marker_exp, marker_tag_mat, data_type))
        if (requireNamespace("GSVA", quietly = TRUE) && exists("ssgseaParam", envir = asNamespace("GSVA"))) {
            param <- GSVA::ssgseaParam(exprData = as.matrix(exp_new), geneSets = paper_marker, normalize = TRUE, verbose = FALSE)
            result <- GSVA::gsva(param)
        }
        else {
            result <- GSVA::gsva(as.matrix(exp_new), paper_marker, method = "ssgsea", ssgsea.norm = TRUE, parallel.sz = 20)
        }
        if (ncol(result) < 3) {
            result[which(result < 0)] <- 0
        }
        else {
            result <- result - apply(result, 1, min)
        }
        compensation_matrix_num <- apply(compensation_matrix, 2, as.numeric)
        if (customer == 0) {
            row.names(compensation_matrix_num) <- row.names(compensation_matrix)
            result_norm <- compensation(result, compensation_matrix_num)
            if (ncol(result_norm) == 1) {
                InfiltrationScore <- sum(result_norm)
            }
            else {
                if (nrow(result_norm) == 24) {
                  InfiltrationScore <- apply(result_norm[c("Bcell", "CD4_T", "CD8_T", "DC", "Macrophage", "Monocyte", "Neutrophil",
                    "NK"), ], 2, sum)
                }
                else {
                  InfiltrationScore <- apply(result_norm, 2, sum)
                }
            }
            InfiltrationScore <- (InfiltrationScore/max(InfiltrationScore)) * 0.9
            result_mat <- rbind(result_norm, InfiltrationScore)
        }
        else {
            result_mat <- result
        }
        if (ncol(result_mat) > 1) {
            result_mat <- apply(result_mat, 1, function(x) round(x, 3))
        }
        else {
            result_mat <- t(round(result_mat, 3))
        }
        result_mat <- as.matrix(result_mat)
        if (ncol(result_mat) == 1) {
            colnames(result_mat) <- colnames(sample)
            T_FRE <- t(result_mat)
        }
        else {
            T_FRE <- result_mat
        }
        T_FRE
    }
    environment(patched_sample_abundance) <- ns
    was_locked <- bindingIsLocked("Sample_abundance_calculation", ns)
    if (was_locked)
        unlockBinding("Sample_abundance_calculation", ns)
    assign("Sample_abundance_calculation", patched_sample_abundance, envir = ns)
    if (was_locked)
        lockBinding("Sample_abundance_calculation", ns)
    invisible(TRUE)
}

select_signature_genes <- function(reference, n = 1000) {
    ref <- log2(as.matrix(reference) + 1)
    state_var <- apply(ref, 1, stats::var, na.rm = TRUE)
    state_max <- apply(ref, 1, max, na.rm = TRUE)
    state_mean <- rowMeans(ref, na.rm = TRUE)
    score <- state_var * (state_max + 1e-08)/(state_mean + 1e-08)
    score[!is.finite(score)] <- 0
    genes <- names(sort(score, decreasing = TRUE))
    genes[seq_len(min(n, length(genes)))]
}
