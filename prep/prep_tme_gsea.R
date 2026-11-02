# =============================================================================
# DLBCL-PRISM - pathway (GSEA) and tumor-microenvironment precomputation
#  1. Pre-ranked GSEA (fgsea; MSigDB Hallmark + Reactome) on the pooled
#     meta-analysis statistic (logHR / SE) of every analysis in the app
#  2. Per-patient microenvironment scores (mean within-cohort z-score of
#     canonical marker genes) and per-cohort Cox models of each score
# Run from the DLBCL-PRISM folder after prep_data.R and prep_extra.R:
#   Rscript prep/prep_tme_gsea.R
# =============================================================================
suppressPackageStartupMessages({ library(fgsea); library(msigdbr); library(survival) })
OUT <- "docs/data"
info <- readRDS(file.path(OUT, "meta_info.rds"))
CURATED_DROP <- c(names(info$overlap), info$non_dlbcl)
KMIN_SURV <- 5; KMIN_EF <- 4          # genes entering the GSEA ranking
set.seed(2026)

# ---- 1. GSEA ---------------------------------------------------------------------
mk <- function(x) lapply(split(toupper(x$gene_symbol), x$gs_name), unique)
sets <- list(Hallmark = mk(msigdbr(species = "Homo sapiens", collection = "H")),
             Reactome = mk(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")))

files <- c(sprintf("pooled_%s_%s_%s", rep(c("OS", "PFS"), each = 4), rep(rep(c("uni", "multi"), each = 2), 2), c("curated", "all")),
           "pooled_EF_curated", "pooled_EF_all")
universe <- character(0); res <- list()
for (f in files) {
  p <- readRDS(file.path(OUT, paste0(f, ".rds")))
  p <- p[p$k >= (if (grepl("_EF_", f)) KMIN_EF else KMIN_SURV) & is.finite(p$se) & p$se > 0, ]
  st <- setNames(p$logHR / p$se, p$key); st <- st[!duplicated(names(st))]
  universe <- union(universe, names(st))
  res[[sub("pooled_", "", f)]] <- do.call(rbind, lapply(names(sets), function(cl) {
    r <- suppressWarnings(fgsea(sets[[cl]], st, minSize = 15, maxSize = 500, eps = 0))
    data.frame(collection = cl, pathway = r$pathway, size = r$size, ES = signif(r$ES, 4), NES = signif(r$NES, 4),
               p = signif(r$pval, 4), fdr = signif(r$padj, 4),
               leading = vapply(r$leadingEdge, paste, "", collapse = ", "), stringsAsFactors = FALSE)
  }))
  cat(sprintf("%-18s genes=%5d  Hallmark FDR<0.05: %2d  Reactome FDR<0.05: %3d\n", f, length(st),
              sum(res[[length(res)]]$fdr < 0.05 & res[[length(res)]]$collection == "Hallmark"),
              sum(res[[length(res)]]$fdr < 0.05 & res[[length(res)]]$collection == "Reactome")))
}
# gene sets stored as integer indices into a shared gene vector (small download)
used <- unique(unlist(lapply(res, `[[`, "pathway")))
universe <- sort(universe)
idx <- lapply(sets, function(s) { s <- s[names(s) %in% used]; lapply(s, function(g) sort(match(intersect(g, universe), universe))) })
saveRDS(list(genes = universe, sets = idx, kmin = c(surv = KMIN_SURV, ef = KMIN_EF),
             msigdb = as.character(packageVersion("msigdbr"))), file.path(OUT, "gsea_sets.rds"), compress = "xz")
saveRDS(res, file.path(OUT, "gsea_results.rds"), compress = "xz")

# ---- 2. Microenvironment scores ---------------------------------------------------------
TME <- list(
  "T cells"                  = c("CD3D", "CD3E", "CD3G", "CD2", "CD5", "CD6", "CD28", "LCK", "TRAT1"),
  "CD8 T / cytotoxic"        = c("CD8A", "CD8B", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "NKG7", "GNLY"),
  "Regulatory T cells"       = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "CCR8"),
  "NK cells"                 = c("KLRD1", "KLRF1", "NCR1", "SH2D1B", "KIR2DL1", "KIR3DL1"),
  "Macrophages"              = c("CD68", "CD14", "CSF1R", "C1QA", "C1QB", "C1QC", "AIF1", "FCGR1A"),
  "M2-like macrophages"      = c("CD163", "MRC1", "MSR1", "CD209", "F13A1", "STAB1", "SIGLEC1", "LYVE1"),
  "Dendritic cells"          = c("CD1C", "CLEC10A", "FCER1A", "CLEC9A", "XCR1"),
  "Neutrophils"              = c("FCGR3B", "CXCR1", "CXCR2", "CSF3R", "FPR1", "CEACAM8"),
  "IFN-gamma / inflammation" = c("IFNG", "STAT1", "CXCL9", "CXCL10", "CXCL11", "IDO1", "GBP1", "GBP5"),
  "Immune checkpoints"       = c("CD274", "PDCD1LG2", "PDCD1", "LAG3", "HAVCR2", "TIGIT"),
  "MHC class II"             = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DMA", "HLA-DMB", "CD74", "CIITA"),
  "Matrix / CAF (stromal-1)" = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A2", "FN1", "SPARC", "LUM", "DCN", "FAP", "PDGFRB", "MMP2"),
  "Endothelium (stromal-2)"  = c("PECAM1", "VWF", "CDH5", "KDR", "ENG", "TEK", "ESAM", "EMCN"),
  "Proliferation"            = c("MKI67", "TOP2A", "CCNB1", "CDK1", "BUB1", "AURKA", "CENPF", "MCM2"),
  "B-cell program (tumor)"  = c("MS4A1", "CD19", "CD79A", "CD79B", "CD22", "PAX5"))

cl <- readRDS(file.path(OUT, "clinical.rds")); gidx <- readRDS(file.path(OUT, "gene_index.rds"))
decode <- function(r) { v <- (as.integer(r) - 128) / 10; v[as.integer(r) == 0L] <- NA; v }
allg <- unique(unlist(TME)); gi <- gidx[gidx$key %in% allg & gidx$expr, ]
Z <- matrix(NA_real_, nrow(cl), length(allg), dimnames = list(NULL, allg))
for (ch in unique(gi$chunk)) {
  M <- readRDS(file.path(OUT, "chunks", sprintf("c%04d.rds", ch)))$expr
  for (g in intersect(gi$key[gi$chunk == ch], rownames(M))) Z[, g] <- decode(M[g, ])
}
cat("Marker genes found:", sum(colSums(!is.na(Z)) > 0), "/", length(allg), "\n")

scores <- matrix(NA_real_, nrow(cl), length(TME), dimnames = list(NULL, names(TME)))
coverage <- list()
for (coh in unique(cl$cohort)) {
  i <- which(cl$cohort == coh)
  present <- colMeans(!is.na(Z[i, , drop = FALSE])) >= 0.8
  for (s in names(TME)) {
    g <- intersect(TME[[s]], allg[present]); coverage[[length(coverage) + 1]] <- data.frame(cohort = coh, cell = s, n_genes = length(g), of = length(TME[[s]]))
    if (length(g) < ceiling(length(TME[[s]]) / 2)) next
    v <- rowMeans(Z[i, g, drop = FALSE], na.rm = TRUE)
    scores[i, s] <- (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  }
}
coverage <- do.call(rbind, coverage)
saveRDS(list(scores = round(scores, 3), genes = TME, coverage = coverage), file.path(OUT, "tme_scores_v2.rds"), compress = "xz")

# per-cohort Cox of each score: endpoint x model x cohort set
cox_rows <- list()
for (set in c("curated", "all")) {
  keep <- if (set == "all") rep(TRUE, nrow(cl)) else !(cl$cohort %in% CURATED_DROP) & !(cl$histology %in% c("BL", "PMBL"))
  for (ep in c("OS", "PFS")) for (md in c("uni", "multi")) for (coh in unique(cl$cohort[keep])) for (s in names(TME)) {
    i <- keep & cl$cohort == coh
    d <- data.frame(t = cl[[if (ep == "OS") "os_m" else "pfs_m"]][i], e = cl[[if (ep == "OS") "os_e" else "pfs_e"]][i],
                    x = scores[i, s], age = cl$age[i], ipi = cl$ipi_risk[i])
    f <- "Surv(t, e) ~ x"
    if (md == "multi") {
      if (sum(!is.na(d$age)) >= 0.8 * nrow(d)) f <- paste(f, "+ age")
      if (sum(!is.na(d$ipi)) >= 0.8 * nrow(d) && length(unique(na.omit(d$ipi))) > 1) f <- paste(f, "+ ipi")
      if (f == "Surv(t, e) ~ x") next
    }
    d <- d[complete.cases(d[, all.vars(as.formula(f))]), ]
    if (nrow(d) < 20 || sum(d$e) < 5) next
    fit <- tryCatch(coxph(as.formula(f), d), error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit)) next
    cox_rows[[length(cox_rows) + 1]] <- data.frame(set = set, endpoint = ep, model = md, cohort = coh, cell = s,
      logHR = unname(coef(fit)["x"]), se = sqrt(vcov(fit)["x", "x"]), n = nrow(d), events = sum(d$e))
  }
}
cx <- do.call(rbind, cox_rows); cx <- cx[is.finite(cx$se) & cx$se < 5, ]
saveRDS(cx, file.path(OUT, "tme_cox_v2.rds"), compress = "xz")
print(aggregate(cohort ~ cell, cx[cx$set == "curated" & cx$endpoint == "OS" & cx$model == "uni", ], length))
cat("Done\n")
