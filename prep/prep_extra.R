# =============================================================================
# DLBCL-PRISM - extra precomputation (run after prep_data.R)
#  1. Pooled OS/PFS results with REML + Hartung-Knapp-Sidik-Jonkman (HKSJ)
#  2. Genome-wide early treatment failure (event < 24 months) in R-CHOP-based
#     patients: per-cohort logistic regression, pooled by REML + HKSJ
# Run from the DLBCL-PRISM folder:  Rscript prep/prep_extra.R
# =============================================================================
IN  <- ".."
OUT <- "docs/data"
info <- readRDS(file.path(OUT, "meta_info.rds"))
CURATED_DROP <- c(names(info$overlap), info$non_dlbcl)

# REML + HKSJ (untruncated, as metafor test = "knha") for many genes at once (g = gene key per row)
reml_hksj_grouped <- function(y, v, g, maxit = 100, tol = 1e-8) {
  g <- factor(g); k <- tabulate(g)
  sumg <- function(x) as.vector(rowsum(as.numeric(x), g, reorder = TRUE))
  w0 <- 1 / v; mu0 <- sumg(w0 * y) / sumg(w0)
  Q <- sumg(w0 * (y - mu0[g])^2)
  tau2 <- pmax(0, (Q - (k - 1)) / (sumg(w0) - sumg(w0^2) / sumg(w0)))
  tau2[k < 2 | !is.finite(tau2)] <- 0
  for (it in seq_len(maxit)) {
    w <- 1 / (v + tau2[g]); mu <- sumg(w * y) / sumg(w)
    new <- pmax(0, sumg(w^2 * ((y - mu[g])^2 - v)) / sumg(w^2) + 1 / sumg(w))
    new[k < 2 | !is.finite(new)] <- 0
    if (max(abs(new - tau2)) < tol) { tau2 <- new; break }
    tau2 <- new
  }
  w <- 1 / (v + tau2[g]); mu <- sumg(w * y) / sumg(w)
  qhk <- ifelse(k > 1, sumg(w * (y - mu[g])^2) / (k - 1), 1)
  se <- sqrt(qhk / sumg(w))                                   # untruncated HKSJ (metafor test = "knha")
  p <- ifelse(k > 1, 2 * pt(-abs(mu / se), pmax(k - 1, 1)), 2 * pnorm(-abs(mu / se)))
  vt <- (k - 1) * sumg(w0) / (sumg(w0)^2 - sumg(w0^2))
  data.frame(key = levels(g), k = k, logHR = mu, se = se, p = p,
             I2 = ifelse(k > 1, 100 * tau2 / (tau2 + vt), NA), tau2 = tau2,
             pQ = ifelse(k > 1, pchisq(Q, k - 1, lower.tail = FALSE), NA),
             cons = sumg(sign(y) == sign(mu[g])) / k, stringsAsFactors = FALSE)
}
finish <- function(r, npat) {
  r$npat <- as.vector(npat[r$key]); r$fdr <- p.adjust(r$p, "BH")
  num <- c("logHR", "se", "p", "I2", "tau2", "pQ", "cons", "fdr"); r[num] <- lapply(r[num], signif, 4); r
}

# ---- 1. OS / PFS pooled with HKSJ ----------------------------------------------
cm <- read.csv(file.path(IN, "meta_results/cohort_level_all_models.csv"), stringsAsFactors = FALSE)
cm$key <- toupper(cm$gene)
cm <- cm[!duplicated(cm[, c("key", "cohort", "endpoint", "model")]) & is.finite(cm$logHR) & is.finite(cm$se) & cm$se > 0, ]
for (ep in c("OS", "PFS")) for (md in c("uni", "multi")) {
  d <- cm[cm$endpoint == ep & cm$model == (if (md == "uni") "univariate" else "multivariate"), ]
  for (set in c("all", "curated")) {
    dd <- if (set == "all") d else d[!d$cohort %in% CURATED_DROP, ]
    r <- finish(reml_hksj_grouped(dd$logHR, dd$se^2, dd$key), tapply(dd$n, dd$key, sum))
    saveRDS(r, file.path(OUT, sprintf("pooled_%s_%s_%s.rds", ep, md, set)), compress = "xz")
    cat(sprintf("%s %s %-8s genes=%d  FDR<0.05 & k>=8: %d\n", ep, md, set, nrow(r), sum(r$fdr < 0.05 & r$k >= 8)))
  }
}

# ---- 2. Early failure on R-CHOP: genome-wide ------------------------------------
cl <- readRDS(file.path(OUT, "clinical.rds"))
gidx <- readRDS(file.path(OUT, "gene_index.rds"))
ef <- !is.na(cl$early_failure) & cl$tx_group %in% c("R-chemo", "R-chemo + bortezomib")   # R-CHOP-based (incl. RB-CHOP)
y_all <- cl$early_failure == "Yes"
print(table(cl$cohort[ef], cl$early_failure[ef]))

fit_set <- function(keep) {
  cohorts <- names(which(tapply(y_all[keep], cl$cohort[keep], function(v) sum(v) >= 5 && sum(!v) >= 5)))
  rows <- list()
  for (ch in sort(unique(gidx$chunk[gidx$expr]))) {
    M <- readRDS(file.path(OUT, "chunks", sprintf("c%04d.rds", ch)))$expr
    for (coh in cohorts) {
      i <- which(keep & cl$cohort == coh); y <- y_all[i]
      Z <- M[, i, drop = FALSE]
      for (g in rownames(M)) {
        r <- as.integer(Z[g, ]); ok <- r != 0L
        if (sum(ok) < 20 || sum(y[ok]) < 5 || sum(!y[ok]) < 5) next
        z <- (r[ok] - 128) / 10; if (sd(z) == 0) next
        f <- suppressWarnings(glm.fit(cbind(1, z), y[ok], family = binomial()))
        if (!f$converged) next
        V <- tryCatch(chol2inv(f$qr$qr[1:2, 1:2]), error = function(e) NULL); if (is.null(V)) next
        rows[[length(rows) + 1]] <- c(f$coefficients[2], sqrt(V[2, 2]), sum(ok))
        names(rows)[length(rows)] <- paste(g, coh, sep = "|")
      }
    }
    if (ch %% 25 == 0) cat("  chunk", ch, "\n")
  }
  m <- do.call(rbind, rows); kk <- do.call(rbind, strsplit(rownames(m), "|", fixed = TRUE))
  d <- data.frame(key = kk[, 1], cohort = kk[, 2], logOR = m[, 1], se = m[, 2], n = m[, 3], stringsAsFactors = FALSE)
  d <- d[is.finite(d$se) & d$se > 0 & d$se < 10, ]
  finish(reml_hksj_grouped(d$logOR, d$se^2, d$key), tapply(d$n, d$key, sum))
}
keep_cur <- ef & !(cl$cohort %in% CURATED_DROP) & !(cl$histology %in% c("BL", "PMBL"))
for (set in c("curated", "all")) {
  r <- fit_set(if (set == "all") ef else keep_cur)
  saveRDS(r, file.path(OUT, sprintf("pooled_EF_%s.rds", set)), compress = "xz")
  cat(sprintf("EF %-8s genes=%d  max k=%d  FDR<0.05: %d\n", set, nrow(r), max(r$k), sum(r$fdr < 0.05)))
}
cat("Done\n")
