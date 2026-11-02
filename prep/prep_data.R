# =============================================================================
# DLBCL-PRISM - data preparation
# Builds the compact, web-ready data used by the Shinylive app (docs/data/).
# Run from the DLBCL-PRISM folder:  Rscript prep/prep_data.R
# Inputs: ../clinico, ../coortes, ../meta_results
# =============================================================================
suppressMessages({ library(SummarizedExperiment) })
set.seed(1)
IN  <- ".."
OUT <- "docs/data"
dir.create(file.path(OUT, "chunks"), recursive = TRUE, showWarnings = FALSE)
CHUNK_SIZE <- 150

# ---- cohorts measured in days (converted to months) --------------------------
DAYS <- c("CGCI-HTMCP-DLBCL","CTSP-DLBCL1","E-MEXP-3488","GSE125966","GSE136971",
          "GSE248835","GSE273595","GSE34171","GSE69049","GSE90953","GSE93984","NCICCR-DLBCL")
# ---- redundant cohorts (patient overlap detected by time/event/age matching) --
#  GSE87371 == GSE136971 (same 221 patients, LNH03-6B); GSE11318 & GSE53786 subsets of GSE10846
OVERLAP <- c(GSE87371 = "Same patients as GSE136971",
             GSE11318 = "Patients largely contained in GSE10846",
             GSE53786 = "Patients largely contained in GSE10846")
NON_DLBCL <- c("CGCI-BLGSP")   # Burkitt Lymphoma Genome Sequencing Project

# =============================================================================
# 1. Clinical harmonisation
# =============================================================================
cl <- read.csv(file.path(IN, "clinico/ALL_COHORTS_clinical.csv"), stringsAsFactors = FALSE)
cl[cl == ""] <- NA
k <- cl$dataset %in% DAYS
tm <- function(x) ifelse(k, x / 30.4375, x)

coo_raw <- toupper(trimws(cl$coo))
coo <- ifelse(grepl("^ABC", coo_raw), "ABC",
       ifelse(grepl("^GC", coo_raw), "GCB",
       ifelse(coo_raw %in% c("UNC","UC","NC","UNCLASSIFIED DLBCL","TYPE III"), "Unclassified", NA)))
histology <- ifelse(coo_raw %in% "BL", "BL", ifelse(coo_raw %in% "PMBL", "PMBL", "DLBCL"))
histology[cl$dataset %in% NON_DLBCL] <- "BL"
histology[!is.na(cl$hist_group) & cl$hist_group == "HGBL"] <- "HGBL"

ipi_num <- suppressWarnings(as.numeric(cl$ipi))
ipi_txt <- cl$ipi
ipi_risk <- ifelse(!is.na(ipi_num), cut(ipi_num, c(-1, 1, 3, 5), labels = FALSE), NA)
ipi_risk <- c("Low (0-1)", "Intermediate (2-3)", "High (4-5)")[ipi_risk]
ipi_risk[ipi_txt %in% "Low"] <- "Low (0-1)"
ipi_risk[ipi_txt %in% c("Low-Intermediate","High-Intermediate","Medium")] <- "Intermediate (2-3)"
ipi_risk[ipi_txt %in% "High"] <- "High (4-5)"

st <- toupper(gsub("STAGE| ", "", cl$stage))
st <- sub("^(IV|III|II|I)[A-Z]*$", "\\1", st)
st <- c("1" = "I", "2" = "II", "3" = "III", "4" = "IV", I = "I", II = "II", III = "III", IV = "IV")[st]
ecog <- suppressWarnings(as.numeric(cl$ecog)); ecog[!ecog %in% 0:4] <- NA

ldh <- suppressWarnings(as.numeric(cl$ldh)); ldh_hi <- rep(NA, nrow(cl))
for (d in unique(cl$dataset[!is.na(ldh)])) {
  i <- cl$dataset == d & !is.na(ldh); mx <- max(ldh[i])
  if (mx <= 1) ldh_hi[i] <- ldh[i] == 1                         # already binary
  else if (median(ldh[i]) < 3) ldh_hi[i] <- ldh[i] > 1          # ratio to ULN
}                                                               # U/L with unknown ULN -> NA

sex <- toupper(substr(cl$gender, 1, 1)); sex <- c(M = "Male", F = "Female")[sex]

tx <- cl$treatment; txs <- cl$treatment_std; arm <- cl$treatment_arm
tx_group <- rep(NA_character_, nrow(cl))
is_R <- grepl("(^|[^A-Z])R-|-R($|[^A-Z])|RITU|Ritux|^R$|^R-|CHOP-R|CVP-R", tx) | txs %in% "R"
tx_group[grepl("CHOP|CVP|COP|ACVBP|^R$", tx) | txs %in% c("R","noR")] <- "Chemo without R"
tx_group[is_R] <- "R-chemo"
tx_group[grepl("RB-CHOP|Bortezomib", tx)] <- "R-chemo + bortezomib"
tx_group[grepl("R-ACVBP", tx)] <- "R-ACVBP"
tx_group[grepl("CODOX|ABVD|DHAP|IDARAM", tx)] <- "Other chemo"
tx_group[grepl("palliat|Supportive|radiotherapy|resection|steroid", tx, ignore.case = TRUE)] <- "Palliative/local"
tx_group[arm %in% "Axicabtagene Ciloleucel"] <- "CAR-T (axi-cel)"
tx_group[arm %in% "Standard of Care Chemotherapy"] <- "2L SOC chemo"
# R-CHOP by study design, treatment not in the GEO annotation:
# GSE31312 (International DLBCL R-CHOP Consortium), GSE21846 (Montes-Moreno, Blood 2011)
tx_group[is.na(tx_group) & cl$dataset %in% c("GSE31312", "GSE21846")] <- "R-chemo"
rituximab <- ifelse(tx_group %in% c("R-chemo","R-chemo + bortezomib","R-ACVBP"), "Yes",
             ifelse(tx_group %in% "Chemo without R", "No", NA))

os_m <- tm(cl$os_time); os_e <- cl$os_event
pfs_m <- tm(cl$pfs_time); pfs_e <- cl$pfs_event; pfs_src <- ifelse(is.na(pfs_m), NA, "PFS")
j <- is.na(pfs_m) & !is.na(cl$efs_time)
pfs_m[j] <- tm(cl$efs_time)[j]; pfs_e[j] <- cl$efs_event[j]; pfs_src[j] <- "EFS"
early <- ifelse(!is.na(pfs_e) & pfs_e == 1 & pfs_m < 24, "Yes", ifelse(!is.na(pfs_m) & pfs_m >= 24, "No", NA))
os24  <- ifelse(!is.na(os_e) & os_e == 1 & os_m < 24, "Yes", ifelse(!is.na(os_m) & os_m >= 24, "No", NA))

clin <- data.frame(
  sample = cl$sample_id, cohort = cl$dataset, platform = cl$platform,
  os_m = round(os_m, 3), os_e = os_e, pfs_m = round(pfs_m, 3), pfs_e = pfs_e, pfs_src = pfs_src,
  coo = coo, histology = histology, ipi_score = ipi_num, ipi_risk = ipi_risk,
  stage = st, adv_stage = ifelse(is.na(st), NA, ifelse(st %in% c("III","IV"), "III-IV", "I-II")),
  ecog = ecog, ecog2 = ifelse(is.na(ecog), NA, ifelse(ecog >= 2, ">=2", "0-1")),
  ldh_high = ifelse(is.na(ldh_hi), NA, ifelse(ldh_hi, "Elevated", "Normal")),
  age = round(cl$age, 1), age60 = ifelse(is.na(cl$age), NA, ifelse(cl$age > 60, ">60", "<=60")),
  sex = unname(sex), tx_group = tx_group, rituximab = rituximab,
  early_failure = early, death24 = os24, stringsAsFactors = FALSE)
clin$ipi_risk <- factor(clin$ipi_risk, c("Low (0-1)", "Intermediate (2-3)", "High (4-5)"))
clin$stage <- factor(clin$stage, c("I","II","III","IV"))
clin$coo <- factor(clin$coo, c("GCB","ABC","Unclassified"))
cat("Clinical:", nrow(clin), "samples\n"); print(table(clin$tx_group, useNA = "ifany"))

# =============================================================================
# 2. Expression -> within-cohort z-scores, int8-coded (z*10, NA = 0)
# =============================================================================
files <- list.files(file.path(IN, "coortes"), pattern = "_SummarizedExperiment.rds$", full.names = TRUE)
zl <- list(); gcount <- list()
for (f in files) {
  coh <- sub("_SummarizedExperiment.rds$", "", basename(f))
  se <- readRDS(f); x <- as.matrix(assay(se)); storage.mode(x) <- "double"
  if (max(x, na.rm = TRUE) > 100) x <- log2(pmax(x, 0) + 1)           # non-log data (GSE34171)
  key <- toupper(rownames(x)); keep <- !duplicated(key); x <- x[keep, , drop = FALSE]; rownames(x) <- key[keep]
  idx <- match(colnames(x), clin$sample[clin$cohort == coh])
  if (anyNA(idx)) stop("Unmatched samples in ", coh)
  mu <- rowMeans(x, na.rm = TRUE); sdv <- apply(x, 1, sd, na.rm = TRUE)
  z <- (x - mu) / sdv; z[!is.finite(z)] <- NA; z <- z[sdv > 0 & !is.na(sdv), , drop = FALSE]
  zl[[coh]] <- z; gcount[[coh]] <- rownames(z)
  cat(sprintf("%-18s %5d samples %6d genes\n", coh, ncol(z), nrow(z)))
}
gtab <- table(unlist(gcount))
expr_genes <- sort(names(gtab)[gtab >= 2])

# =============================================================================
# 3. Cohort-level Cox estimates
# =============================================================================
cm <- read.csv(file.path(IN, "meta_results/cohort_level_all_models.csv"), stringsAsFactors = FALSE)
cm$key <- toupper(cm$gene)
cm <- cm[!duplicated(cm[, c("key","cohort","endpoint","model")]), ]
cm <- data.frame(key = cm$key, cohort = cm$cohort, endpoint = cm$endpoint,
                 model = ifelse(cm$model == "univariate", "uni", "multi"),
                 logHR = signif(cm$logHR, 5), se = signif(cm$se, 5), n = cm$n, nevt = cm$nevt,
                 stringsAsFactors = FALSE)
cm <- cm[is.finite(cm$logHR) & is.finite(cm$se) & cm$se > 0, ]

# =============================================================================
# 4. Pooled results: provided (all cohorts) + curated (overlap removed), REML
# =============================================================================
reml_grouped <- function(y, v, g, maxit = 100, tol = 1e-8) {
  g <- factor(g); k <- tabulate(g)
  sumg <- function(x) as.vector(rowsum(as.numeric(x), g, reorder = TRUE))
  w0 <- 1 / v; mu0 <- sumg(w0 * y) / sumg(w0)
  Q <- sumg(w0 * (y - mu0[g])^2)
  C <- sumg(w0) - sumg(w0^2) / sumg(w0)
  tau2 <- pmax(0, (Q - (k - 1)) / C)                                   # DL start
  tau2[k < 2 | !is.finite(tau2)] <- 0
  for (it in seq_len(maxit)) {
    w <- 1 / (v + tau2[g]); mu <- sumg(w * y) / sumg(w)
    new <- pmax(0, sumg(w^2 * ((y - mu[g])^2 - v)) / sumg(w^2) + 1 / sumg(w))
    new[k < 2 | !is.finite(new)] <- 0
    if (max(abs(new - tau2)) < tol) { tau2 <- new; break }
    tau2 <- new
  }
  w <- 1 / (v + tau2[g]); mu <- sumg(w * y) / sumg(w); se <- sqrt(1 / sumg(w))
  vt <- (k - 1) * sumg(w0) / (sumg(w0)^2 - sumg(w0^2))                 # typical within-study var
  I2 <- ifelse(k > 1, 100 * tau2 / (tau2 + vt), NA)
  cons <- sumg(sign(y) == sign(mu[g])) / k
  data.frame(key = levels(g), k = k, logHR = mu, se = se, p = 2 * pnorm(-abs(mu / se)),
             I2 = I2, tau2 = tau2, pQ = ifelse(k > 1, pchisq(Q, k - 1, lower.tail = FALSE), NA),
             cons = cons, npat = sumg(rep(1, length(y))), stringsAsFactors = FALSE)
}
npat_sum <- function(d) tapply(d$n, d$key, sum)

for (ep in c("OS","PFS")) for (md in c("uni","multi")) {
  d <- cm[cm$endpoint == ep & cm$model == md, ]
  for (set in c("all","curated")) {
    dd <- if (set == "all") d else d[!d$cohort %in% c(names(OVERLAP), NON_DLBCL), ]
    r <- reml_grouped(dd$logHR, dd$se^2, dd$key)
    r$npat <- as.vector(npat_sum(dd)[r$key])
    r$fdr <- p.adjust(r$p, "BH")
    num <- c("logHR","se","p","I2","tau2","pQ","cons","fdr")
    r[num] <- lapply(r[num], signif, 4)
    saveRDS(r, file.path(OUT, sprintf("pooled_%s_%s_%s.rds", ep, md, set)), compress = "xz")
    if (set == "all") {                                                # validate vs provided file
      pv <- read.csv(file.path(IN, sprintf("meta_results/%s_%s_pooled.csv", ep,
                     if (md == "uni") "univariate" else "multivariate")))
      pv$key <- toupper(pv$gene); m <- merge(r, pv, by = "key")
      m <- m[m$k.x == m$k.y & m$k.x > 1, ]
      cat(sprintf("[validate %s %s] n=%d  max|dlogHR|=%.2e  cor(logHR)=%.6f\n", ep, md, nrow(m),
                  max(abs(m$logHR.x - m$logHR.y)), cor(m$logHR.x, m$logHR.y)))
    }
  }
}

# =============================================================================
# 5. Gene index + chunks (expression + cohort-level estimates)
# =============================================================================
keys <- sort(union(expr_genes, unique(cm$key)))
chunk <- ceiling(seq_along(keys) / CHUNK_SIZE)
gindex <- data.frame(key = keys, chunk = chunk, expr = keys %in% expr_genes,
                     n_expr = as.integer(ifelse(keys %in% names(gtab), gtab[keys], 0L)),
                     stringsAsFactors = FALSE)
saveRDS(gindex, file.path(OUT, "gene_index.rds"), compress = "xz")

N <- nrow(clin)
col_idx <- lapply(names(zl), function(coh) match(paste(coh, colnames(zl[[coh]])), paste(clin$cohort, clin$sample)))
names(col_idx) <- names(zl)
cm_split <- split(cm[, -1], factor(cm$key, levels = keys))
for (ch in unique(chunk)) {
  kk <- keys[chunk == ch]; ek <- kk[kk %in% expr_genes]
  M <- matrix(as.raw(0), nrow = length(ek), ncol = N, dimnames = list(ek, NULL))
  for (coh in names(zl)) {
    z <- zl[[coh]]; r <- intersect(ek, rownames(z)); if (!length(r)) next
    v <- round(z[r, , drop = FALSE] * 10) + 128; v[v < 1] <- 1; v[v > 255] <- 255; v[is.na(v)] <- 0
    M[r, col_idx[[coh]]] <- as.raw(v)
  }
  res <- cm_split[kk]; res <- res[vapply(res, nrow, 1L) > 0]
  saveRDS(list(expr = M, res = res), file.path(OUT, "chunks", sprintf("c%04d.rds", ch)), compress = "xz")
  if (ch %% 25 == 0) cat("chunk", ch, "/", max(chunk), "\n")
}

# =============================================================================
# 6. Cohort catalogue
# =============================================================================
ct <- do.call(rbind, lapply(split(clin, clin$cohort), function(s) data.frame(
  cohort = s$cohort[1], platform = s$platform[1], n = nrow(s),
  genes = if (s$cohort[1] %in% names(zl)) nrow(zl[[s$cohort[1]]]) else NA,
  os_n = sum(!is.na(s$os_m)), os_events = sum(s$os_e, na.rm = TRUE),
  pfs_n = sum(!is.na(s$pfs_m)), pfs_events = sum(s$pfs_e, na.rm = TRUE),
  pfs_type = paste(unique(na.omit(s$pfs_src)), collapse = "/"),
  coo = sum(!is.na(s$coo)), ipi = sum(!is.na(s$ipi_risk)), stage = sum(!is.na(s$stage)),
  ecog = sum(!is.na(s$ecog)), ldh = sum(!is.na(s$ldh_high)), age = sum(!is.na(s$age)),
  sex = sum(!is.na(s$sex)), treatment = sum(!is.na(s$tx_group)),
  time_unit_source = if (s$cohort[1] %in% DAYS) "days" else "months",
  note = unname(c(OVERLAP, `CGCI-BLGSP` = "Burkitt lymphoma (BLGSP)")[s$cohort[1]]),
  stringsAsFactors = FALSE)))
ct$note[is.na(ct$note)] <- ""
rownames(ct) <- NULL
saveRDS(clin, file.path(OUT, "clinical.rds"), compress = "xz")
saveRDS(ct, file.path(OUT, "cohorts.rds"), compress = "xz")
saveRDS(list(overlap = OVERLAP, non_dlbcl = NON_DLBCL, chunk_size = CHUNK_SIZE,
             built = format(Sys.time(), "%Y-%m-%d")), file.path(OUT, "meta_info.rds"))
cat("Done. Files:", length(list.files(OUT, recursive = TRUE)), "\n")
