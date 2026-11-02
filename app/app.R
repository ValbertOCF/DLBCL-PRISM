# =============================================================================
# DLBCL-PRISM
# Prognostic Resource for Integrated Survival Meta-analysis in DLBCL
# Shiny app designed to run fully in the browser via Shinylive (webR).
# =============================================================================
library(shiny)
library(bslib)
library(survival)
library(DT)

# Optional absolute URL of the data folder (leave "" to auto-detect, e.g.
# "https://<user>.github.io/<repo>/data/")
DATA_URL <- ""

# ---- palette ------------------------------------------------------------------
COL <- list(main = "#1f5f7a", risk = "#c0392b", prot = "#2471a3", grey = "#8c959f",
            light = "#d0d7de", ink = "#24292f", hi = "#b03a2e", lo = "#2e6da4")
GROUP_COLS <- c("#2e6da4", "#c0392b", "#e09f3e", "#5a9e6f", "#7d5ba6", "#8c959f")

# =============================================================================
# Data access: local folder when running in desktop R, HTTP fetch in webR
# =============================================================================
LOCAL_DIR <- Filter(dir.exists, c("data", "../docs/data", "docs/data"))[1]
.cache <- new.env()
.chunk_order <- character(0)

fetch_rds <- function(rel, base) {
  if (!is.null(.cache[[rel]])) return(.cache[[rel]])
  if (!is.na(LOCAL_DIR) && file.exists(file.path(LOCAL_DIR, rel))) {
    obj <- readRDS(file.path(LOCAL_DIR, rel))
  } else {
    tf <- tempfile(fileext = ".rds")
    on.exit(unlink(tf))
    utils::download.file(paste0(base, rel), tf, mode = "wb", quiet = TRUE)
    obj <- readRDS(tf)
  }
  if (startsWith(rel, "chunks/")) {                 # keep at most 40 chunks in memory
    .chunk_order <<- c(.chunk_order, rel)
    if (length(.chunk_order) > 40) { rm(list = .chunk_order[1], envir = .cache); .chunk_order <<- .chunk_order[-1] }
  }
  assign(rel, obj, envir = .cache)
  obj
}

# unannotated ("?|id") and multi-gene Affymetrix ("A///B") keys are hidden from the app
clean_key <- function(k) !grepl("?|", k, fixed = TRUE) & !grepl("///", k, fixed = TRUE)

decode_z <- function(r) { v <- (as.integer(r) - 128) / 10; v[as.integer(r) == 0L] <- NA; v }

# =============================================================================
# Statistics helpers
# =============================================================================
reml <- function(y, v, hksj = TRUE) {
  k <- length(y); w0 <- 1 / v
  if (k == 1) return(list(est = y, se = sqrt(v), lo = y - 1.96 * sqrt(v), hi = y + 1.96 * sqrt(v),
                          p = 2 * pnorm(-abs(y / sqrt(v))), tau2 = 0, I2 = NA, Q = NA, pQ = NA,
                          k = 1, pi_lo = NA, pi_hi = NA, w = 100))
  mu0 <- sum(w0 * y) / sum(w0); Q <- sum(w0 * (y - mu0)^2)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w0) - sum(w0^2) / sum(w0)))
  for (i in 1:200) {
    w <- 1 / (v + tau2); mu <- sum(w * y) / sum(w)
    new <- max(0, sum(w^2 * ((y - mu)^2 - v)) / sum(w^2) + 1 / sum(w))
    if (abs(new - tau2) < 1e-10) { tau2 <- new; break }
    tau2 <- new
  }
  w <- 1 / (v + tau2); mu <- sum(w * y) / sum(w); se <- sqrt(1 / sum(w)); crit <- 1.96
  if (hksj) {
    se <- sqrt(sum(w * (y - mu)^2) / (k - 1) / sum(w)); crit <- qt(0.975, k - 1)
    p <- 2 * pt(-abs(mu / se), k - 1)
  } else p <- 2 * pnorm(-abs(mu / se))
  vt <- (k - 1) * sum(w0) / (sum(w0)^2 - sum(w0^2))
  pis <- if (k >= 3) qt(0.975, k - 2) * sqrt(tau2 + se^2) else NA
  list(est = mu, se = se, lo = mu - crit * se, hi = mu + crit * se, p = p, tau2 = tau2,
       I2 = 100 * tau2 / (tau2 + vt), Q = Q, pQ = pchisq(Q, k - 1, lower.tail = FALSE), k = k,
       pi_lo = mu - pis, pi_hi = mu + pis, w = 100 * w / sum(w))
}

fmt_p <- function(p) vapply(p, function(x) if (is.na(x)) "NA" else if (x < 1e-4) formatC(x, format = "e", digits = 1)
                            else formatC(x, format = "f", digits = if (x < 0.01) 4 else 3), "")
fmt_hr <- function(e, l, h) sprintf("%.2f (%.2f-%.2f)", exp(e), exp(l), exp(h))

# pooled table -> HR with HKSJ CI (t, k-1 df); FDR recomputed among genes with k >= kmin
add_ci_fdr <- function(p, kmin) {
  p <- p[p$k >= kmin, ]; crit <- ifelse(p$k > 1, qt(0.975, pmax(p$k - 1, 1)), 1.96)
  p$HR <- exp(p$logHR); p$lo <- exp(p$logHR - crit * p$se); p$hi <- exp(p$logHR + crit * p$se)
  p$fdr <- p.adjust(p$p, "BH"); p
}

cox_row <- function(fit, term) {
  s <- summary(fit)$coefficients; ci <- summary(fit)$conf.int
  if (!term %in% rownames(s)) return(c(HR = NA, lo = NA, hi = NA, p = NA))
  c(HR = ci[term, 1], lo = ci[term, 3], hi = ci[term, 4], p = s[term, ncol(s)])
}

# within-cohort quantile groups
cut_within <- function(z, cohort, how) {
  probs <- switch(how, median = c(0, .5, 1), tertile = c(0, 1/3, 2/3, 1), quartile = c(0, .25, .5, .75, 1))
  labs <- switch(how, median = c("Low", "High"), tertile = c("Low", "Mid", "High"),
                 quartile = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"))
  out <- rep(NA_character_, length(z))
  for (c in unique(cohort)) {
    i <- which(cohort == c & !is.na(z)); if (length(i) < 10) next
    q <- quantile(z[i], probs, na.rm = TRUE)
    if (any(duplicated(q))) q <- quantile(z[i] + rnorm(length(i), 0, 1e-6), probs)
    out[i] <- as.character(cut(z[i], q, labels = labs, include.lowest = TRUE))
  }
  factor(out, labs)
}

auc <- function(score, y) {
  ok <- !is.na(score) & !is.na(y); s <- score[ok]; y <- y[ok]
  n1 <- sum(y); n0 <- sum(!y); if (n1 < 2 || n0 < 2) return(NA)
  (sum(rank(s)[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# =============================================================================
# Plot helpers (base graphics: fast in webR)
# =============================================================================
base_par <- function(...) par(family = "sans", las = 1, bty = "n", mgp = c(2.2, 0.6, 0),
                              tcl = -0.3, col.axis = COL$ink, col.lab = COL$ink, cex.axis = 0.85, ...)

# width of the widest string in inches; falls back to a character estimate when no plot is open
text_in <- function(x, cex = 1, font = 1) {
  x <- as.character(x); x <- x[!is.na(x)]; if (!length(x)) return(0)
  tryCatch(max(strwidth(x, "inches", cex = cex, font = font)), error = function(e) max(nchar(x)) * 0.085 * cex)
}

grp_cols <- function(n) if (n == 2) c(COL$lo, COL$hi) else if (n == 3) c(COL$lo, COL$grey, COL$hi) else c(COL$lo, "#7fa7c9", "#e3a19a", COL$hi)

empty_plot <- function(msg) {
  par(mar = c(0, 0, 0, 0)); plot.new(); text(0.5, 0.5, msg, col = "black", cex = 1.1)
}

forest_plot <- function(df, pooled = NULL, title = "", xlab = "Hazard ratio (per 1 SD increase)",
                        het = NULL, cols = c("Cohort", "N", "Events"), ref = 1, effect = "HR",
                        favors = c("favors low expression", "favors high expression")) {
  # df: label, n, nevt, est, lo, hi, w (est on log scale)
  k <- nrow(df); ytop <- k + 3.6
  rng <- range(c(df$lo, df$hi, pooled$lo, pooled$hi, pooled$pi_lo, pooled$pi_hi, 0), na.rm = TRUE, finite = TRUE)
  rng <- pmax(pmin(rng, log(20)), log(0.05)); rng <- rng + c(-0.08, 0.08) * diff(rng)
  w_in <- par("fin")[1]; base_par(mai = c(1.25, 0.29 * w_in, 0.55, 0.22 * w_in))
  # text columns placed from measured string widths so they never overlap
  sw <- function(x, cex = 0.85) text_in(x, cex, font = 2)
  x2 <- max(0.20 * w_in, 0.01 * w_in + sw(c(df$label, pooled$label)) + 0.12 + sw(c(cols[2], df$n, pooled$n)))
  x3 <- x2 + 0.15 + sw(c(cols[3], df$nevt, pooled$nevt))
  if (length(cols) < 3) x3 <- x2
  par(mai = c(1.25, max(0.29 * w_in, x3 + 0.2), 0.55, 0.22 * w_in))
  plot(NA, xlim = rng, ylim = c(0.2, ytop), axes = FALSE, xlab = "", ylab = "", yaxs = "i")
  mtext(xlab, side = 1, line = 2.1, cex = 0.9)
  ticks <- pretty(exp(rng), 6); ticks <- ticks[ticks > 0]
  if (diff(rng) > 2.5) ticks <- c(0.1, 0.25, 0.5, 1, 2, 4, 8)[c(0.1, 0.25, 0.5, 1, 2, 4, 8) >= exp(rng[1]) & c(0.1, 0.25, 0.5, 1, 2, 4, 8) <= exp(rng[2])]
  axis(1, at = log(ticks), labels = ticks, col = COL$grey)
  abline(v = log(ref), lty = 2, col = COL$grey)
  if (!is.null(pooled)) abline(v = pooled$est, lty = 3, col = COL$main)
  y <- rev(seq_len(k)) + 2
  xs <- grconvertX(c(0.01, x2 / w_in, x3 / w_in, 0.785, 0.985), "nfc", "user")
  hdr <- ytop - 0.35
  text(xs[1], hdr, cols[1], adj = 0, font = 2, xpd = NA, cex = 0.85)
  if (length(cols) > 1) text(xs[2], hdr, cols[2], adj = 1, font = 2, xpd = NA, cex = 0.85)
  if (length(cols) > 2) text(xs[3], hdr, cols[3], adj = 1, font = 2, xpd = NA, cex = 0.85)
  text(xs[4], hdr, paste0(effect, " (95% CI)"), adj = 0, font = 2, xpd = NA, cex = 0.85)
  if (!is.null(df$w)) text(xs[5], hdr, "Weight", adj = 1, font = 2, xpd = NA, cex = 0.85)
  segments(pmax(df$lo, rng[1]), y, pmin(df$hi, rng[2]), y, col = COL$ink, lwd = 1.2)
  sz <- if (!is.null(df$w)) 0.6 + 1.6 * sqrt(df$w / max(df$w, na.rm = TRUE)) else 1.2
  colp <- ifelse(df$est > 0, COL$risk, COL$prot)
  points(df$est, y, pch = 15, cex = sz, col = colp)
  text(xs[1], y, df$label, adj = 0, xpd = NA, cex = 0.8)
  if (!is.null(df$n)) text(xs[2], y, df$n, adj = 1, xpd = NA, cex = 0.8)
  if (!is.null(df$nevt) && length(cols) > 2) text(xs[3], y, df$nevt, adj = 1, xpd = NA, cex = 0.8)
  text(xs[4], y, fmt_hr(df$est, df$lo, df$hi), adj = 0, xpd = NA, cex = 0.8)
  if (!is.null(df$w)) text(xs[5], y, sprintf("%.1f%%", df$w), adj = 1, xpd = NA, cex = 0.8)
  if (!is.null(pooled)) {
    yp <- 1.1
    if (!is.na(pooled$pi_lo)) segments(pooled$pi_lo, yp, pooled$pi_hi, yp, col = COL$main, lwd = 2.5, lend = 1)
    polygon(c(pooled$lo, pooled$est, pooled$hi, pooled$est), yp + c(0, 0.38, 0, -0.38), col = COL$main, border = COL$main)
    text(xs[1], yp, pooled$label, adj = 0, font = 2, xpd = NA, cex = 0.82)
    if (!is.null(pooled$n)) text(xs[2], yp, pooled$n, adj = 1, font = 2, xpd = NA, cex = 0.8)
    if (!is.null(pooled$nevt)) text(xs[3], yp, pooled$nevt, adj = 1, font = 2, xpd = NA, cex = 0.8)
    text(xs[4], yp, fmt_hr(pooled$est, pooled$lo, pooled$hi), adj = 0, font = 2, xpd = NA, cex = 0.82)
    if (!is.null(df$w)) text(xs[5], yp, "100%", adj = 1, font = 2, xpd = NA, cex = 0.8)
    abline(h = 1.9, col = COL$light)
  }
  if (!is.null(het)) mtext(het, side = 1, line = 4.3, adj = 0, at = grconvertX(0.01, "nfc", "user"), cex = 0.72, col = "black")
  mtext(title, side = 3, line = 1.2, adj = 0, at = grconvertX(0.01, "nfc", "user"), font = 2, cex = 1.05)
  if (!is.null(favors)) mtext(paste(c("←", ""), favors, c("", "→")), side = 1, line = 3.2,
        at = c(rng[1], rng[2]), adj = c(0, 1), cex = 0.62, col = "black")
}

km_plot <- function(time, event, group, title = "", xlab = "Months", cols = GROUP_COLS, strata = NULL,
                    tmax = NULL, sub = NULL, ylab = "Survival probability", mark_t = NULL) {
  ok <- !is.na(time) & !is.na(event) & !is.na(group)
  time <- time[ok]; event <- event[ok]; group <- droplevels(as.factor(group[ok]))
  st <- if (!is.null(strata)) strata[ok]
  if (length(time) < 10 || nlevels(group) < 2) return(empty_plot("Not enough data for this analysis"))
  fit <- survfit(Surv(time, event) ~ group)
  if (is.null(tmax)) tmax <- min(120, quantile(time, 0.95))
  base_par(mar = c(8.5, 4.5, 3, 1))
  plot(fit, col = cols[seq_len(nlevels(group))], lwd = 2.2, xlim = c(0, tmax), xlab = xlab, ylab = ylab,
       mark.time = length(time) <= 800, pch = "|", cex = 0.6, axes = FALSE, conf.int = FALSE)
  axis(1, col = COL$grey); axis(2, col = COL$grey)
  if (!is.null(mark_t) && mark_t < tmax) abline(v = mark_t, lty = 3, col = COL$grey)
  lt <- if (is.null(st)) survdiff(Surv(time, event) ~ group) else survdiff(Surv(time, event) ~ group + strata(st))
  p <- pchisq(lt$chisq, length(lt$n) - 1, lower.tail = FALSE)
  legend("topright", legend = paste0(levels(group), " (n=", as.vector(table(group)), ")"),
         col = cols[seq_len(nlevels(group))], lwd = 2.2, bty = "n", cex = 0.85)
  mtext(title, side = 3, line = 1.2, adj = 0, font = 2)
  mtext(paste0(if (is.null(strata)) "Log-rank" else "Stratified log-rank", " p = ", fmt_p(p), if (!is.null(sub)) paste0("   |   ", sub)),
        side = 3, line = 0.1, adj = 0, cex = 0.8, col = "black")
  # numbers at risk
  tt <- pretty(c(0, tmax), 6); tt <- tt[tt <= tmax]
  s <- summary(fit, times = tt, extend = TRUE)
  stl <- if (is.null(s$strata)) rep(levels(group)[1], length(s$time)) else sub("group=", "", s$strata)
  mtext("At risk", side = 1, line = 3.4, adj = 1, at = grconvertX(0.06, "nfc", "user"), cex = 0.7, font = 2)
  for (i in seq_len(nlevels(group))) {
    g <- levels(group)[i]; n <- s$n.risk[stl == g]
    mtext(n, side = 1, line = 3.4 + i * 0.95, at = tt[seq_along(n)], cex = 0.68, col = cols[i])
  }
  invisible(p)
}

# =============================================================================
# UI
# =============================================================================
theme <- bs_theme(version = 5, bootswatch = "flatly", primary = COL$main, base_font = font_collection("system-ui", "-apple-system", "Segoe UI", "Roboto", "sans-serif"),
                  "navbar-bg" = COL$main, font_scale = 0.92)

# PNG / PDF buttons for a plot registered with reg_plot(id, ...)
dl_buttons <- function(id) div(class = "dl-bar",
  downloadButton(paste0(id, "_png"), "PNG", class = "btn-sm btn-dl"),
  downloadButton(paste0(id, "_pdf"), "PDF", class = "btn-sm btn-dl"))
dl_csv <- function(id, label) downloadButton(id, label, class = "btn-sm btn-dl mt-2")

stat_tile <- function(ic, out, label, sub) div(class = "stat-tile",
  div(class = "stat-icon", icon(ic)),
  div(div(class = "stat-value", textOutput(out, inline = TRUE)), div(class = "stat-label", label),
      div(class = "stat-sub", textOutput(sub, inline = TRUE))))

# scatter plot with a floating tooltip on hover and gene selection on click
hover_plot <- function(id, pfx, height = "560px") div(style = "position:relative",
  plotOutput(id, height = height, click = paste0(pfx, "_click"),
             hover = hoverOpts(paste0(pfx, "_hover"), delay = 50, delayType = "throttle", nullOutside = TRUE)),
  uiOutput(paste0(pfx, "_tip")))

CLIN_VARS <- c("Cell of origin (COO)" = "coo", "IPI risk group" = "ipi_risk", "IPI score (0-5)" = "ipi_score",
               "Ann Arbor stage" = "stage", "Advanced stage (III-IV)" = "adv_stage", "ECOG PS >= 2" = "ecog2",
               "LDH" = "ldh_high", "Age (continuous)" = "age", "Age > 60" = "age60", "Sex" = "sex",
               "Histology" = "histology", "Treatment group" = "tx_group", "Rituximab" = "rituximab",
               "Early failure on R-CHOP-based therapy (event < 24 mo)" = "ef_rchop", "Early failure, any treatment (event < 24 mo)" = "early_failure",
               "Death < 24 mo" = "death24")
TME_CELLS <- c("T cells", "CD8 T / cytotoxic", "Regulatory T cells", "NK cells", "Macrophages", "M2-like macrophages",
               "Dendritic cells", "Neutrophils", "IFN-gamma / inflammation", "Immune checkpoints", "MHC class II",
               "Matrix / CAF (stromal-1)", "Endothelium (stromal-2)", "Proliferation", "B-cell program (tumor)")
SUBGROUP_VARS <- c("coo", "ipi_risk", "adv_stage", "ecog2", "ldh_high", "age60", "sex", "rituximab")

ui <- page_navbar(
  title = tags$span(tags$b("DLBCL-PRISM"), tags$small(class = "ms-2 opacity-75 d-none",
                                                     "Prognostic Resource for Integrated Survival Meta-analysis")),
  theme = theme, fillable = FALSE, id = "nav",
  header = tags$head(tags$style(HTML("
    .card{border-radius:10px;border-color:#e3e7eb} .bslib-value-box .value-box-value{font-size:1.6rem}
    .small-note{font-size:.8rem;color:#57606a} .sidebar-title{font-weight:600}
    .dataTables_wrapper{font-size:.82rem} .nav-link{font-weight:500}
    .shiny-output-error-validation{color:#57606a;font-style:italic;padding:1rem}
    .plot-tip{position:absolute;pointer-events:none;z-index:20;background:rgba(255,255,255,.97);border:1px solid #d0d7de;
              border-radius:6px;padding:6px 9px;font-size:.78rem;line-height:1.35;box-shadow:0 3px 10px rgba(0,0,0,.12);white-space:nowrap}
    .dl-bar{display:flex;gap:.5rem;justify-content:flex-end;margin-top:.5rem}
    .btn-dl{background:#1f5f7a;color:#fff;border:1px solid #1f5f7a;font-weight:600;padding:.25rem .8rem}
    .btn-dl:hover,.btn-dl:focus{background:#174a60;color:#fff;border-color:#174a60}
    .stat-strip{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:1rem;margin-bottom:1rem}
    @media (max-width:900px){.stat-strip{grid-template-columns:repeat(2,minmax(0,1fr))}}
    .stat-tile{display:flex;align-items:center;gap:.9rem;background:#fff;border:1px solid #e3e7eb;border-radius:10px;
               padding:.85rem 1.1rem}
    .stat-icon{flex:0 0 44px;height:44px;border-radius:50%;background:#e8f1f5;color:#1f5f7a;display:flex;
               align-items:center;justify-content:center;font-size:1.15rem}
    .stat-value{font-size:1.55rem;font-weight:700;line-height:1.1;color:#000;font-variant-numeric:tabular-nums}
    .stat-label{font-size:.72rem;text-transform:uppercase;letter-spacing:.06em;color:#000;font-weight:600}
    .stat-sub{font-size:.72rem;color:#000}
  "))),
  sidebar = sidebar(width = 290, open = "desktop",
    selectizeInput("gene", "Gene", choices = NULL, options = list(placeholder = "Type a gene symbol...", maxOptions = 30)),
    radioButtons("endpoint", "Endpoint", c("Overall survival" = "OS", "Progression-free survival" = "PFS"), inline = FALSE),
    radioButtons("model", "Cox model", c("Univariable" = "uni", "Multivariable (adjusted for age + IPI)" = "multi"))
  ),

  # ---------------------------------------------------------------- Overview
  nav_panel("Overview", icon = icon("house"),
    div(class = "stat-strip",
      stat_tile("layer-group", "vb_coh", "Cohorts", "vb_coh_sub"),
      stat_tile("users", "vb_pat", "Patients", "vb_pat_sub"),
      stat_tile("heart-pulse", "vb_evt", "OS events", "vb_evt_sub"),
      stat_tile("dna", "vb_gen", "Genes meta-analyzed", "vb_gen_sub")),
    layout_columns(col_widths = c(7, 5),
      card(card_header("Cohort catalog (non-redundant DLBCL cohorts)"), DTOutput("cohort_tbl")),
      card(card_header("Clinical data availability (% of patients)"), plotOutput("avail_plot", height = "560px"), dl_buttons("avail_plot")))),

  # ---------------------------------------------------------------- Forest
  nav_panel("Meta forest", icon = icon("chart-gantt"),
    layout_columns(col_widths = c(9, 3),
      card(full_screen = TRUE, card_header(textOutput("forest_title", inline = TRUE)),
           plotOutput("forest", height = "620px"), dl_buttons("forest")),
      card(card_header("Options"),
           div(class = "small-note mb-2", "Random-effects meta-analysis: REML τ² with Hartung-Knapp-Sidik-Jonkman confidence interval."),
           selectInput("sort_by", "Sort cohorts by", c("Weight" = "w", "Hazard ratio" = "est", "Name" = "label", "Sample size" = "n")),
           selectizeInput("excl", "Exclude cohorts", choices = NULL, multiple = TRUE),
           hr(), uiOutput("forest_summary"),
           dl_csv("dl_forest_csv", "Cohort estimates (CSV)"))),
    card(card_header("Leave-one-cohort-out sensitivity"), plotOutput("loo", height = "560px"), dl_buttons("loo"))
  ),

  # ---------------------------------------------------------------- Discovery
  nav_panel("Discovery", icon = icon("magnifying-glass-chart"),
    layout_columns(col_widths = c(3, 9),
      card(card_header("Filters"),
           sliderInput("kmin", "Minimum number of cohorts (k)", 2, 25, 8, step = 1),
           selectInput("fdr", "FDR threshold", c("0.001", "0.01", "0.05", "0.10"), "0.05"),
           sliderInput("hrmin", "Minimum effect |HR - 1| per SD", 0, 0.4, 0.05, step = 0.01),
           sliderInput("i2max", "Maximum I² (%)", 0, 100, 50, step = 5),
           sliderInput("consmin", "Minimum direction consistency (%)", 50, 100, 70, step = 5),
           radioButtons("dir", "Direction", c("Both" = "both", "Adverse (HR > 1)" = "risk", "Favorable (HR < 1)" = "prot"), inline = FALSE),
           hr(), uiOutput("disc_n"),
           dl_csv("dl_hits", "Download hits (CSV)")),
      navset_card_tab(full_screen = TRUE,
        nav_panel("Volcano", hover_plot("volcano", "volc"),
                  div(class = "small-note", "Hover over a point to see the gene; click to select it (all tabs update)."), dl_buttons("volcano")),
        nav_panel("Top hits forest", sliderInput("ntop", "Number of genes", 5, 40, 20, width = "300px"),
                  plotOutput("top_forest", height = "720px"), dl_buttons("top_forest")),
        nav_panel("Table", DTOutput("hits_tbl"))))
  ),

  # ---------------------------------------------------------------- Clinical
  nav_panel("Clinical associations", icon = icon("stethoscope"),
    layout_columns(col_widths = c(3, 9),
      card(card_header("Variable"),
           selectInput("cvar", NULL, CLIN_VARS, "coo"),
           checkboxGroupInput("cplot", "Display", c("Violin" = "violin", "Points" = "pts"), c("violin", "pts"), inline = TRUE),
           div(class = "small-note", "Expression is the within-cohort z-score, so values are comparable across platforms. Unadjusted tests: Wilcoxon / Kruskal-Wallis / Spearman. Cohort-adjusted: linear model z ~ variable + cohort (F-test).")),
      card(full_screen = TRUE, card_header(textOutput("clin_title", inline = TRUE)),
           layout_columns(col_widths = c(7, 5),
             plotOutput("clin_plot", height = "460px"),
             plotOutput("clin_smd", height = "460px")),
           dl_buttons("clin_card"))),
    card(card_header("All clinical variables at a glance"), DTOutput("clin_all"))
  ),

  # ---------------------------------------------------------------- Survival
  nav_panel("Survival", icon = icon("chart-line"),
    layout_columns(col_widths = c(3, 9),
      card(card_header("Options"),
           selectInput("cut", "Expression groups (within cohort)", c("Median" = "median", "Tertiles" = "tertile", "Quartiles" = "quartile")),
           selectizeInput("km_coh", "Cohort", choices = c("All (stratified by cohort)" = "ALL")),
           selectizeInput("km_excl", "Exclude cohorts", choices = NULL, multiple = TRUE,
                          options = list(placeholder = "None (applies to KM, Cox and subgroups)")),
           checkboxGroupInput("adj", "Adjust Cox model for", c("IPI risk" = "ipi_risk", "COO" = "coo", "Age" = "age", "Sex" = "sex", "Stage III-IV" = "adv_stage"), inline = FALSE),
           div(class = "small-note", "Pooled analyses use Cox models stratified by cohort (each cohort keeps its own baseline hazard). Adjusted models use complete cases.")),
      card(full_screen = TRUE, card_header(textOutput("km_title", inline = TRUE)),
           plotOutput("km", height = "520px"), tableOutput("cox_tbl"), dl_buttons("km"))),
    card(full_screen = TRUE, card_header("Subgroup analysis: gene HR (per SD) within clinical subgroups, with interaction tests"),
         plotOutput("subgroup", height = "640px"), dl_buttons("subgroup"))
  ),

  # ---------------------------------------------------------------- Treatment
  nav_panel("Treatment & response", icon = icon("syringe"),
    div(class = "alert alert-light py-2 small",
        tags$b("Early treatment failure on R-CHOP-based therapy"), " = progression, relapse or death (PFS, or EFS when PFS is unavailable) within 24 months of diagnosis in patients treated with first-line R-CHOP (including R-CHOP + bortezomib, REMoDL-B; R-CHOP assigned by study design for GSE31312 and GSE21846). Patients censored before 24 months are excluded."),
    card(full_screen = TRUE, card_header(textOutput("ef_title", inline = TRUE)),
         layout_columns(col_widths = c(4, 8),
           plotOutput("ef_box", height = "460px"),
           plotOutput("ef_forest", height = "460px")),
         uiOutput("ef_stats"), dl_buttons("ef_card")),
    layout_columns(col_widths = c(3, 9),
      card(card_header("Kaplan-Meier options"),
           selectInput("tx_reg", "Regimen", c("R-CHOP-based (all)" = "all", "R-CHOP" = "R-chemo", "R-CHOP + bortezomib (RB-CHOP)" = "R-chemo + bortezomib")),
           radioButtons("tx_ep", "Endpoint", c("Progression-free survival" = "PFS", "Overall survival" = "OS")),
           selectInput("tx_cut", "Expression groups (within cohort)", c("Median" = "median", "Tertiles" = "tertile", "Quartiles" = "quartile")),
           selectInput("tx_coo", "Cell of origin", c("All" = "all", "ABC" = "ABC", "GCB" = "GCB", "Unclassified" = "Unclassified")),
           selectInput("tx_ipi", "IPI risk", c("All" = "all", "Low (0-1)" = "Low (0-1)", "Intermediate (2-3)" = "Intermediate (2-3)", "High (4-5)" = "High (4-5)")),
           selectInput("tx_tmax", "Follow-up shown", c("24 months" = 24, "60 months" = 60, "120 months" = 120), 60),
           div(class = "small-note", "Groups are formed within each cohort; log-rank and Cox models are stratified by cohort. The dotted line marks 24 months.")),
      card(full_screen = TRUE, card_header(textOutput("tx_km_title", inline = TRUE)),
           plotOutput("tx_km", height = "520px"), dl_buttons("tx_km"),
           tableOutput("tx_cox"), tableOutput("tx_rates")))
  ),

  # ---------------------------------------------------------------- Pathways
  nav_panel("Pathways", icon = icon("diagram-project"),
    layout_columns(col_widths = c(3, 9),
      card(card_header("Pre-ranked GSEA"),
           radioButtons("gs_src", "Ranking", c("Survival (follows sidebar)" = "surv", "Early failure on R-CHOP-based therapy" = "ef")),
           radioButtons("gs_coll", "Gene-set collection", c("Hallmark" = "Hallmark", "Reactome" = "Reactome"), inline = TRUE),
           selectInput("gs_fdr", "FDR threshold", c("0.01", "0.05", "0.10", "0.25"), "0.05"),
           hr(), uiOutput("gs_n"),
           dl_csv("dl_gsea", "Download results (CSV)"),
           div(class = "small-note mt-2", "Genes are ranked by the pooled meta-analysis statistic (log HR / SE; log OR / SE for early failure) among genes measured in k >= 5 cohorts (k >= 4 for early failure). fgsea on MSigDB gene sets of 15-500 genes. Positive NES = pathway genes associated with worse outcome.")),
      navset_card_tab(id = "gs_tabs", full_screen = TRUE,
        nav_panel("Top pathways", plotOutput("gs_bar", height = "700px"), dl_buttons("gs_bar")),
        nav_panel("Enrichment plot",
                  selectizeInput("gs_path", "Pathway", choices = NULL, width = "100%"),
                  plotOutput("gs_enrich", height = "480px"), dl_buttons("gs_enrich"), uiOutput("gs_le")),
        nav_panel("Table", DTOutput("gs_tbl"),
                  div(class = "small-note", "Click a row to open its enrichment plot."))))
  ),

  # ---------------------------------------------------------------- Microenvironment
  nav_panel("Microenvironment", icon = icon("cubes"),
    div(class = "alert alert-light py-2 small", tags$b("Microenvironment scores"),
        " = mean within-cohort z-score of canonical marker genes for each cell type or program, re-standardized within cohort and computed when >= 50% of the markers are measured. They reflect relative abundance, not absolute cell fractions; the B-cell program and MHC class II largely reflect the tumor cells."),
    layout_columns(col_widths = c(7, 5),
      card(full_screen = TRUE, card_header(textOutput("tme_prog_title", inline = TRUE)),
           plotOutput("tme_prog", height = "560px"), dl_buttons("tme_prog")),
      card(full_screen = TRUE, card_header(textOutput("tme_cor_title", inline = TRUE)),
           plotOutput("tme_cor", height = "560px"),
           div(class = "small-note", "Pearson correlation of within-cohort z-scores in each cohort, Fisher-transformed and pooled by REML-HKSJ. * = the gene is one of the markers of that score, so the correlation is partly circular."),
           dl_buttons("tme_cor"))),
    layout_columns(col_widths = c(3, 9),
      card(card_header("Cell type / program"),
           selectInput("tme_cell", NULL, TME_CELLS, "M2-like macrophages"),
           uiOutput("tme_markers"), hr(),
           tags$b("Kaplan-Meier", class = "small"),
           selectInput("tme_cut", "Score groups (within cohort)", c("Median" = "median", "Tertiles" = "tertile", "Quartiles" = "quartile")),
           selectInput("tme_pop", "Patients", c("All curated patients" = "all", "R-CHOP-based therapy only" = "rchop")),
           selectInput("tme_coo_f", "Cell of origin", c("All" = "all", "ABC" = "ABC", "GCB" = "GCB")),
           div(class = "small-note", "Endpoint and Cox model follow the sidebar.")),
      navset_card_tab(full_screen = TRUE, title = textOutput("tme_sel_title", inline = TRUE),
        nav_panel("Associations",
          layout_columns(col_widths = c(6, 6),
            plotOutput("tme_scatter", height = "420px"),
            plotOutput("tme_coo", height = "420px")),
          dl_buttons("tme_assoc")),
        nav_panel("Treatment response",
          div(class = "small-note mb-2", "Early failure = progression, relapse or death within 24 months in patients treated with R-CHOP-based therapy (patients censored before 24 months excluded)."),
          layout_columns(col_widths = c(5, 7),
            plotOutput("tme_ef", height = "420px"),
            div(class = "pt-4", uiOutput("tme_ef_stats"))),
          plotOutput("tme_ef_forest", height = "480px"), dl_buttons("tme_txr")),
        nav_panel("Kaplan-Meier", plotOutput("tme_km", height = "520px"), dl_buttons("tme_km"), tableOutput("tme_km_cox")),
        nav_panel("Forest by cohort", plotOutput("tme_forest", height = "620px"), dl_buttons("tme_forest"))))
  ),

  # ---------------------------------------------------------------- Signature
  nav_panel("Signature", icon = icon("layer-group"),
    layout_columns(col_widths = c(3, 9),
      card(card_header("Build a multi-gene score"),
           textAreaInput("sig_genes", "Genes (one per line or comma-separated, max 60)", rows = 7, placeholder = "e.g. ACP5\nFMOD\nMS4A4A"),
           actionButton("sig_top", "Fill with top Discovery hits", class = "btn-sm btn-outline-primary mb-2"),
           numericInput("sig_n", "Number of top hits", 15, 2, 60),
           radioButtons("sig_w", "Weights", c("Sign of pooled log-HR (±1)" = "sign", "Pooled log-HR" = "beta")),
           checkboxInput("sig_loco", "Leave-one-cohort-out weights (recommended)", TRUE),
           actionButton("sig_run", "Run", class = "btn-primary w-100"),
           div(class = "small-note mt-2", "Score = weighted mean of available within-cohort z-scores (>= 50% of genes required). With leave-one-cohort-out, the weights applied to each cohort are re-estimated from a meta-analysis that excludes that cohort, giving an honest internal-external validation.")),
      navset_card_tab(full_screen = TRUE,
        nav_panel("Validation forest", plotOutput("sig_forest", height = "600px"), dl_buttons("sig_forest")),
        nav_panel("Kaplan-Meier", plotOutput("sig_km", height = "520px"), dl_buttons("sig_km")),
        nav_panel("C-index vs IPI", plotOutput("sig_cidx", height = "460px"), dl_buttons("sig_cidx"), tableOutput("sig_cidx_tbl")),
        nav_panel("Genes used", tableOutput("sig_used"))))
  ),

  # ---------------------------------------------------------------- About
  nav_panel("Methods", icon = icon("book"),
    card(card_body(markdown("
### DLBCL-PRISM: methods

**Data.** Patient-level expression profiles and harmonized clinical annotation from publicly available DLBCL cohorts (GEO, ArrayExpress, NCI GDC). Expression was standardized **within each cohort** (gene-wise z-score; non-log data log2-transformed first), making effects interpretable per 1 SD across platforms (Affymetrix, Illumina, RNA-seq, NanoString). Survival times were harmonized to months.

**Cohort-level models.** In each cohort, Cox proportional hazards models of OS and PFS on gene z-score: univariable, and multivariable adjusted for age and IPI (high vs low risk); in cohorts without IPI the multivariable model is adjusted for age only.

**Meta-analysis.** Random-effects inverse-variance pooling with REML estimation of τ² and Hartung-Knapp-Sidik-Jonkman (HKSJ) confidence intervals and p-values (t distribution, k-1 df; untruncated, as metafor test = knha). I², Cochran's Q and leave-one-cohort-out sensitivity are reported. Benjamini-Hochberg FDR across genes. Direction consistency = % of cohorts whose estimate shares the sign of the pooled estimate.

**Early treatment failure on R-CHOP-based therapy.** Event (PFS, or EFS when PFS unavailable) within 24 months among patients treated with first-line R-CHOP, including R-CHOP + bortezomib (REMoDL-B); GSE31312 (International DLBCL R-CHOP Consortium) and GSE21846 were assigned R-CHOP by study design; ibrutinib-treated GSE93984 is not included; excluding patients censored before 24 months. Per-cohort logistic regression on gene z-score (cohorts with >= 5 failures and >= 5 non-failures), pooled by REML-HKSJ.

**Non-redundant cohort set.** Patient overlap between public series was detected by matching survival time, status and age: GSE87371 and GSE136971 contain the same 221 patients, and GSE11318 and GSE53786 are largely contained in GSE10846. All analyses keep one copy (GSE136971, GSE10846); patient-level analyses additionally exclude BL and PMBL samples.

**Kaplan-Meier.** Expression (or microenvironment score) groups are formed by quantiles within each cohort; log-rank tests and Cox models are stratified by cohort. In the Treatment tab, analyses are restricted to patients treated with R-CHOP-based therapy, optionally by regimen, cell of origin and IPI risk group; 24- and 60-month estimates are pooled Kaplan-Meier estimates.

**Patient-level analyses.** Cox models stratified by cohort; Kaplan-Meier by within-cohort quantiles; subgroup interaction tests by likelihood-ratio tests of product terms; clinical associations by non-parametric tests, cohort-adjusted linear models and per-cohort Hedges' g pooled by REML-HKSJ.

**Pathways.** Pre-ranked gene-set enrichment analysis (fgsea, multilevel) of the pooled meta-analysis statistic (log HR / SE, or log OR / SE for early failure) among genes measured in at least 5 cohorts (4 for early failure), using MSigDB Hallmark and Reactome gene sets of 15-500 genes; FDR within each collection.

**Microenvironment.** For each patient, cell-type and program scores were computed as the mean within-cohort z-score of canonical marker genes (T cells, CD8/cytotoxic, regulatory T, NK, macrophages, M2-like macrophages, dendritic cells, neutrophils, IFN-γ/inflammation, immune checkpoints, MHC class II, matrix/CAF, endothelium, proliferation, B-cell program), re-standardized within cohort and computed when at least half of the markers were measured in that cohort. Scores reflect relative abundance, not absolute fractions. Per-cohort Cox models of each score were pooled by REML-HKSJ; gene-score associations are per-cohort Pearson correlations pooled after Fisher transformation.

**Signature.** Weighted mean z-score; leave-one-cohort-out re-estimation of weights; Harrell's C-index (stratified by cohort) compared with the IPI.

*Runs entirely in your browser (Shinylive / webR). No data leave your computer.*
"))))
)

# =============================================================================
# Server
# =============================================================================
server <- function(input, output, session) {

  base_url <- reactive({
    if (nzchar(DATA_URL)) return(DATA_URL)
    cd <- session$clientData
    path <- sub("[^/]*$", "", cd$url_pathname)
    path <- sub("app_[^/]*/$", "", path)          # shinylive serves the app under app_<id>/
    paste0(cd$url_protocol, "//", cd$url_hostname, if (nzchar(cd$url_port)) paste0(":", cd$url_port), path, "data/")
  })
  F <- function(rel) fetch_rds(rel, base_url())

  # PNG (300 dpi) / PDF download handlers <id>_png, <id>_pdf; w and h in inches (number or function)
  reg_dl <- function(id, draw, w = 10, h = 7, name = function() id) {
    for (f in c("png", "pdf")) local({
      ff <- f
      output[[paste0(id, "_", ff)]] <- downloadHandler(
        filename = function() sprintf("DLBCL-PRISM_%s.%s", gsub("[^A-Za-z0-9_.-]+", "_", name()), ff),
        content = function(file) {
          wi <- if (is.function(w)) w() else w; hi <- if (is.function(h)) h() else h
          if (ff == "png") png(file, wi, hi, units = "in", res = 300) else pdf(file, wi, hi)
          on.exit(dev.off()); draw()
        })
    })
  }
  reg_plot <- function(id, draw, ...) { output[[id]] <- renderPlot(draw(), res = 96); reg_dl(id, draw, ...) }
  safe_draw <- function(f) tryCatch(f(), error = function(e) empty_plot(conditionMessage(e)))
  model_txt <- function() if (input$model == "uni") "univariable" else "multivariable"

  D <- reactive({
    withProgress(message = "Loading cohorts...", value = 0.3, {
      clin <- F("clinical.rds")
      clin$ef_rchop <- ifelse(clin$tx_group %in% c("R-chemo", "R-chemo + bortezomib"), clin$early_failure, NA)
      gidx <- F("gene_index.rds"); gidx <- gidx[clean_key(gidx$key), ]
      list(clin = clin, gidx = gidx, coh = F("cohorts.rds"), info = F("meta_info.rds"))
    })
  })

  observeEvent(D(), {
    g <- D()$gidx$key
    pooled_os <- F("pooled_OS_uni_curated.rds")
    first <-pooled_os$key[pooled_os$k >= 8][which.min(pooled_os$p[pooled_os$k >= 8])]
    updateSelectizeInput(session, "gene", choices = g, selected = first, server = TRUE)
    updateSelectizeInput(session, "km_coh", choices = c("All (stratified by cohort)" = "ALL", sort(unique(D()$clin$cohort))))
  }, once = TRUE)

  gene <- reactive({ req(input$gene, input$gene %in% D()$gidx$key); input$gene })


  # ---- sample filter (patient-level): non-redundant DLBCL cohorts, BL / PMBL samples removed ----
  excluded_cohorts <- reactive(c(names(D()$info$overlap), D()$info$non_dlbcl))
  keep <- reactive({ cl <- D()$clin; !(cl$cohort %in% excluded_cohorts()) & !(cl$histology %in% c("BL", "PMBL")) })

  # ---- gene data -----------------------------------------------------------------
  gdat <- reactive({
    g <- gene(); gi <- D()$gidx[D()$gidx$key == g, ]
    ch <- withProgress(message = paste("Loading", g, "..."), value = 0.5, F(sprintf("chunks/c%04d.rds", gi$chunk)))
    z <- if (g %in% rownames(ch$expr)) decode_z(ch$expr[g, ]) else rep(NA_real_, nrow(D()$clin))
    list(z = z, res = ch$res[[g]])
  })
  pdat <- reactive({                         # patient-level frame for current gene & filter
    d <- cbind(D()$clin, z = gdat()$z)[keep(), ]
    d[!is.na(d$z), ]
  })
  surv_cols <- reactive(if (input$endpoint == "OS") c("os_m", "os_e") else c("pfs_m", "pfs_e"))

  pooled_raw <- reactive(F(sprintf("pooled_%s_%s_curated.rds", input$endpoint, input$model)))   # ranking used by fgsea
  pooled <- reactive({ p <- pooled_raw(); p[clean_key(p$key), ] })

  # ============================== Overview =======================================
  cl_keep <- reactive(D()$clin[keep(), ])
  fmt_n <- function(x) format(x, big.mark = ",")
  output$vb_coh <- renderText(length(unique(cl_keep()$cohort)))
  output$vb_coh_sub <- renderText(sprintf("%d platforms", length(unique(D()$coh$platform[D()$coh$cohort %in% cl_keep()$cohort]))))
  output$vb_pat <- renderText(fmt_n(nrow(cl_keep())))
  output$vb_pat_sub <- renderText(sprintf("%s with OS data", fmt_n(sum(!is.na(cl_keep()$os_e)))))
  output$vb_evt <- renderText(fmt_n(sum(cl_keep()$os_e, na.rm = TRUE)))
  output$vb_evt_sub <- renderText(sprintf("%s PFS/EFS events", fmt_n(sum(cl_keep()$pfs_e, na.rm = TRUE))))
  output$vb_gen <- renderText(fmt_n(sum(pooled()$k >= 2)))
  output$vb_gen_sub <- renderText(sprintf("%s in >= 10 cohorts", fmt_n(sum(pooled()$k >= 10))))

  coh_keep <- reactive({                     # catalog without the redundant / non-DLBCL cohorts
    ct <- D()$coh[!D()$coh$cohort %in% excluded_cohorts(), ]
    ct$n_dlbcl <- as.vector(table(factor(cl_keep()$cohort, ct$cohort)))
    ct
  })
  output$cohort_tbl <- renderDT({
    ct <- coh_keep()[, c("cohort", "platform", "n_dlbcl", "genes", "os_n", "os_events", "pfs_n", "pfs_events", "pfs_type")]
    names(ct) <- c("Cohort", "Platform", "N (DLBCL)", "Genes", "OS n", "OS events", "PFS n", "PFS events", "PFS/EFS")
    datatable(ct, rownames = FALSE, options = list(pageLength = 30, dom = "ft", scrollY = "520px", scrollCollapse = TRUE),
              selection = "none")
  })

  reg_plot("avail_plot", function() {
    ct <- coh_keep(); v <- c("os_n", "pfs_n", "coo", "ipi", "stage", "ecog", "ldh", "age", "sex", "treatment")
    m <- as.matrix(ct[, v]) / ct$n * 100; rownames(m) <- ct$cohort
    m <- m[order(-ct$n), ]
    labs <- c("OS", "PFS/EFS", "COO", "IPI", "Stage", "ECOG", "LDH", "Age", "Sex", "Tx")
    base_par(mar = c(4, 9.5, 1, 1))
    pal <- colorRampPalette(c("#f6f8fa", "#9ecae1", COL$main))(101)
    image(seq_along(v), seq_len(nrow(m)), t(m[nrow(m):1, ]), col = pal, zlim = c(0, 100), axes = FALSE, xlab = "", ylab = "")
    axis(1, at = seq_along(v), labels = labs, las = 2, tick = FALSE, cex.axis = 0.8)
    axis(2, at = seq_len(nrow(m)), labels = rev(rownames(m)), tick = FALSE, cex.axis = 0.72)
    abline(h = seq_len(nrow(m)) + 0.5, v = seq_along(v) + 0.5, col = "white", lwd = 1.5)
  }, 6, 8, function() "clinical_data_availability")

  # ============================== Meta forest ====================================
  observe({
    r <- gdat()$res; req(r)
    cohs <- sort(unique(r$cohort))
    updateSelectizeInput(session, "excl", choices = cohs, selected = intersect(isolate(input$excl), cohs))
  })

  forest_data <- reactive({
    r <- gdat()$res
    validate(need(!is.null(r), "No cohort-level estimates for this gene."))
    r <- r[r$endpoint == input$endpoint & r$model == input$model & !r$cohort %in% c(excluded_cohorts(), input$excl), ]
    validate(need(nrow(r) >= 1, "No cohort-level estimates for this gene / endpoint / model."))
    m <- reml(r$logHR, r$se^2)
    df <- data.frame(label = r$cohort, n = r$n, nevt = r$nevt, est = r$logHR,
                     lo = r$logHR - 1.96 * r$se, hi = r$logHR + 1.96 * r$se, w = m$w, se = r$se)
    o <- switch(input$sort_by, w = order(-df$w), est = order(-df$est), n = order(-df$n), label = order(df$label))
    list(df = df[o, ], m = m)
  })

  forest_title_txt <- reactive(sprintf("%s: %s, %s", gene(), input$endpoint, if (input$model == "uni") "univariable" else "multivariable"))
  output$forest_title <- renderText(forest_title_txt())

  draw_forest <- function() {
    fd <- forest_data(); m <- fd$m; df <- fd$df
    pooled <- list(est = m$est, lo = m$lo, hi = m$hi, label = "Pooled (REML-HKSJ)",
                   n = sum(df$n), nevt = sum(df$nevt),
                   pi_lo = NA, pi_hi = NA)
    het <- if (m$k > 1) sprintf("Heterogeneity: I² = %.0f%%, τ² = %.4f, Q = %.1f (p = %s)   |   Overall effect: p = %s",
                                m$I2, m$tau2, m$Q, fmt_p(m$pQ), fmt_p(m$p)) else NULL
    forest_plot(df, pooled, title = forest_title_txt(), het = het)
  }
  reg_plot("forest", draw_forest, 11, function() 2.2 + 0.3 * nrow(forest_data()$df),
           function() paste("forest", gene(), input$endpoint, input$model, sep = "_"))

  output$forest_summary <- renderUI({
    m <- forest_data()$m
    tagList(
      tags$div(tags$b("Pooled HR: "), sprintf("%.2f (%.2f-%.2f)", exp(m$est), exp(m$lo), exp(m$hi))),
      tags$div(tags$b("p: "), fmt_p(m$p), " | ", tags$b("k: "), m$k),
      if (m$k > 1) tags$div(tags$b("I²: "), sprintf("%.0f%%", m$I2), " | ", tags$b("τ²: "), sprintf("%.4f", m$tau2))
    )
  })

  reg_plot("loo", function() {
    df <- forest_data()$df; validate(need(nrow(df) >= 3, "Need at least 3 cohorts."))
    l <- lapply(seq_len(nrow(df)), function(i) { m <- reml(df$est[-i], df$se[-i]^2); c(m$est, m$lo, m$hi, m$I2) })
    l <- as.data.frame(do.call(rbind, l)); names(l) <- c("est", "lo", "hi", "I2")
    l$label <- paste("without", df$label); l$n <- NULL
    m <- forest_data()$m
    forest_plot(data.frame(label = l$label, est = l$est, lo = l$lo, hi = l$hi),
                list(est = m$est, lo = m$lo, hi = m$hi, label = "All cohorts", pi_lo = NA, pi_hi = NA),
                cols = "Omitted cohort", title = sprintf("%s: leave-one-cohort-out (%s, %s)", gene(), input$endpoint, model_txt()))
  }, 11, function() 2.2 + 0.3 * nrow(forest_data()$df), function() paste("leave_one_out", gene(), input$endpoint, input$model, sep = "_"))
  output$dl_forest_csv <- downloadHandler(
    filename = function() sprintf("cohort_estimates_%s_%s_%s.csv", gene(), input$endpoint, input$model),
    content = function(file) {
      df <- forest_data()$df
      write.csv(data.frame(cohort = df$label, n = df$n, events = df$nevt, logHR = df$est, se = df$se,
                           HR = exp(df$est), lo = exp(df$lo), hi = exp(df$hi), weight = df$w), file, row.names = FALSE)
    })

  # ============================== Discovery ======================================
  hits_all <- reactive({
    p <- add_ci_fdr(pooled(), input$kmin)
    p$pass <- p$k >= input$kmin & p$fdr < as.numeric(input$fdr) & abs(p$HR - 1) >= input$hrmin &
              (is.na(p$I2) | p$I2 <= input$i2max) & p$cons * 100 >= input$consmin &
              (input$dir == "both" | (input$dir == "risk" & p$logHR > 0) | (input$dir == "prot" & p$logHR < 0))
    p
  })
  hits <- reactive({ h <- hits_all(); h <- h[h$pass, ]; h[order(h$p), ] })

  output$disc_n <- renderUI({
    h <- hits(); tagList(tags$h4(nrow(h), class = "mb-0"),
                         div(class = "small-note", sprintf("genes pass (%d adverse, %d favorable) out of %s with k >= %d",
                             sum(h$logHR > 0), sum(h$logHR < 0), format(sum(hits_all()$k >= input$kmin), big.mark = ","), input$kmin)))
  })

  volc_df <- reactive({ h <- hits_all(); h <- h[h$k >= input$kmin, ]; h$y <- -log10(pmax(h$p, 1e-300)); h })
  reg_plot("volcano", function() {
    h <- volc_df(); validate(need(nrow(h) > 0, "No genes with this minimum number of cohorts."))
    base_par(mar = c(4.5, 4.5, 2.5, 1))
    col <- ifelse(!h$pass, adjustcolor(COL$grey, 0.25), ifelse(h$logHR > 0, adjustcolor(COL$risk, 0.8), adjustcolor(COL$prot, 0.8)))
    o <- order(h$pass)
    plot(h$logHR[o], h$y[o], ylim = c(0, max(h$y) * 1.1), pch = 19, cex = ifelse(h$pass[o], 0.75, 0.45), col = col[o], bty = "l", fg = COL$grey,
         xlab = "Pooled log HR per SD", ylab = expression(-log[10]~italic(p)))
    abline(v = 0, col = COL$grey, lty = 3)
    if (any(h$fdr < as.numeric(input$fdr))) abline(h = min(h$y[h$fdr < as.numeric(input$fdr)]), lty = 2, col = COL$grey)
    top <- head(h[h$pass, ][order(h$p[h$pass]), ], 15)
    if (nrow(top)) text(top$logHR, top$y, top$key, pos = 3, cex = 0.7, offset = 0.35, col = COL$ink, xpd = NA)
    g <- isolate(input$gene); if (!is.null(g) && g %in% h$key) {
      s <- h[h$key == g, ]; points(s$logHR, s$y, pch = 21, cex = 1.3, lwd = 1.2, col = "white", bg = COL$main)
      text(s$logHR, s$y, g, pos = 4, font = 2, cex = 0.75, col = COL$main, offset = 0.4)
    }
    mtext(sprintf("%s, %s, k >= %d", input$endpoint, model_txt(), input$kmin),
          side = 3, line = 0.6, adj = 0, cex = 0.9, col = "black")
  }, 9, 7, function() paste("volcano", input$endpoint, input$model, "k", input$kmin, sep = "_"))

  near_gene <- function(ev, df, x, y) { if (is.null(ev)) return(NULL); np <- nearPoints(df, ev, xvar = x, yvar = y, maxpoints = 1, threshold = 8); if (nrow(np)) np$key else NULL }
  set_gene <- function(g) if (!is.null(g)) updateSelectizeInput(session, "gene", choices = D()$gidx$key, selected = g, server = TRUE)
  observeEvent(input$volc_click, set_gene(near_gene(input$volc_click, volc_df(), "logHR", "y")))
  # floating tooltip next to the cursor (flips left on the right-hand side of the plot)
  plot_tip <- function(ev, html) {
    if (is.null(ev) || is.null(html)) return(NULL)
    right <- ev$coords_img$x > 0.65 * ev$range$right
    style <- sprintf("left:%.0fpx;top:%.0fpx;%s", ev$coords_css$x + if (right) -14 else 14, ev$coords_css$y + 12,
                     if (right) "transform:translateX(-100%);" else "")
    div(class = "plot-tip", style = style, HTML(html))
  }
  gene_tip <- function(s) sprintf("<b>%s</b><br>HR %.2f (%.2f-%.2f)<br>p = %s &nbsp; FDR = %s<br>k = %d &nbsp; I² = %.0f%% &nbsp; consistency %.0f%%",
                                  s$key, s$HR, s$lo, s$hi, fmt_p(s$p), fmt_p(s$fdr), s$k, ifelse(is.na(s$I2), 0, s$I2), 100 * s$cons)
  output$volc_tip <- renderUI({
    g <- near_gene(input$volc_hover, volc_df(), "logHR", "y"); if (is.null(g)) return(NULL)
    plot_tip(input$volc_hover, gene_tip(volc_df()[volc_df()$key == g, ]))
  })

  draw_top_forest <- function() {
    h <- head(hits(), input$ntop); validate(need(nrow(h) > 0, "No genes pass the filters."))
    df <- data.frame(label = h$key, n = h$k, nevt = h$npat, est = h$logHR, lo = log(h$lo), hi = log(h$hi))
    forest_plot(df, NULL, title = sprintf("Top %d prognostic genes (%s, %s)", nrow(h), input$endpoint, if (input$model == "uni") "univariable" else "multivariable"),
                cols = c("Gene", "k", "N"))
  }
  reg_plot("top_forest", draw_top_forest, 11, function() 2.2 + 0.3 * min(nrow(hits()), input$ntop),
           function() paste("top_hits", input$endpoint, input$model, sep = "_"))

  hits_view <- reactive({
    h <- hits(); data.frame(Gene = h$key, k = h$k, Patients = h$npat, HR = round(h$HR, 3), `95% CI` = sprintf("%.2f-%.2f", h$lo, h$hi),
                            p = signif(h$p, 3), FDR = signif(h$fdr, 3), `I2 (%)` = round(h$I2, 1), `Consistency (%)` = round(100 * h$cons),
                            check.names = FALSE)
  })
  output$hits_tbl <- renderDT(datatable(hits_view(), rownames = FALSE, selection = "single",
                                        options = list(pageLength = 15, scrollX = TRUE)))
  observeEvent(input$hits_tbl_rows_selected, set_gene(hits_view()$Gene[input$hits_tbl_rows_selected]))
  output$dl_hits <- downloadHandler(filename = function() sprintf("DLBCL-PRISM_hits_%s_%s.csv", input$endpoint, input$model),
                                    content = function(file) write.csv(hits_view(), file, row.names = FALSE))

  # ============================== Clinical =======================================
  clin_df <- reactive({
    d <- pdat(); v <- input$cvar; d$x <- d[[v]]
    d <- d[!is.na(d$x), ]
    validate(need(nrow(d) >= 10, "Not enough patients with this variable and expression of the gene."))
    if (is.character(d$x)) d$x <- factor(d$x)
    if (is.factor(d$x)) d$x <- droplevels(d$x)
    d
  })
  output$clin_title <- renderText(sprintf("%s z-score by %s", gene(), names(CLIN_VARS)[CLIN_VARS == input$cvar]))

  clin_test <- function(d) {
    if (is.numeric(d$x)) {
      ct <- suppressWarnings(cor.test(d$z, d$x, method = "spearman"))
      fa <- anova(lm(z ~ x + cohort, d))
      list(type = "num", txt = sprintf("Spearman rho = %.3f, p = %s  |  cohort-adjusted (lm) p = %s  |  n = %d, %d cohorts",
                                       ct$estimate, fmt_p(ct$p.value), fmt_p(fa["x", "Pr(>F)"]), nrow(d), length(unique(d$cohort))),
           p = ct$p.value, padj = fa["x", "Pr(>F)"], eff = ct$estimate)
    } else {
      tst <- if (nlevels(d$x) == 2) wilcox.test(z ~ x, d) else kruskal.test(z ~ x, d)
      fa <- if (length(unique(d$cohort)) > 1) anova(lm(z ~ x + cohort, d)) else anova(lm(z ~ x, d))
      eff <- if (nlevels(d$x) == 2) diff(tapply(d$z, d$x, mean)) else NA
      list(type = "cat", txt = sprintf("%s p = %s  |  cohort-adjusted (lm F-test) p = %s  |  n = %d, %d cohorts",
                                       if (nlevels(d$x) == 2) "Wilcoxon" else "Kruskal-Wallis", fmt_p(tst$p.value), fmt_p(fa["x", "Pr(>F)"]),
                                       nrow(d), length(unique(d$cohort))),
           p = tst$p.value, padj = fa["x", "Pr(>F)"], eff = eff)
    }
  }

  draw_clin <- function() {
    d <- clin_df(); base_par(mar = c(5.5, 4.5, 2.5, 1))
    if (is.numeric(d$x)) {
      plot(d$x, d$z, pch = 19, cex = 0.5, col = adjustcolor(COL$main, 0.35), xlab = names(CLIN_VARS)[CLIN_VARS == input$cvar], ylab = paste(gene(), "z-score"), bty = "l", fg = COL$grey)
      lo <- lowess(d$x, d$z); lines(lo, col = COL$risk, lwd = 2.5)
    } else {
      lv <- levels(d$x); n <- length(lv); cols <- GROUP_COLS[seq_len(n)]
      plot(NA, xlim = c(0.5, n + 0.5), ylim = quantile(d$z, c(0.002, 0.998)), axes = FALSE, xlab = "", ylab = paste(gene(), "z-score"))
      axis(2, col = COL$grey); axis(1, at = seq_len(n), labels = paste0(lv, "\n(n=", as.vector(table(d$x)), ")"), tick = FALSE, padj = 0.6, cex.axis = 0.8)
      for (i in seq_len(n)) {
        zi <- d$z[d$x == lv[i]]; if (length(zi) < 3) next
        if ("violin" %in% input$cplot) {
          de <- density(zi, adjust = 1.1); w <- de$y / max(de$y) * 0.38
          polygon(c(i - w, rev(i + w)), c(de$x, rev(de$x)), col = adjustcolor(cols[i], 0.25), border = cols[i])
        }
        if ("pts" %in% input$cplot) points(jitter(rep(i, length(zi)), amount = 0.12), zi, pch = 19, cex = 0.35, col = adjustcolor(cols[i], 0.35))
        q <- quantile(zi, c(.25, .5, .75)); rect(i - 0.07, q[1], i + 0.07, q[3], col = "white", border = COL$ink)
        segments(i - 0.07, q[2], i + 0.07, q[2], lwd = 2.5)
      }
    }
    mtext(clin_test(d)$txt, side = 3, line = 0.6, adj = 0, cex = 0.78, col = "black")
  }
  output$clin_plot <- renderPlot(draw_clin(), res = 96)
  output$clin_stats <- renderUI(div(class = "small-note", clin_test(clin_df())$txt))

  draw_clin_smd <- function() {
    d <- clin_df()
    validate(need(!is.numeric(d$x), "Per-cohort SMD meta-analysis applies to categorical variables."))
    lv <- levels(d$x); a <- lv[1]; b <- if (input$cvar == "coo" && "ABC" %in% lv) "ABC" else lv[length(lv)]
    rows <- lapply(split(d, d$cohort), function(s) {
      x1 <- s$z[s$x == b]; x0 <- s$z[s$x == a]; n1 <- length(x1); n0 <- length(x0)
      if (n1 < 3 || n0 < 3) return(NULL)
      sp <- sqrt(((n1 - 1) * var(x1) + (n0 - 1) * var(x0)) / (n1 + n0 - 2)); J <- 1 - 3 / (4 * (n1 + n0) - 9)
      g <- J * (mean(x1) - mean(x0)) / sp; v <- (n1 + n0) / (n1 * n0) + g^2 / (2 * (n1 + n0))
      data.frame(label = s$cohort[1], n = n1 + n0, g = g, v = v)
    })
    r <- do.call(rbind, rows); validate(need(!is.null(r) && nrow(r) >= 2, "Fewer than 2 cohorts with both groups (>=3 patients each)."))
    m <- reml(r$g, r$v)
    # plot on SMD scale: reuse forest_plot by passing exp-free values -> custom simple forest
    base_par(mar = c(4.5, 8.5, 3, 1)); k <- nrow(r); y <- rev(seq_len(k)) + 1
    xr <- range(c(r$g - 1.96 * sqrt(r$v), r$g + 1.96 * sqrt(r$v), 0)); xr <- pmax(pmin(xr, 3), -3)
    plot(NA, xlim = xr, ylim = c(0.3, k + 1.7), axes = FALSE, xlab = sprintf("Hedges' g (%s vs %s)", b, a), ylab = "")
    axis(1, col = COL$grey); abline(v = 0, lty = 2, col = COL$grey)
    segments(r$g - 1.96 * sqrt(r$v), y, r$g + 1.96 * sqrt(r$v), y)
    points(r$g, y, pch = 15, cex = 0.6 + 1.4 * sqrt(m$w / max(m$w)), col = COL$main)
    axis(2, at = y, labels = paste0(r$label, " (", r$n, ")"), tick = FALSE, cex.axis = 0.72)
    polygon(c(m$lo, m$est, m$hi, m$est), 0.8 + c(0, 0.3, 0, -0.3), col = COL$risk, border = COL$risk)
    axis(2, at = 0.8, labels = "Pooled (REML-HKSJ)", tick = FALSE, font = 2, cex.axis = 0.78)
    mtext(sprintf("g = %.2f (%.2f to %.2f), p = %s, I² = %.0f%%", m$est, m$lo, m$hi, fmt_p(m$p), ifelse(is.na(m$I2), 0, m$I2)),
          side = 3, line = 0.6, adj = 0, cex = 0.8)
  }
  output$clin_smd <- renderPlot(draw_clin_smd(), res = 96)
  reg_dl("clin_card", function() {
    par(fig = c(0, 0.58, 0, 1)); draw_clin()
    par(fig = c(0.58, 1, 0, 1), new = TRUE); safe_draw(draw_clin_smd)
  }, 13, 6, function() paste("clinical", gene(), input$cvar, sep = "_"))

  output$clin_all <- renderDT({
    d0 <- pdat()
    out <- lapply(CLIN_VARS, function(v) {
      d <- d0; d$x <- d[[v]]; d <- d[!is.na(d$x), ]
      if (nrow(d) < 20) return(NULL)
      if (is.character(d$x)) d$x <- factor(d$x); if (is.factor(d$x)) d$x <- droplevels(d$x)
      if (is.factor(d$x) && nlevels(d$x) < 2) return(NULL)
      t <- tryCatch(clin_test(d), error = function(e) NULL); if (is.null(t)) return(NULL)
      hi <- if (is.factor(d$x)) { mm <- tapply(d$z, d$x, mean); paste0(names(which.max(mm)), " (", sprintf("%+.2f", max(mm)), ")") } else sprintf("rho %+.2f", t$eff)
      data.frame(Variable = names(CLIN_VARS)[CLIN_VARS == v], n = nrow(d), Cohorts = length(unique(d$cohort)),
                 `Highest group (mean z)` = hi, p = signif(t$p, 3), `Cohort-adjusted p` = signif(t$padj, 3), check.names = FALSE)
    })
    tb <- do.call(rbind, out); tb$`FDR (adj p)` <- signif(p.adjust(tb$`Cohort-adjusted p`, "BH"), 3)
    datatable(tb, rownames = FALSE, selection = "none", options = list(dom = "t", pageLength = 20)) |>
      formatStyle("Cohort-adjusted p", fontWeight = styleInterval(0.05, c("bold", "normal")), color = styleInterval(0.05, c(COL$risk, COL$ink)))
  })


  # ============================== Survival =======================================
  surv_all <- reactive({                     # patients with the endpoint, before cohort exclusion
    d <- pdat(); sc <- surv_cols(); d$time <- d[[sc[1]]]; d$event <- d[[sc[2]]]
    d[!is.na(d$time) & !is.na(d$event) & d$time > 0, ]
  })
  observe({
    cohs <- sort(unique(surv_all()$cohort))
    updateSelectizeInput(session, "km_excl", choices = cohs, selected = intersect(isolate(input$km_excl), cohs))
  })
  surv_base <- reactive({ d <- surv_all(); d[!d$cohort %in% input$km_excl, ] })

  km_df <- reactive({
    d <- surv_base()
    if (input$km_coh != "ALL") d <- d[d$cohort == input$km_coh, ]
    validate(need(nrow(d) >= 20, "Not enough patients with this endpoint and gene expression."))
    d$grp <- cut_within(d$z, d$cohort, input$cut); d
  })
  output$km_title <- renderText(sprintf("%s: %s by %s expression (%s)", gene(), input$endpoint, input$cut,
    if (input$km_coh != "ALL") input$km_coh else if (length(input$km_excl)) sprintf("excluding %d cohort(s)", length(input$km_excl)) else "all cohorts"))
  draw_km <- function() {
    d <- km_df(); pooled <- input$km_coh == "ALL"
    cols <- if (nlevels(d$grp) == 2) c(COL$lo, COL$hi) else if (nlevels(d$grp) == 3) c(COL$lo, COL$grey, COL$hi) else c(COL$lo, "#7fa7c9", "#e3a19a", COL$hi)
    km_plot(d$time, d$event, d$grp, title = paste(gene(), "-", input$endpoint), cols = cols,
            strata = if (pooled) d$cohort else NULL,
            sub = sprintf("%d patients, %d events, %d cohort(s)", nrow(d), sum(d$event), length(unique(d$cohort))),
            ylab = if (input$endpoint == "OS") "Overall survival" else "Progression-free survival")
  }
  reg_plot("km", draw_km, 9, 7, function() paste("KM", gene(), input$endpoint, input$cut, sep = "_"))

  # Cox table for a frame with time, event, z (per SD) and grp; stratified by cohort when > 1 cohort
  cox_summary <- function(d, adj = list(), unit = "per SD") {
    st <- if (length(unique(d$cohort)) > 1) " + strata(cohort)" else ""
    rows <- list()
    f1 <- coxph(as.formula(paste0("Surv(time, event) ~ z", st)), d)
    rows[[1]] <- c(Model = sprintf("Continuous (%s), unadjusted", unit), N = f1$n, Events = f1$nevent, cox_row(f1, "z"))
    d$hi <- as.integer(d$grp == levels(d$grp)[nlevels(d$grp)]); dd <- d[d$grp %in% levels(d$grp)[c(1, nlevels(d$grp))], ]
    f2 <- coxph(as.formula(paste0("Surv(time, event) ~ hi", st)), dd)
    rows[[2]] <- c(Model = sprintf("%s vs %s, unadjusted", levels(d$grp)[nlevels(d$grp)], levels(d$grp)[1]), N = f2$n, Events = f2$nevent, cox_row(f2, "hi"))
    for (a in adj) {
      fa <- tryCatch(coxph(as.formula(paste0("Surv(time, event) ~ z + ", paste(a, collapse = " + "), st)), d), error = function(e) NULL)
      if (!is.null(fa)) rows[[length(rows) + 1]] <- c(Model = sprintf("Continuous (%s), adjusted for %s", unit, paste(a, collapse = ", ")),
                                                      N = fa$n, Events = fa$nevent, cox_row(fa, "z"))
    }
    tb <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    for (c in c("HR", "lo", "hi")) tb[[c]] <- sprintf("%.2f", as.numeric(tb[[c]]))
    tb$p <- fmt_p(as.numeric(tb$p)); tb$`95% CI` <- paste0(tb$lo, "-", tb$hi); tb$lo <- tb$hi <- NULL
    tb[, c("Model", "N", "Events", "HR", "95% CI", "p")]
  }
  # Kaplan-Meier estimates (pooled, unstratified) at fixed times by group
  km_rates <- function(d, times, lab) {
    s <- summary(survfit(Surv(time, event) ~ grp, d), times = times, extend = TRUE)
    est <- ifelse(s$n.risk > 0, sprintf("%.0f (%.0f-%.0f)", 100 * s$surv, 100 * s$lower, 100 * s$upper), "-")
    setNames(data.frame(sub("grp=", "", s$strata), sprintf("%g mo", s$time), est, as.integer(s$n.risk)), c("Group", "Time", lab, "At risk"))
  }
  output$cox_tbl <- renderTable({
    d <- km_df(); cox_summary(d, if (length(input$adj)) list(input$adj) else list())
  }, striped = TRUE, spacing = "s", width = "100%")

  draw_sub <- function() {
    d <- surv_base()
    validate(need(nrow(d) >= 50, "Not enough patients."))
    labs <- c(coo = "COO", ipi_risk = "IPI", adv_stage = "Stage", ecog2 = "ECOG", ldh_high = "LDH", age60 = "Age", sex = "Sex", rituximab = "Rituximab")
    rows <- list(); f0 <- coxph(Surv(time, event) ~ z + strata(cohort), d)
    rows[[1]] <- data.frame(label = "All patients", n = f0$n, nevt = f0$nevent, est = coef(f0), lo = confint(f0)[1], hi = confint(f0)[2], pint = NA, head = FALSE)
    for (v in SUBGROUP_VARS) {
      dv <- d[!is.na(d[[v]]), ]; dv$x <- droplevels(factor(dv[[v]]))
      if (nlevels(dv$x) < 2) next
      pint <- tryCatch({ fi <- coxph(Surv(time, event) ~ z * x + strata(cohort), dv); fm <- coxph(Surv(time, event) ~ z + x + strata(cohort), dv)
                         anova(fm, fi)[2, "Pr(>|Chi|)"] }, error = function(e) NA)
      rows[[length(rows) + 1]] <- data.frame(label = labs[v], n = NA, nevt = NA, est = NA, lo = NA, hi = NA, pint = pint, head = TRUE)
      for (l in levels(dv$x)) {
        s <- dv[dv$x == l, ]; if (sum(s$event) < 10) next
        f <- tryCatch(coxph(Surv(time, event) ~ z + strata(cohort), s), error = function(e) NULL); if (is.null(f)) next
        rows[[length(rows) + 1]] <- data.frame(label = paste0("   ", l), n = f$n, nevt = f$nevent, est = coef(f), lo = confint(f)[1], hi = confint(f)[2], pint = NA, head = FALSE)
      }
    }
    r <- do.call(rbind, rows)
    k <- nrow(r); y <- rev(seq_len(k))
    rng <- range(c(r$lo, r$hi, 0), na.rm = TRUE); rng <- rng + c(-0.1, 0.1) * diff(rng)
    w_in <- par("fin")[1]; base_par(mai = c(0.8, 0.29 * w_in, 0.55, 0.22 * w_in))
    plot(NA, xlim = rng, ylim = c(0.5, k + 1), axes = FALSE, xlab = sprintf("HR per SD of %s (%s, cohort-stratified Cox)", gene(), input$endpoint), ylab = "")
    tk <- pretty(exp(rng), 5); tk <- tk[tk > 0]; axis(1, at = log(tk), labels = tk, col = COL$grey); abline(v = 0, lty = 2, col = COL$grey)
    xs <- grconvertX(c(0.01, 0.2, 0.80, 0.985), "nfc", "user")
    text(xs[1], k + 0.8, "Subgroup", adj = 0, font = 2, xpd = NA, cex = 0.85); text(xs[2], k + 0.8, "N / events", adj = 1, font = 2, xpd = NA, cex = 0.85)
    text(xs[3], k + 0.8, "HR (95% CI)", adj = 0, font = 2, xpd = NA, cex = 0.85); text(xs[4], k + 0.8, "P interaction", adj = 1, font = 2, xpd = NA, cex = 0.85)
    for (i in seq_len(k)) {
      rr <- r[i, ]
      text(xs[1], y[i], rr$label, adj = 0, xpd = NA, cex = 0.8, font = ifelse(rr$head | i == 1, 2, 1))
      if (!is.na(rr$est)) {
        segments(rr$lo, y[i], rr$hi, y[i]); points(rr$est, y[i], pch = if (i == 1) 18 else 15, cex = if (i == 1) 2 else 1.2, col = ifelse(rr$est > 0, COL$risk, COL$prot))
        text(xs[2], y[i], paste0(rr$n, " / ", rr$nevt), adj = 1, xpd = NA, cex = 0.78)
        text(xs[3], y[i], fmt_hr(rr$est, rr$lo, rr$hi), adj = 0, xpd = NA, cex = 0.78)
      }
      if (!is.na(rr$pint)) text(xs[4], y[i], fmt_p(rr$pint), adj = 1, xpd = NA, cex = 0.78, font = ifelse(rr$pint < 0.05, 2, 1), col = ifelse(rr$pint < 0.05, COL$risk, COL$ink))
    }
  }
  reg_plot("subgroup", draw_sub, 11, 8, function() paste("subgroups", gene(), input$endpoint, sep = "_"))

  # ============================== Treatment ======================================
  ef_pooled <- reactive(F("pooled_EF_curated.rds"))     # genome-wide early-failure ranking (used by Pathways)

  # Kaplan-Meier within R-CHOP-based therapy
  tx_ep_lab <- reactive(if (input$tx_ep == "OS") "Overall survival" else "Progression-free survival")
  tx_lab <- reactive(paste0(c(all = "R-CHOP-based therapy", "R-chemo" = "R-CHOP", "R-chemo + bortezomib" = "RB-CHOP")[[input$tx_reg]],
                            if (input$tx_coo != "all") paste0(", ", input$tx_coo), if (input$tx_ipi != "all") paste0(", IPI ", input$tx_ipi)))
  tx_df <- reactive({
    d <- pdat(); reg <- if (input$tx_reg == "all") c("R-chemo", "R-chemo + bortezomib") else input$tx_reg
    d <- d[d$tx_group %in% reg, ]
    if (input$tx_coo != "all") d <- d[d$coo %in% input$tx_coo, ]
    if (input$tx_ipi != "all") d <- d[d$ipi_risk %in% input$tx_ipi, ]
    sc <- if (input$tx_ep == "OS") c("os_m", "os_e") else c("pfs_m", "pfs_e")
    d$time <- d[[sc[1]]]; d$event <- d[[sc[2]]]
    d <- d[!is.na(d$time) & !is.na(d$event) & d$time > 0, ]
    validate(need(nrow(d) >= 30, "Fewer than 30 patients with this regimen / subgroup, the endpoint and expression of the gene."))
    d$grp <- cut_within(d$z, d$cohort, input$tx_cut); d <- d[!is.na(d$grp), ]
    validate(need(nrow(d) >= 30, "Not enough patients per cohort to form expression groups."))
    d
  })
  output$tx_km_title <- renderText(sprintf("%s: %s by %s expression, %s", gene(), tx_ep_lab(), input$tx_cut, tx_lab()))
  draw_tx_km <- function() {
    d <- tx_df()
    km_plot(d$time, d$event, d$grp, title = sprintf("%s: %s, %s", gene(), input$tx_ep, tx_lab()), cols = grp_cols(nlevels(d$grp)),
            strata = d$cohort, tmax = as.numeric(input$tx_tmax), mark_t = 24, ylab = tx_ep_lab(),
            sub = sprintf("%d patients, %d events, %d cohort(s)", nrow(d), sum(d$event), length(unique(d$cohort))))
  }
  reg_plot("tx_km", draw_tx_km, 9, 7, function() paste("KM", gene(), input$tx_ep, input$tx_reg, input$tx_coo, input$tx_ipi, input$tx_cut, sep = "_"))
  output$tx_cox <- renderTable({
    adj <- list(if (input$tx_ipi == "all") "ipi_risk", if (input$tx_ipi == "all" && input$tx_coo == "all") c("ipi_risk", "coo"))
    cox_summary(tx_df(), Filter(length, adj))
  }, striped = TRUE, spacing = "s", width = "100%")
  output$tx_rates <- renderTable(km_rates(tx_df(), c(24, 60), if (input$tx_ep == "OS") "Survival, % (95% CI)" else "Progression-free, % (95% CI)"),
                                 striped = TRUE, spacing = "s", width = "100%")
  # selected gene: early failure on R-CHOP
  output$ef_title <- renderText(sprintf("%s and early treatment failure on R-CHOP-based therapy", gene()))
  ef_df <- reactive({
    d <- pdat(); d <- d[!is.na(d$ef_rchop), ]
    validate(need(nrow(d) >= 30, "Not enough R-CHOP-based patients with 24-month follow-up and expression of this gene."))
    d$y <- d$ef_rchop == "Yes"; d
  })
  draw_ef_box <- function() {
    d <- ef_df(); d$x <- factor(ifelse(d$y, "Early failure\n(< 24 mo)", "No event\nat 24 mo"), c("No event\nat 24 mo", "Early failure\n(< 24 mo)"))
    base_par(mar = c(5, 4.5, 3, 1))
    boxplot(z ~ x, d, col = adjustcolor(c(COL$lo, COL$hi), 0.25), border = c(COL$lo, COL$hi), outline = FALSE, xlab = "", ylab = paste(gene(), "z-score"),
            frame = FALSE, xaxt = "n", whisklty = 1)
    mtext(sprintf("%s\n(n=%d)", c("No event", "Early failure"), as.vector(table(d$x))), side = 1, at = 1:2, line = 1.8, cex = 0.78)
    points(jitter(as.numeric(d$x), amount = 0.15), d$z, pch = 19, cex = 0.3, col = adjustcolor(c(COL$lo, COL$hi)[as.numeric(d$x)], 0.3))
    mtext(sprintf("Wilcoxon rank-sum test, p = %s", fmt_p(wilcox.test(z ~ y, d)$p.value)), side = 3, line = 0.6, adj = 0, cex = 0.85)
  }
  # per-cohort logistic forest of early failure on z (frame with cohort, y, z)
  ef_forest_of <- function(d, title, xlab) {
    r <- do.call(rbind, lapply(split(d, d$cohort), function(s) {
      if (sum(s$y) < 5 || sum(!s$y) < 5) return(NULL)
      f <- glm(y ~ z, binomial, s); co <- summary(f)$coefficients
      data.frame(label = s$cohort[1], n = nrow(s), nevt = sum(s$y), est = co["z", 1], se = co["z", 2])
    }))
    validate(need(!is.null(r) && nrow(r) >= 1, "No cohort with >= 5 patients in each outcome group."))
    m <- reml(r$est, r$se^2)
    df <- data.frame(label = r$label, n = r$n, nevt = r$nevt, est = r$est, lo = r$est - 1.96 * r$se, hi = r$est + 1.96 * r$se, w = m$w)
    forest_plot(df, list(est = m$est, lo = m$lo, hi = m$hi, label = "Pooled (REML-HKSJ)", n = sum(df$n), nevt = sum(df$nevt), pi_lo = NA, pi_hi = NA),
                title = title, xlab = xlab, effect = "OR", favors = NULL,
                het = if (m$k > 1) sprintf("I² = %.0f%%, overall p = %s", m$I2, fmt_p(m$p)) else NULL, cols = c("Cohort", "N", "Failures"))
  }
  # cohort-adjusted (and IPI-adjusted) logistic OR text
  ef_or_text <- function(d, unit) {
    f1 <- glm(y ~ z + cohort, binomial, d); co <- summary(f1)$coefficients["z", ]
    txt <- sprintf("Cohort-adjusted logistic model: OR %s = %.2f (%.2f-%.2f), p = %s; n = %d (%d early failures).", unit,
                   exp(co[1]), exp(co[1] - 1.96 * co[2]), exp(co[1] + 1.96 * co[2]), fmt_p(co[4]), nrow(d), sum(d$y))
    di <- d[!is.na(d$ipi_risk), ]
    if (nrow(di) > 50 && length(unique(di$ipi_risk)) > 1) {
      f2 <- glm(y ~ z + ipi_risk + cohort, binomial, di); c2 <- summary(f2)$coefficients["z", ]
      txt <- paste(txt, sprintf("Adjusted also for IPI risk: OR = %.2f (%.2f-%.2f), p = %s (n = %d).",
                                exp(c2[1]), exp(c2[1] - 1.96 * c2[2]), exp(c2[1] + 1.96 * c2[2]), fmt_p(c2[4]), nrow(di)))
    }
    div(class = "alert alert-light mt-2 mb-1", txt)
  }
  draw_ef_forest <- function() ef_forest_of(ef_df(), "Odds of early failure per SD (logistic, per cohort)", "Odds ratio per 1 SD")
  output$ef_box <- renderPlot(draw_ef_box(), res = 96)
  output$ef_forest <- renderPlot(draw_ef_forest(), res = 96)
  reg_dl("ef_card", function() {
    par(fig = c(0, 0.33, 0, 1)); draw_ef_box()
    par(fig = c(0.33, 1, 0, 1), new = TRUE); safe_draw(draw_ef_forest)
  }, 14, 6, function() paste("early_failure", gene(), sep = "_"))
  output$ef_stats <- renderUI(ef_or_text(ef_df(), "per SD"))

  # ============================== Pathways =======================================
  pretty_path <- function(x) {
    x <- gsub("_", " ", sub("^(HALLMARK|REACTOME)_", "", x))
    ifelse(nchar(x) > 62, paste0(substr(x, 1, 60), "..."), x)
  }
  gs_sets <- reactive(F("gsea_sets.rds"))
  gs_key <- reactive(if (input$gs_src == "ef") "EF_curated" else paste(input$endpoint, input$model, "curated", sep = "_"))
  gs_label <- reactive(if (input$gs_src == "ef") "early failure on R-CHOP-based therapy" else sprintf("%s, %s", input$endpoint, model_txt()))
  gs_res <- reactive({
    r <- F("gsea_results.rds")[[gs_key()]]; r <- r[r$collection == input$gs_coll & !is.na(r$p) & !is.na(r$fdr), ]
    r[order(r$p, -abs(r$NES)), ]
  })
  observe({
    req(input$nav == "Pathways"); r <- gs_res()
    sel <- isolate(if (isTRUE(input$gs_path %in% r$pathway)) input$gs_path else r$pathway[1])
    updateSelectizeInput(session, "gs_path", choices = setNames(r$pathway, pretty_path(r$pathway)), selected = sel, server = TRUE)
  })
  output$gs_n <- renderUI({
    r <- gs_res(); s <- r[r$fdr < as.numeric(input$gs_fdr), ]
    tagList(tags$h4(nrow(s), class = "mb-0"),
            div(class = "small-note", sprintf("pathways at FDR < %s (%d adverse, %d favorable) of %d tested", input$gs_fdr, sum(s$NES > 0), sum(s$NES < 0), nrow(r))))
  })
  draw_gs_bar <- function() {
    r <- gs_res(); thr <- as.numeric(input$gs_fdr)
    up <- head(r[r$NES > 0, ], 15); dn <- head(r[r$NES < 0, ], 15)
    b <- rbind(dn[order(dn$NES), ], up[order(up$NES), ])
    validate(need(nrow(b) > 0, "No pathways for this analysis."))
    sig <- b$fdr < thr; lab <- pretty_path(b$pathway)
    col <- ifelse(b$NES > 0, COL$risk, COL$prot); col <- ifelse(sig, col, adjustcolor(col, 0.3))
    # symmetric axis in whole NES units spanning the full plotting range, with light gridlines
    lim <- ceiling(max(abs(b$NES)) + 0.6); ticks <- seq(-lim, lim, 1)
    base_par(); par(mai = c(0.9, text_in(lab, 0.72) + 0.25, 0.45, 0.3))
    y <- barplot(b$NES, horiz = TRUE, plot = FALSE, space = 0.3)
    plot(NA, xlim = c(-lim, lim), ylim = c(0, max(y) + 0.8), xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
    abline(v = ticks, col = "#eef1f4", lwd = 1)
    rect(pmin(0, b$NES), y - 0.5, pmax(0, b$NES), y + 0.5, col = col, border = NA)
    abline(v = 0, col = COL$grey)
    axis(1, at = ticks, col = COL$grey); mtext("Normalized enrichment score (NES)", side = 1, line = 2.2, cex = 0.9)
    axis(2, at = y, labels = lab, tick = FALSE, cex.axis = 0.72, line = -0.6)
    text(b$NES, y, ifelse(sig, paste("FDR", fmt_p(b$fdr)), ""), pos = ifelse(b$NES > 0, 4, 2), cex = 0.6, col = "black", xpd = NA)
    mtext(sprintf("%s: %s  (faded bars: FDR >= %s)", input$gs_coll, gs_label(), input$gs_fdr), side = 3, line = 0.8, adj = 0, cex = 0.85, col = "black")
    mtext(c("← favorable", "adverse →"), side = 1, line = 3.4, at = c(-lim, lim), adj = c(0, 1), cex = 0.75, col = c(COL$prot, COL$risk))
  }
  reg_plot("gs_bar", draw_gs_bar, 11, 9, function() paste("GSEA", gs_key(), input$gs_coll, sep = "_"))
  output$dl_gsea <- downloadHandler(
    filename = function() sprintf("DLBCL-PRISM_GSEA_%s_%s.csv", gs_key(), input$gs_coll),
    content = function(file) write.csv(gs_res(), file, row.names = FALSE))

  gs_stats <- reactive({                      # ranking used by fgsea, rebuilt from the pooled table
    ef <- input$gs_src == "ef"; p <- if (ef) ef_pooled() else pooled_raw()
    p <- p[p$k >= gs_sets()$kmin[if (ef) "ef" else "surv"] & is.finite(p$se) & p$se > 0, ]
    st <- setNames(p$logHR / p$se, p$key); sort(st[!duplicated(names(st))], decreasing = TRUE)
  })
  reg_plot("gs_enrich", function() {
    req(input$gs_path); s <- gs_res(); s <- s[s$pathway == input$gs_path, ]; req(nrow(s) == 1)
    gs <- gs_sets(); members <- gs$genes[gs$sets[[input$gs_coll]][[input$gs_path]]]
    st <- gs_stats(); hit <- names(st) %in% members; N <- length(st); nh <- sum(hit)
    validate(need(nh > 0, "No genes of this pathway in the ranking."))
    a <- abs(st) * hit; es <- cumsum(a) / sum(a) - cumsum(!hit) / (N - nh)
    layout(matrix(1:2), heights = c(4, 0.9)); par(cex = 1)
    base_par(mar = c(0.3, 4.5, 3.2, 1))
    plot(seq_len(N), es, type = "l", lwd = 2, col = COL$main, xaxt = "n", xlab = "", ylab = "Running enrichment score", xaxs = "i")
    abline(h = 0, col = COL$grey, lty = 3); abline(v = which.max(abs(es)), col = COL$grey, lty = 2)
    mtext(pretty_path(input$gs_path), side = 3, line = 1.4, adj = 0, font = 2, cex = 0.95)
    mtext(sprintf("NES = %.2f, p = %s, FDR = %s, %d genes  |  %s", s$NES, fmt_p(s$p), fmt_p(s$fdr), s$size, gs_label()),
          side = 3, line = 0.3, adj = 0, cex = 0.78, col = "black")
    par(mar = c(2.2, 4.5, 0, 1)); plot(NA, xlim = c(1, N), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "", xaxs = "i")
    segments(which(hit), 0.05, which(hit), 0.95, col = adjustcolor(COL$ink, 0.6), lwd = 0.6)
    mtext(c("← genes associated with worse outcome", "genes associated with better outcome →"), side = 1, line = 0.7,
          at = c(1, N), adj = c(0, 1), cex = 0.8, col = c(COL$risk, COL$prot))
  }, 9, 5.5, function() paste("GSEA_enrichment", sub("^(HALLMARK|REACTOME)_", "", input$gs_path), gs_key(), sep = "_"))
  output$gs_le <- renderUI({
    req(input$gs_path); s <- gs_res(); s <- s[s$pathway == input$gs_path, ]; req(nrow(s) == 1)
    g <- head(strsplit(s$leading, ", ", fixed = TRUE)[[1]], 60)
    links <- lapply(g, function(x) tags$a(href = "#", x, class = "me-1",
      onclick = sprintf("Shiny.setInputValue('gs_le_gene', '%s', {priority: 'event'}); return false;", x)))
    div(class = "small-note mt-2", tags$b(sprintf("Leading-edge genes (%d%s): ", length(g), if (length(g) == 60) "+" else "")), links,
        div("Click a gene to select it in all tabs."))
  })
  observeEvent(input$gs_le_gene, set_gene(input$gs_le_gene))
  output$gs_tbl <- renderDT({
    r <- gs_res()
    d <- data.frame(Pathway = pretty_path(r$pathway), Genes = r$size, NES = round(r$NES, 2), p = signif(r$p, 3), FDR = signif(r$fdr, 3),
                    `Leading edge` = ifelse(nchar(r$leading) > 90, paste0(substr(r$leading, 1, 88), "..."), r$leading), check.names = FALSE)
    datatable(d, rownames = FALSE, selection = "single", options = list(pageLength = 15, scrollX = TRUE)) |>
      formatStyle("NES", color = styleInterval(0, c(COL$prot, COL$risk)))
  })
  observeEvent(input$gs_tbl_rows_selected, {
    r <- gs_res(); sel <- r$pathway[input$gs_tbl_rows_selected]
    updateSelectizeInput(session, "gs_path", choices = setNames(r$pathway, pretty_path(r$pathway)), selected = sel, server = TRUE)
    nav_select("gs_tabs", "Enrichment plot")
  })

  # ============================== Microenvironment ================================
  tme <- reactive(F("tme_scores_v2.rds"))
  model_lab <- reactive(if (input$model == "uni") "univariable" else "adjusted for age + IPI")
  output$tme_prog_title <- renderText(sprintf("Prognostic value of microenvironment scores (%s, %s)", input$endpoint, model_lab()))
  tme_cx <- reactive({ cx <- F("tme_cox_v2.rds"); cx[cx$set == "curated" & cx$endpoint == input$endpoint & cx$model == input$model, ] })
  tme_pool <- reactive({
    cx <- tme_cx()
    do.call(rbind, lapply(TME_CELLS, function(cel) {
      s <- cx[cx$cell == cel, ]; if (nrow(s) < 2) return(NULL)
      m <- reml(s$logHR, s$se^2)
      data.frame(cell = cel, n = sum(s$n), nevt = sum(s$events), k = m$k, est = m$est, lo = m$lo, hi = m$hi, p = m$p, I2 = m$I2)
    }))
  })
  reg_plot("tme_prog", function() {
    t <- tme_pool(); validate(need(!is.null(t), "No microenvironment data for this selection."))
    t <- t[order(t$est), ]
    df <- data.frame(label = sprintf("%s (%d)%s", t$cell, t$k, ifelse(t$p < 0.05, " *", "")), n = t$n, nevt = t$nevt, est = t$est, lo = t$lo, hi = t$hi)
    forest_plot(df, NULL, title = "", xlab = "Pooled hazard ratio per 1 SD of score", cols = c("Score (cohorts)", "N", "Events"),
                het = "REML-HKSJ pooled per-cohort Cox models.  * p < 0.05", favors = NULL)
  }, 10, 7, function() paste("TME_prognosis", input$endpoint, input$model, sep = "_"))

  tme_idx <- reactive(which(keep() & !is.na(gdat()$z)))         # patients with the current gene measured
  output$tme_cor_title <- renderText(sprintf("%s vs microenvironment: pooled within-cohort correlation", gene()))
  tme_cor_df <- reactive({
    i <- tme_idx(); z <- gdat()$z[i]; coh <- D()$clin$cohort[i]; S <- tme()$scores[i, , drop = FALSE]
    do.call(rbind, lapply(TME_CELLS, function(cel) {
      r <- do.call(rbind, lapply(split(seq_along(i), coh), function(j) {
        ok <- j[!is.na(S[j, cel])]; if (length(ok) < 20) return(NULL)
        rr <- suppressWarnings(cor(z[ok], S[ok, cel])); if (!is.finite(rr)) NULL else data.frame(r = rr, n = length(ok))
      }))
      if (is.null(r) || nrow(r) < 2) return(NULL)
      m <- reml(atanh(pmin(pmax(r$r, -0.995), 0.995)), 1 / (r$n - 3))
      data.frame(cell = cel, r = tanh(m$est), lo = tanh(m$lo), hi = tanh(m$hi), p = m$p, k = m$k, I2 = m$I2,
                 marker = gene() %in% tme()$genes[[cel]])
    }))
  })
  reg_plot("tme_cor", function() {
    t <- tme_cor_df(); validate(need(!is.null(t) && nrow(t) > 0, "Not enough patients with this gene and the microenvironment scores."))
    t <- t[order(t$r), ]; y <- seq_len(nrow(t)); lab <- paste0(t$cell, ifelse(t$marker, " *", ""))
    base_par(); par(mai = c(0.9, text_in(lab, 0.8) + 0.25, 0.45, 1.6))
    sig <- t$p < 0.05; col <- ifelse(sig, ifelse(t$r > 0, COL$risk, COL$prot), COL$grey)
    plot(t$r, y, xlim = range(c(-0.1, 0.1, t$lo, t$hi)), ylim = c(0.5, nrow(t) + 0.5), pch = 19, cex = 1.2, col = col,
         yaxt = "n", xlab = "Pooled correlation (r, 95% CI)", ylab = "")
    abline(v = 0, col = COL$grey, lty = 3); segments(t$lo, y, t$hi, y, col = col, lwd = 1.6)
    axis(2, at = y, labels = lab, tick = FALSE, cex.axis = 0.8)
    text(par("usr")[2], y, sprintf("r = %.2f  p = %s", t$r, fmt_p(t$p)), adj = 0, pos = 4, cex = 0.68, col = "black", xpd = NA)
    mtext(sprintf("%s; k = %d-%d cohorts", gene(), min(t$k), max(t$k)), side = 3, line = 0.6, adj = 0, cex = 0.85, col = "black")
  }, 8, 7, function() paste("TME_correlation", gene(), sep = "_"))

  output$tme_markers <- renderUI({
    tm <- tme(); cv <- tm$coverage; cv <- cv[cv$cell == input$tme_cell & cv$cohort %in% unique(D()$clin$cohort[keep()]), ]
    div(class = "small-note", tags$b("Markers: "), paste(tm$genes[[input$tme_cell]], collapse = ", "),
        div(class = "mt-1", sprintf("Score available in %d of %d cohorts (>= 50%% of markers measured).", sum(cv$n_genes >= ceiling(cv$of / 2)), nrow(cv))))
  })
  output$tme_sel_title <- renderText(sprintf("%s score", input$tme_cell))
  tme_sel <- reactive({
    i <- which(keep()); cl <- D()$clin[i, ]
    d <- data.frame(cohort = cl$cohort, coo = cl$coo, ef = cl$ef_rchop, ipi_risk = cl$ipi_risk, z = gdat()$z[i], s = tme()$scores[i, input$tme_cell])
    d[!is.na(d$s), ]
  })
  draw_tme_scatter <- function() {
    d <- tme_sel(); d <- d[!is.na(d$z), ]; validate(need(nrow(d) >= 30, "Not enough patients with both measurements."))
    base_par(mar = c(4.5, 4.5, 3, 1))
    plot(d$s, d$z, pch = 19, cex = 0.35, col = adjustcolor(COL$main, 0.3), xlab = paste(input$tme_cell, "score"), ylab = paste(gene(), "z-score"), bty = "l", fg = COL$grey)
    abline(lm(z ~ s, d), col = COL$risk, lwd = 2)
    t <- tme_cor_df(); t <- t[t$cell == input$tme_cell, ]
    if (nrow(t)) mtext(c(sprintf("%s: pooled r = %.2f (%.2f to %.2f)", gene(), t$r, t$lo, t$hi), sprintf("p = %s (REML-HKSJ)", fmt_p(t$p))), side = 3, line = c(1.3, 0.3), adj = 0, cex = 0.78)
  }
  draw_tme_coo <- function() {
    d <- tme_sel(); d <- d[!is.na(d$coo), ]; validate(need(nrow(d) >= 30 && length(unique(d$coo)) > 1, "Cell of origin not available."))
    d$coo <- droplevels(d$coo); p <- anova(lm(s ~ coo + cohort, d))["coo", "Pr(>F)"]
    base_par(mar = c(4.5, 4.5, 3, 1)); cols <- GROUP_COLS[seq_along(levels(d$coo))]
    boxplot(s ~ coo, d, col = adjustcolor(cols, 0.25), border = cols, outline = FALSE, xlab = "", ylab = paste(input$tme_cell, "score"), frame = FALSE, xaxt = "n", whisklty = 1)
    mtext(sprintf("%s\n(n=%d)", sub("Unclassified", "Unclass.", levels(d$coo)), as.vector(table(d$coo))), side = 1, at = seq_along(levels(d$coo)), line = 1.8, cex = 0.72)
    mtext(c("By cell of origin", sprintf("cohort-adjusted p = %s", fmt_p(p))), side = 3, line = c(1.3, 0.3), adj = 0, cex = 0.78)
  }
  tme_ef_df <- reactive({
    d <- tme_sel(); d <- d[!is.na(d$ef), ]
    validate(need(nrow(d) >= 30 && length(unique(d$ef)) > 1, "Not enough R-CHOP-based patients with 24-month follow-up and this score."))
    d$y <- d$ef == "Yes"; d$z <- d$s; d
  })
  draw_tme_ef <- function() {
    d <- tme_ef_df()
    d$x <- factor(ifelse(d$y, "Early failure", "No event"), c("No event", "Early failure"))
    base_par(mar = c(4.5, 4.5, 3, 1))
    boxplot(s ~ x, d, col = adjustcolor(c(COL$lo, COL$hi), 0.25), border = c(COL$lo, COL$hi), outline = FALSE, xlab = "", ylab = paste(input$tme_cell, "score"),
            frame = FALSE, xaxt = "n", whisklty = 1)
    points(jitter(as.numeric(d$x), amount = 0.15), d$s, pch = 19, cex = 0.3, col = adjustcolor(c(COL$lo, COL$hi)[as.numeric(d$x)], 0.3))
    mtext(sprintf("%s
(n=%d)", c("No event", "Early failure"), as.vector(table(d$x))), side = 1, at = 1:2, line = 1.8, cex = 0.78)
    mtext(sprintf("Wilcoxon rank-sum test, p = %s", fmt_p(wilcox.test(s ~ y, d)$p.value)), side = 3, line = 0.6, adj = 0, cex = 0.85)
  }
  draw_tme_ef_forest <- function() ef_forest_of(tme_ef_df(), sprintf("%s score: odds of early failure on R-CHOP-based therapy", input$tme_cell),
                                                "Odds ratio per 1 SD of score")
  output$tme_scatter <- renderPlot(draw_tme_scatter(), res = 96)
  output$tme_coo <- renderPlot(draw_tme_coo(), res = 96)
  output$tme_ef <- renderPlot(draw_tme_ef(), res = 96)
  output$tme_ef_forest <- renderPlot(draw_tme_ef_forest(), res = 96)
  output$tme_ef_stats <- renderUI(ef_or_text(tme_ef_df(), "per SD of score"))
  reg_dl("tme_txr", function() {
    par(fig = c(0, 0.33, 0, 1)); draw_tme_ef()
    par(fig = c(0.33, 1, 0, 1), new = TRUE); safe_draw(draw_tme_ef_forest)
  }, 14, 6, function() paste("TME_early_failure", input$tme_cell, sep = "_"))
  reg_dl("tme_assoc", function() { par(mfrow = c(1, 2)); safe_draw(draw_tme_scatter); safe_draw(draw_tme_coo) },
         10, 5, function() paste("TME_associations", input$tme_cell, gene(), sep = "_"))

  # Kaplan-Meier by score groups (endpoint and adjustment follow the sidebar)
  tme_km_df <- reactive({
    i <- which(keep()); d <- D()$clin[i, ]; d$z <- tme()$scores[i, input$tme_cell]
    sc <- surv_cols(); d$time <- d[[sc[1]]]; d$event <- d[[sc[2]]]
    if (input$tme_pop == "rchop") d <- d[d$tx_group %in% c("R-chemo", "R-chemo + bortezomib"), ]
    if (input$tme_coo_f != "all") d <- d[d$coo %in% input$tme_coo_f, ]
    d <- d[!is.na(d$z) & !is.na(d$time) & !is.na(d$event) & d$time > 0, ]
    validate(need(nrow(d) >= 30, "Not enough patients with this score and endpoint."))
    d$grp <- cut_within(d$z, d$cohort, input$tme_cut); d <- d[!is.na(d$grp), ]
    d
  })
  tme_km_lab <- reactive(paste0(input$tme_cell, " score, ", input$endpoint,
                                if (input$tme_pop == "rchop") ", R-CHOP-based" else "", if (input$tme_coo_f != "all") paste0(", ", input$tme_coo_f) else ""))
  reg_plot("tme_km", function() {
    d <- tme_km_df()
    km_plot(d$time, d$event, d$grp, title = tme_km_lab(), cols = grp_cols(nlevels(d$grp)), strata = d$cohort,
            sub = sprintf("%d patients, %d events, %d cohorts", nrow(d), sum(d$event), length(unique(d$cohort))),
            ylab = if (input$endpoint == "OS") "Overall survival" else "Progression-free survival")
  }, 9, 7, function() paste("TME_KM", input$tme_cell, input$endpoint, input$tme_pop, input$tme_coo_f, input$tme_cut, sep = "_"))
  output$tme_km_cox <- renderTable({
    adj <- if (input$model == "multi") list(c("age", "ipi_risk")) else list()
    cox_summary(tme_km_df(), adj, unit = "per SD of score")
  }, striped = TRUE, spacing = "s", width = "100%")

  reg_plot("tme_forest", function() {
    s <- tme_cx(); s <- s[s$cell == input$tme_cell, ]; validate(need(nrow(s) >= 2, "Score available in fewer than 2 cohorts."))
    m <- reml(s$logHR, s$se^2)
    df <- data.frame(label = s$cohort, n = s$n, nevt = s$events, est = s$logHR, lo = s$logHR - 1.96 * s$se, hi = s$logHR + 1.96 * s$se, w = m$w)
    df <- df[order(-df$w), ]
    forest_plot(df, list(est = m$est, lo = m$lo, hi = m$hi, label = "Pooled (REML-HKSJ)", n = sum(df$n), nevt = sum(df$nevt), pi_lo = NA, pi_hi = NA),
                title = sprintf("%s score and %s by cohort (%s)", input$tme_cell, input$endpoint, model_lab()), xlab = "Hazard ratio per 1 SD of score",
                het = sprintf("I² = %.0f%%, τ² = %.4f, overall p = %s", m$I2, m$tau2, fmt_p(m$p)), favors = c("favors low score", "favors high score"))
  }, 11, function() 2.6 + 0.3 * sum(tme_cx()$cell == input$tme_cell), function() paste("TME_forest", input$tme_cell, input$endpoint, input$model, sep = "_"))
  # ============================== Signature ======================================
  observeEvent(input$sig_top, {
    h <- head(hits(), input$sig_n)
    if (nrow(h) == 0) showNotification("No Discovery hits with the current filters.", type = "warning")
    else updateTextAreaInput(session, "sig_genes", value = paste(h$key, collapse = "\n"))
  })

  sig <- eventReactive(input$sig_run, {
    g <- unique(toupper(trimws(unlist(strsplit(input$sig_genes, "[,;\n\t ]+")))))
    g <- g[nzchar(g)]; g <- head(g[g %in% D()$gidx$key[D()$gidx$expr]], 60)
    validate(need(length(g) >= 2, "Enter at least 2 genes with expression data."))
    ep <- input$endpoint; md <- input$model; excl <- excluded_cohorts()
    withProgress(message = "Building signature...", value = 0, {
      Z <- matrix(NA_real_, nrow(D()$clin), length(g), dimnames = list(NULL, g)); est <- list()
      for (i in seq_along(g)) {
        incProgress(0.7 / length(g))
        gi <- D()$gidx[D()$gidx$key == g[i], ]; ch <- F(sprintf("chunks/c%04d.rds", gi$chunk))
        Z[, i] <- decode_z(ch$expr[g[i], ])
        r <- ch$res[[g[i]]]; est[[g[i]]] <- if (is.null(r)) NULL else r[r$endpoint == ep & r$model == md & !r$cohort %in% excl, ]
      }
      cl <- D()$clin; kp <- keep(); cohorts <- unique(cl$cohort[kp])
      weight_for <- function(drop) sapply(g, function(x) {
        r <- est[[x]]; if (is.null(r)) return(NA); r <- r[r$cohort != drop, ]; if (nrow(r) < 2) return(NA)
        b <- reml(r$logHR, r$se^2)$est; if (input$sig_w == "sign") sign(b) else b })
      w_all <- weight_for("__none__")
      score <- rep(NA_real_, nrow(cl))
      for (c in cohorts) {
        i <- which(cl$cohort == c & kp); w <- if (isTRUE(input$sig_loco)) weight_for(c) else w_all
        Zi <- Z[i, , drop = FALSE]; ok <- !is.na(w)
        if (sum(ok) < 2) next
        Zi <- Zi[, ok, drop = FALSE]; wi <- w[ok]
        avail <- colSums(!is.na(Zi)) > 0; if (mean(avail) < 0.5) next
        s <- apply(Zi, 1, function(z) { a <- !is.na(z); if (!any(a)) NA else sum(z[a] * wi[a]) / sum(abs(wi[a])) })
        score[i] <- (s - mean(s, na.rm = TRUE)) / sd(s, na.rm = TRUE)
      }
      incProgress(0.2)
      d <- cl; d$score <- score; sc <- if (ep == "OS") c("os_m", "os_e") else c("pfs_m", "pfs_e")
      d$time <- d[[sc[1]]]; d$event <- d[[sc[2]]]
      d <- d[kp & !is.na(d$score) & !is.na(d$time) & !is.na(d$event) & d$time > 0, ]
      list(d = d, genes = g, w = w_all, ep = ep)
    })
  })

  draw_sig_forest <- function() {
    s <- sig(); d <- s$d
    r <- do.call(rbind, lapply(split(d, d$cohort), function(x) {
      if (sum(x$event) < 5) return(NULL)
      f <- coxph(Surv(time, event) ~ score, x)
      data.frame(label = x$cohort[1], n = f$n, nevt = f$nevent, est = coef(f), se = sqrt(vcov(f)[1]))
    }))
    validate(need(!is.null(r), "No cohort with enough events."))
    m <- reml(r$est, r$se^2)
    df <- data.frame(label = r$label, n = r$n, nevt = r$nevt, est = r$est, lo = r$est - 1.96 * r$se, hi = r$est + 1.96 * r$se, w = m$w)
    df <- df[order(-df$w), ]
    forest_plot(df, list(est = m$est, lo = m$lo, hi = m$hi, label = "Pooled (REML-HKSJ)", n = sum(df$n), nevt = sum(df$nevt), pi_lo = NA, pi_hi = NA),
                title = sprintf("%d-gene score, %s (%s weights)", length(s$genes), s$ep, if (isTRUE(input$sig_loco)) "leave-one-cohort-out" else "full-data"),
                het = if (m$k > 1) sprintf("I² = %.0f%%, τ² = %.4f, overall p = %s", m$I2, m$tau2, fmt_p(m$p)) else NULL,
                xlab = "Hazard ratio per 1 SD of score")
  }
  reg_plot("sig_forest", draw_sig_forest, 11, 8, function() "signature_forest")
  reg_plot("sig_km", function() {
    s <- sig(); d <- s$d; d$grp <- cut_within(d$score, d$cohort, "tertile")
    km_plot(d$time, d$event, d$grp, title = sprintf("%d-gene score tertiles (%s)", length(s$genes), s$ep),
            cols = c(COL$lo, COL$grey, COL$hi), strata = d$cohort, sub = sprintf("%d patients, %d cohorts", nrow(d), length(unique(d$cohort))))
  }, 9, 7, function() "signature_KM")
  cidx <- reactive({
    d <- sig()$d; d <- d[!is.na(d$ipi_risk), ]
    validate(need(nrow(d) >= 50, "Not enough patients with IPI."))
    d$ipi <- as.numeric(d$ipi_risk)
    fits <- list(IPI = coxph(Surv(time, event) ~ ipi + strata(cohort), d),
                 Score = coxph(Surv(time, event) ~ score + strata(cohort), d),
                 `IPI + score` = coxph(Surv(time, event) ~ ipi + score + strata(cohort), d))
    data.frame(Model = names(fits), C = sapply(fits, function(f) f$concordance["concordance"]),
               SE = sapply(fits, function(f) f$concordance["std"]), AIC = sapply(fits, AIC), n = nrow(d), cohorts = length(unique(d$cohort)),
               LRT_p = c(NA, NA, anova(fits[[1]], fits[[3]])[2, "Pr(>|Chi|)"]))
  })
  reg_plot("sig_cidx", function() {
    t <- cidx(); base_par(mar = c(4.5, 8, 3, 1))
    y <- barplot(t$C, names.arg = t$Model, horiz = TRUE, xlim = c(0.5, max(t$C + 2 * t$SE, 0.75)), xpd = FALSE,
                 col = c(COL$grey, COL$main, COL$risk), border = NA, xlab = "Harrell's C-index (stratified by cohort)")
    arrows(t$C - 1.96 * t$SE, y, t$C + 1.96 * t$SE, y, angle = 90, code = 3, length = 0.05)
    abline(v = 0.5, lty = 3); mtext(sprintf("n = %d patients with IPI in %d cohorts; LRT for adding score to IPI p = %s", t$n[1], t$cohorts[1], fmt_p(t$LRT_p[3])), side = 3, line = 0.6, adj = 0, cex = 0.85)
  }, 8, 5, function() "signature_Cindex")
  output$sig_cidx_tbl <- renderTable({ t <- cidx(); data.frame(Model = t$Model, `C-index` = sprintf("%.3f (%.3f-%.3f)", t$C, t$C - 1.96 * t$SE, t$C + 1.96 * t$SE), AIC = round(t$AIC, 1), check.names = FALSE) }, striped = TRUE, spacing = "s")
  output$sig_used <- renderTable({ s <- sig(); data.frame(Gene = s$genes, `Weight (full data)` = round(s$w, 3), check.names = FALSE) }, striped = TRUE, spacing = "s")
}

shinyApp(ui, server)
