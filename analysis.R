# Attribution of licensed open access to Polish national transformative agreements and to the
# open-access mandate of the National Science Centre (NCN).
#
# Usage: Rscript analysis.R AGGREGATES OUTPUT
#   AGGREGATES    folder with the aggregate count tables (data/aggregates)
#   OUTPUT        new folder for the result tables, figures and the manifest
#
# Every interval is computed with the R package pnsbib (exact linear program through the native
# backend, Li-Pearl and Shu-Wang-Li formulas where they apply). Identification intervals at point
# inputs, assumption regions and the sampling projection are written to separate tables.
# No random numbers are used.

args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
agg <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
if (file.exists(out)) stop("Output directory exists; choose a new one.")
dir.create(out, recursive = TRUE)
suppressMessages(library(pnsbib))
stopifnot(as.character(packageVersion("pnsbib")) %in% c("0.2.0", "0.2.1"))
BACKEND <- "zig"
ALPHA <- 0.05            # simultaneous level of the sampling projection
TOL <- 1e-9

read_tsv <- function(path) utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, name) utils::write.table(x, file.path(out, name), sep = "\t", quote = FALSE,
                                                  row.names = FALSE, na = "NA")
checks <- data.frame(check = character(), passed = logical(), detail = character())
record <- function(name, ok, detail = "") {
  checks[nrow(checks) + 1L, ] <<- list(name, isTRUE(ok), detail)
  invisible(ok)
}

# ---------------------------------------------------------------- pnsbib layer
# Target population: units exposed to the policy (X = policy). The untreated row
# of the factual table is zero; P(Y_none) is the counterfactual margin of the
# exposed units, supplied by the comparison design.
STATES <- c("none", "policy")
Q_PNS <- poc_query(c(none = "other", policy = "event"))
Q_PN <- poc_query(c(none = "other"), observed_x = "policy", observed_y = "event", conditional = TRUE)
Q_KEEP <- poc_query(c(none = "event", policy = "event"))

binary_observed <- function(pL) matrix(c(0, 0, 1 - pL, pL), 2, 2, byrow = TRUE,
                                       dimnames = list(STATES, c("other", "event")))
binary_model <- function(pL, qL) {
  a <- matrix(c(1 - qL, qL, 1 - pL, pL), 2, 2, byrow = TRUE, dimnames = list(STATES, c("other", "event")))
  poc_model(binary_observed(pL), a, lp_backend = BACKEND)
}
band_model <- function(pL, lo, hi) {
  poc_partial_model(binary_observed(pL),
    list(poc_constraint(poc_query(c(none = "event")), max(0, lo), min(1, hi), label = "counterfactual_event_share")),
    lp_backend = BACKEND)
}
interval <- function(r) c(lower = r$lower[1], upper = r$upper[1])
partial <- function(model, query) interval(poc_partial_bounds(model, query))

category_model <- function(p, q) {
  o <- rbind(none = 0 * p, policy = p); a <- rbind(none = q, policy = p)
  dimnames(o) <- dimnames(a) <- list(STATES, names(p))
  poc_model(o, a, lp_backend = BACKEND)
}
monotone_binary <- function(pL, qL) tryCatch(poc_partial_model(binary_observed(pL), list(
  poc_constraint(poc_query(c(none = "event")), qL, label = "counterfactual_event_share"),
  poc_constraint(poc_query(c(none = "event", policy = "other")), 0, label = "no_removal")),
  lp_backend = BACKEND), error = function(e) NULL)
monotone_category <- function(p, q, event) {
  cats <- names(p)
  o <- rbind(none = 0 * p, policy = p); dimnames(o) <- list(STATES, cats)
  cons <- lapply(cats, function(k) poc_constraint(poc_query(c(none = k)), q[[k]], label = paste0("q_", k)))
  for (j in cats) for (k in setdiff(cats, event)) if (j != k)
    cons[[length(cons) + 1L]] <- poc_constraint(poc_query(c(none = j, policy = k)), 0,
                                                label = paste0("no_move_", j, "_", k))
  tryCatch(poc_partial_model(o, cons, lp_backend = BACKEND), error = function(e) NULL)
}

# Sampling inputs: Clopper-Pearson for the exposed share, Hoeffding for the
# linear counterfactual estimator; Bonferroni over the two inputs.
clopper <- function(x, n, level) c(lo = if (x == 0) 0 else stats::qbeta((1 - level) / 2, x, n - x + 1),
                                   hi = if (x == n) 1 else stats::qbeta(1 - (1 - level) / 2, x + 1, n - x))
hoeffding <- function(v, level) sqrt(v * log(2 / (1 - level)) / 2)
# Outer interval over p in [p_lo, p_hi] and a counterfactual band [q_lo, q_hi].
# Each endpoint function is monotone in p, so both ends of the p-interval suffice.
region <- function(p_lo, p_hi, q_lo, q_hi) {
  ends <- lapply(unique(c(p_lo, p_hi)), function(pv) {
    m <- band_model(pv, q_lo, q_hi)
    c(partial(m, Q_PNS), partial(m, Q_PN))
  })
  e <- do.call(rbind, ends)
  c(PNS_lower = min(e[, 1]), PNS_upper = max(e[, 2]), PN_lower = min(e[, 3]), PN_upper = max(e[, 4]))
}

# Complete-margin calculations for one exposed population and one counterfactual margin.
attribution <- function(p, q, event, labels) {
  pL <- p[[event]]; qL <- q[[event]]
  m2 <- binary_model(pL, qL)
  res <- list()
  for (est in c("PNS", "PN")) {
    query <- if (est == "PNS") Q_PNS else Q_PN
    lp <- poc_exact(m2, query); li <- poc_bounds(m2, query)
    shu <- tryCatch(poc_shu2026(m2, query), error = function(e) NULL)
    res[[est]] <- data.frame(labels, estimand = est,
      exact_lower = lp$lower, exact_upper = lp$upper, exact_sharpness = lp$sharpness,
      lipearl_lower = li$lower, lipearl_upper = li$upper, lipearl_method = as.character(li$method),
      lipearl_sharpness = li$sharpness,
      shu_lower = if (is.null(shu)) NA else shu$lower, shu_upper = if (is.null(shu)) NA else shu$upper,
      stringsAsFactors = FALSE)
  }
  out_b <- do.call(rbind, res)
  m4 <- category_model(p, q)
  flows <- do.call(rbind, lapply(setdiff(names(p), event), function(k) {
    query <- poc_query(setNames(c(k, event), STATES))
    lp <- poc_exact(m4, query); li <- poc_bounds(m4, query)
    data.frame(labels, from = k, exact_lower = lp$lower, exact_upper = lp$upper,
               lipearl_lower = li$lower, lipearl_upper = li$upper, lipearl_method = as.character(li$method),
               lipearl_sharpness = li$sharpness, stringsAsFactors = FALSE)
  }))
  mb <- monotone_binary(pL, qL)
  mc <- monotone_category(p, q, event)
  mono <- data.frame(labels, binary_compatible = !is.null(mb),
    PNS = if (is.null(mb)) NA else partial(mb, Q_PNS)[["lower"]],
    PNS_upper = if (is.null(mb)) NA else partial(mb, Q_PNS)[["upper"]],
    PN = if (is.null(mb)) NA else partial(mb, Q_PN)[["lower"]],
    PN_upper = if (is.null(mb)) NA else partial(mb, Q_PN)[["upper"]],
    categories_compatible = !is.null(mc), stringsAsFactors = FALSE)
  for (k in setdiff(names(p), event)) {
    mono[[paste0("from_", k)]] <- if (is.null(mc)) NA else
      partial(mc, poc_query(setNames(c(k, event), STATES)))[["lower"]]
    mono[[paste0("change_", k)]] <- q[[k]] - p[[k]]
  }
  mono$retained <- if (is.null(mc)) NA else partial(mc, poc_query(setNames(c(event, event), STATES)))[["lower"]]
  list(bounds = out_b, flows = flows, monotone = mono)
}

# ========================================================= national agreements
K <- c("closed", "green", "bronze", "licensed")
ac <- read_tsv(file.path(agg, "agreement_counts.tsv"))
counts <- function(pub, corr, years) {
  r <- ac[ac$publisher == pub & ac$corresponding_author == corr & ac$year %in% years, K, drop = FALSE]
  if (!nrow(r)) NULL else colSums(r)
}
shares <- function(v) v / sum(v)
agreement_years <- function(pub)
  sort(unique(ac$year[ac$publisher == pub & ac$corresponding_author == "PL" & ac$agreement_year == 1]))
GROUPS <- list(
  "2020" = list(label = "Agreements from 2020", publishers = c("ACS", "CUP", "IOP", "OUP"),
                pre = 2016:2019, years = 2020:2023),
  "2019" = list(label = "Agreements from 2019", publishers = c("Elsevier", "Springer Nature"),
                pre = 2016:2018, years = 2019:2023))
record("IEEE lacks foreign corresponding-author articles before its agreement",
       is.null(counts("IEEE", "foreign", 2016:2019)),
       "IEEE is excluded because the trend comparison needs these articles")

publisher_counterfactual <- function(pub, t, pre) {
  pl_pre <- counts(pub, "PL", pre); f_pre <- counts(pub, "foreign", pre); f_t <- counts(pub, "foreign", t)
  if (is.null(pl_pre) || is.null(f_pre) || is.null(f_t)) return(NULL)
  qa <- shares(pl_pre)
  raw <- qa + shares(f_t) - shares(f_pre)
  qb <- pmin(pmax(raw, 0), 1); qb <- qb / sum(qb)
  list(qa = qa, qb = qb, n_pl_pre = sum(pl_pre), n_for_pre = sum(f_pre), n_for_t = sum(f_t),
       clipped = any(raw < 0 | raw > 1), d_for = shares(f_t)[["licensed"]] - shares(f_pre)[["licensed"]],
       n_for_t_weight = sum(f_t))
}
pooled <- function(group, t, pre, covered_only = TRUE) {
  g <- GROUPS[[group]]; parts <- list()
  for (pub in g$publishers) {
    if (covered_only && !(t %in% agreement_years(pub))) next
    pt <- counts(pub, "PL", t); cf <- publisher_counterfactual(pub, t, pre)
    if (is.null(pt) || is.null(cf)) next
    parts[[pub]] <- list(n = sum(pt), counts = pt, cf = cf)
  }
  n <- sum(vapply(parts, function(x) x$n, 0)); w <- vapply(parts, function(x) x$n / n, 0)
  p <- Reduce(`+`, lapply(parts, function(x) x$counts)) / n
  qa <- Reduce(`+`, Map(function(x, wi) wi * x$cf$qa, parts, w))
  qb <- Reduce(`+`, Map(function(x, wi) wi * x$cf$qb, parts, w))
  va <- sum(mapply(function(x, wi) wi^2 / x$cf$n_pl_pre, parts, w))
  vb <- sum(mapply(function(x, wi) wi^2 * (1 / x$cf$n_pl_pre + 1 / x$cf$n_for_pre + 1 / x$cf$n_for_t), parts, w))
  nf <- vapply(parts, function(x) x$cf$n_for_t_weight, 0)
  d_for <- sum(nf * vapply(parts, function(x) x$cf$d_for, 0)) / sum(nf)
  list(publishers = names(parts), n = n, x = sum(vapply(parts, function(x) x$counts[["licensed"]], 0)),
       p = p, q = list(a = qa, b = qb), v = c(a = va, b = vb), d_for = d_for,
       clipped = any(vapply(parts, function(x) x$cf$clipped, TRUE)))
}

# Negative-control year: the last year before the agreements, predicted from earlier years.
margins <- list()
for (g in names(GROUPS)) {
  ctrl <- max(GROUPS[[g]]$pre)
  r <- pooled(g, ctrl, 2016:(ctrl - 1), covered_only = FALSE)
  margins[[g]] <- c(a = abs(r$p[["licensed"]] - r$q$a[["licensed"]]),
                    b = abs(r$p[["licensed"]] - r$q$b[["licensed"]]), year = ctrl)
}

# Publisher choice: change in the Polish-minus-foreign gap in the share of
# agreement publishers, relative to 2016-2018, divided by the Polish share.
vc <- read_tsv(file.path(agg, "venue_counts.tsv"))
venue <- function(corr, y) { r <- vc[vc$corresponding_author == corr & vc$year == y, ]
  c(share = (r$agreement_2019 + r$agreement_2020) / r$n, n = r$n) }
gap <- function(y) { a <- venue("PL", y); b <- venue("foreign", y)
  c(gap = a[["share"]] - b[["share"]],
    se = sqrt(a[["share"]] * (1 - a[["share"]]) / a[["n"]] + b[["share"]] * (1 - b[["share"]]) / b[["n"]])) }
base <- sapply(2016:2018, gap)
base_gap <- mean(base["gap", ]); base_se <- sqrt(sum(base["se", ]^2)) / 3
vo <- read_tsv(file.path(agg, "venue_outside_open.tsv"))
choice <- do.call(rbind, lapply(2020:2023, function(y) {
  gy <- gap(y); change <- gy[["gap"]] - base_gap
  data.frame(year = y, polish_share = venue("PL", y)[["share"]], gap = gy[["gap"]], baseline_gap = base_gap,
             change = change, z = change / sqrt(gy[["se"]]^2 + base_se^2),
             attracted = max(0, change / venue("PL", y)[["share"]]),
             open_elsewhere = vo$open_at_publisher[vo$year == y] / vo$n_polish_outside_agreement_publishers[vo$year == y])
}))
write_tsv(choice, "agreement_publisher_choice.tsv")

inputs <- bounds <- flows <- mono <- regions <- breakdown <- list()
for (g in names(GROUPS)) for (t in GROUPS[[g]]$years) {
  r <- pooled(g, t, GROUPS[[g]]$pre)
  attracted <- if (t %in% choice$year) choice$attracted[choice$year == t] else NA
  inputs[[length(inputs) + 1L]] <- data.frame(group = g, year = t, publishers = paste(r$publishers, collapse = ", "),
    n = r$n, licensed_count = r$x,
    p_closed = r$p[["closed"]], p_green = r$p[["green"]], p_bronze = r$p[["bronze"]], p_licensed = r$p[["licensed"]],
    qa_closed = r$q$a[["closed"]], qa_green = r$q$a[["green"]], qa_bronze = r$q$a[["bronze"]], qa_licensed = r$q$a[["licensed"]],
    qb_closed = r$q$b[["closed"]], qb_green = r$q$b[["green"]], qb_bronze = r$q$b[["bronze"]], qb_licensed = r$q$b[["licensed"]],
    margin_a = margins[[g]][["a"]], margin_b = margins[[g]][["b"]], control_year = margins[[g]][["year"]],
    foreign_change = r$d_for, attracted = attracted, qb_clipped = r$clipped, stringsAsFactors = FALSE)
  p_ci <- clopper(r$x, r$n, 1 - ALPHA / 2)
  for (m in c("a", "b")) {
    lab <- data.frame(study = "agreements", group = g, year = t, method = m, stringsAsFactors = FALSE)
    q <- r$q[[m]]; qL <- q[["licensed"]]; pL <- r$p[["licensed"]]; d <- margins[[g]][[m]]
    a <- attribution(r$p, q, "licensed", lab)
    bounds[[length(bounds) + 1L]] <- a$bounds; flows[[length(flows) + 1L]] <- a$flows
    mono[[length(mono) + 1L]] <- a$monotone
    eps <- hoeffding(r$v[[m]], 1 - ALPHA / 2)
    reg <- rbind(
      point = region(pL, pL, qL, qL),
      margin = region(pL, pL, qL - d, qL + d),
      margin_choice = if (is.na(attracted)) rep(NA, 4) else
        region(pL, pL, (1 - attracted) * (qL - d), qL + d + attracted * (1 - qL)),
      sampling = region(p_ci[["lo"]], p_ci[["hi"]], qL - eps, qL + eps),
      sampling_margin_choice = if (is.na(attracted)) region(p_ci[["lo"]], p_ci[["hi"]], qL - eps - d, qL + eps + d) else
        region(p_ci[["lo"]], p_ci[["hi"]], (1 - attracted) * (qL - eps - d), qL + eps + d + attracted * (1 - qL)))
    regions[[length(regions) + 1L]] <- data.frame(lab, region = rownames(reg), reg, row.names = NULL,
                                                  stringsAsFactors = FALSE)
    breakdown[[length(breakdown) + 1L]] <- data.frame(lab, margin = d, breakdown_margin = max(0, pL - qL),
      ratio = max(0, pL - qL) / d, attracted = attracted,
      breakdown_attracted = max(0, (pL - qL - d) / (1 - qL)),
      p_lo = p_ci[["lo"]], p_hi = p_ci[["hi"]], q_halfwidth = eps, stringsAsFactors = FALSE)
  }
}
agreement_inputs <- do.call(rbind, inputs)
write_tsv(agreement_inputs, "agreement_inputs.tsv")

# Pooled exposed populations: years as disjoint strata in one partial model.
pooled_rows <- list()
pool_strata <- function(ps, ns, los, his, study, group, method, region_name) {
  N <- sum(ns)
  obs <- lapply(seq_along(ps), function(i) binary_observed(ps[i]) * ns[i] / N)
  names(obs) <- names(ps)
  cons <- lapply(seq_along(ps), function(i) poc_constraint(poc_query(c(none = "event")), max(0, los[i]),
                 min(1, his[i]), stratum = names(ps)[i], label = paste0("q_", names(ps)[i])))
  model <- poc_partial_model(obs, cons, lp_backend = BACKEND)
  pns <- partial(model, Q_PNS); pn <- partial(model, Q_PN)
  data.frame(study = study, group = group, method = method, region = region_name, articles = N,
             years = paste(range(as.integer(names(ps))), collapse = "-"),
             PNS_lower = pns[["lower"]], PNS_upper = pns[["upper"]], PN_lower = pn[["lower"]], PN_upper = pn[["upper"]],
             opened_lower = pns[["lower"]] * N, opened_upper = pns[["upper"]] * N, stringsAsFactors = FALSE)
}
for (g in names(GROUPS)) for (m in c("a", "b")) {
  s <- agreement_inputs[agreement_inputs$group == g, ]
  ps <- setNames(s$p_licensed, s$year); qs <- s[[paste0("q", m, "_licensed")]]; d <- s[[paste0("margin_", m)]]
  pooled_rows[[length(pooled_rows) + 1L]] <- pool_strata(ps, s$n, qs, qs, "agreements", g, m, "point")
  pooled_rows[[length(pooled_rows) + 1L]] <- pool_strata(ps, s$n, qs - d, qs + d, "agreements", g, m, "margin")
  ok <- !is.na(s$attracted)
  f <- ifelse(ok, s$attracted, 0)
  pooled_rows[[length(pooled_rows) + 1L]] <- pool_strata(ps[ok], s$n[ok], ((1 - f) * (qs - d))[ok],
    (qs + d + f * (1 - qs))[ok], "agreements", g, m, "margin_choice")
}

# Checks fixed in advance: first-year jump, foreign comparison, lapse.
ag_checks <- list()
for (g in names(GROUPS)) {
  first <- min(GROUPS[[g]]$years); s <- agreement_inputs[agreement_inputs$group == g, ]
  jump <- s$p_licensed[s$year == first] - s$qa_licensed[s$year == first]
  ag_checks[[length(ag_checks) + 1L]] <- data.frame(check = "S1 first-year jump at least three margins", group = g,
    period = first, value_1 = jump, value_2 = 3 * margins[[g]][["a"]], passed = jump >= 3 * margins[[g]][["a"]])
  for (i in seq_len(nrow(s)))
    ag_checks[[length(ag_checks) + 1L]] <- data.frame(check = "S3 Polish change exceeds foreign change", group = g,
      period = s$year[i], value_1 = s$p_licensed[i] - s$qa_licensed[i], value_2 = s$foreign_change[i],
      passed = s$p_licensed[i] - s$qa_licensed[i] > s$foreign_change[i])
}
for (lapse in list(list("CUP", 2021:2023), list("OUP", 2021:2022), list("IEEE", 2023))) {
  pub <- lapse[[1]]; on <- intersect(agreement_years(pub), 2016:2023)
  d_on <- shares(counts(pub, "PL", on)); d_off <- shares(counts(pub, "PL", lapse[[2]]))
  d_pre <- shares(counts(pub, "PL", 2016:2019))
  ag_checks[[length(ag_checks) + 1L]] <- data.frame(check = paste("S2 lapse", pub), group = "2020",
    period = paste0(paste(range(on), collapse = "-"), " vs ", paste(range(lapse[[2]]), collapse = "-"),
                    "; n=", sum(counts(pub, "PL", on)), "/", sum(counts(pub, "PL", lapse[[2]]))),
    value_1 = d_on[["licensed"]], value_2 = d_off[["licensed"]], passed = d_off[["licensed"]] < d_on[["licensed"]])
  ag_checks[[length(ag_checks) + 1L]] <- data.frame(check = paste("S2 pre-agreement level", pub), group = "2020",
    period = "2016-2019", value_1 = d_pre[["licensed"]], value_2 = NA, passed = NA)
}
write_tsv(do.call(rbind, ag_checks), "agreement_checks.tsv")

# =============================================================== NCN mandate
K5 <- c("closed", "green", "bronze", "ccby", "other_licence")
mc <- read_tsv(file.path(agg, "mandate_counts.tsv"))
standardize <- function(rows0, rows1) {
  w1 <- stats::aggregate(n ~ panel + agreement_covered, rows1, sum); n1 <- sum(rows1$n)
  q <- setNames(numeric(length(K5)), K5); coef <- setNames(numeric(nrow(rows0)),
    paste(rows0$panel, rows0$agreement_covered)); fallback <- 0
  for (i in seq_len(nrow(w1))) {
    use <- rows0$panel == w1$panel[i] & rows0$agreement_covered == w1$agreement_covered[i]
    if (sum(rows0$n[use]) < 10) {
      use <- rows0$agreement_covered == w1$agreement_covered[i]
      if (!sum(rows0$n[use])) use <- rep(TRUE, nrow(rows0))
      fallback <- fallback + w1$n[i]
    }
    weight <- w1$n[i] / n1
    q <- q + weight * colSums(rows0[use, K5, drop = FALSE]) / sum(rows0$n[use])
    coef[use] <- coef[use] + weight / sum(rows0$n[use])
  }
  list(q = q, v_terms = list(n = rows0$n, coef = coef), fallback = fallback)
}
years5 <- 2022:2025; m_inputs <- list(); per <- list()
for (t in years5) {
  r1 <- mc[mc$year == t & mc$cohort == "M", ]; r0 <- mc[mc$year == t & mc$cohort == "C", ]
  rp <- mc[mc$year == t & mc$cohort == "P", ]
  p <- colSums(r1[, K5]) / sum(r1$n)
  sa <- standardize(r0, r1); sp <- standardize(rp, r1)
  raw <- 2 * sa$q - sp$q; qb <- pmin(pmax(raw, 0), 1); qb <- qb / sum(qb)
  va <- sum(sa$v_terms$n * sa$v_terms$coef^2)
  vb <- sum(sa$v_terms$n * (2 * sa$v_terms$coef)^2) + sum(sp$v_terms$n * sp$v_terms$coef^2)
  per[[as.character(t)]] <- list(n = sum(r1$n), x = sum(r1$ccby), p = p, q = list(a = sa$q, b = qb),
                                 placebo = sp$q, v = c(a = va, b = vb), fallback = sa$fallback,
                                 n0 = sum(r0$n), np = sum(rp$n), clipped = any(raw < 0 | raw > 1))
}
m_margin <- max(vapply(per, function(x) abs(x$q$a[["ccby"]] - x$placebo[["ccby"]]), 0))
m_bounds <- m_flows <- m_mono <- m_regions <- m_breakdown <- list()
for (t in years5) {
  r <- per[[as.character(t)]]
  m_inputs[[length(m_inputs) + 1L]] <- data.frame(year = t, n = r$n, n_comparison = r$n0, n_placebo = r$np,
    ccby_count = r$x, fallback_articles = r$fallback,
    setNames(as.list(r$p), paste0("p_", K5)), setNames(as.list(r$q$a), paste0("qa_", K5)),
    setNames(as.list(r$q$b), paste0("qb_", K5)), placebo_ccby = r$placebo[["ccby"]],
    margin = m_margin, qb_clipped = r$clipped, check.names = FALSE)
  p_ci <- clopper(r$x, r$n, 1 - ALPHA / 2)
  for (m in c("a", "b")) {
    lab <- data.frame(study = "mandate", group = "NCN", year = t, method = m, stringsAsFactors = FALSE)
    q <- r$q[[m]]; qL <- q[["ccby"]]; pL <- r$p[["ccby"]]
    a <- attribution(r$p, q, "ccby", lab)
    m_bounds[[length(m_bounds) + 1L]] <- a$bounds; m_flows[[length(m_flows) + 1L]] <- a$flows
    m_mono[[length(m_mono) + 1L]] <- a$monotone
    eps <- hoeffding(r$v[[m]], 1 - ALPHA / 2)
    reg <- rbind(point = region(pL, pL, qL, qL), margin = region(pL, pL, qL - m_margin, qL + m_margin),
                 margin_choice = rep(NA, 4),
                 sampling = region(p_ci[["lo"]], p_ci[["hi"]], qL - eps, qL + eps),
                 sampling_margin_choice = region(p_ci[["lo"]], p_ci[["hi"]], qL - eps - m_margin, qL + eps + m_margin))
    m_regions[[length(m_regions) + 1L]] <- data.frame(lab, region = rownames(reg), reg, row.names = NULL,
                                                      stringsAsFactors = FALSE)
    m_breakdown[[length(m_breakdown) + 1L]] <- data.frame(lab, margin = m_margin, breakdown_margin = max(0, pL - qL),
      ratio = max(0, pL - qL) / m_margin, attracted = NA, breakdown_attracted = NA,
      p_lo = p_ci[["lo"]], p_hi = p_ci[["hi"]], q_halfwidth = eps, stringsAsFactors = FALSE)
  }
}
mandate_inputs <- do.call(rbind, m_inputs)
write_tsv(mandate_inputs, "mandate_inputs.tsv")
for (m in c("a", "b")) {
  ps <- setNames(mandate_inputs$p_ccby, mandate_inputs$year); qs <- mandate_inputs[[paste0("q", m, "_ccby")]]
  pooled_rows[[length(pooled_rows) + 1L]] <- pool_strata(ps, mandate_inputs$n, qs, qs, "mandate", "NCN", m, "point")
  pooled_rows[[length(pooled_rows) + 1L]] <- pool_strata(ps, mandate_inputs$n, qs - m_margin, qs + m_margin,
                                                         "mandate", "NCN", m, "margin")
}

# Mandate checks: OPUS discontinuity, agreement layer, licence choice.
opus <- read_tsv(file.path(agg, "mandate_opus.tsv")); sh <- setNames(opus$ccby / opus$n, opus$call)
d1 <- sh[["OPUS 20"]] - sh[["OPUS 19"]]; d0 <- sh[["OPUS 19"]] - sh[["OPUS 18"]]
se1 <- sqrt(sum(sh[c("OPUS 20", "OPUS 19")] * (1 - sh[c("OPUS 20", "OPUS 19")]) / opus$n[match(c("OPUS 20", "OPUS 19"), opus$call)]))
m_checks <- list(data.frame(check = "S1 OPUS 20 minus OPUS 19 exceeds twice the OPUS 19 minus OPUS 18 change",
  value_1 = d1, value_2 = d0, value_3 = se1, passed = d1 > 0 && abs(d1) >= 2 * abs(d0)))
layer_std <- function(rows0, rows1) {
  w <- stats::aggregate(n ~ year + panel, rows1, sum); total <- sum(rows1$n); q <- 0
  for (i in seq_len(nrow(w))) {
    use <- rows0$year == w$year[i] & rows0$panel == w$panel[i]
    if (sum(rows0$n[use]) < 10) use <- rows0$year == w$year[i]
    q <- q + w$n[i] / total * sum(rows0$ccby[use]) / sum(rows0$n[use])
  }
  q
}
for (layer in c(0, 1)) {
  r1 <- mc[mc$cohort == "M" & mc$agreement_covered == layer, ]; r0 <- mc[mc$cohort == "C" & mc$agreement_covered == layer, ]
  rp <- mc[mc$cohort == "P" & mc$agreement_covered == layer, ]
  m_checks[[length(m_checks) + 1L]] <- data.frame(
    check = paste("S2 layer", if (layer) "covered by an agreement" else "outside agreements"),
    value_1 = sum(r1$ccby) / sum(r1$n), value_2 = layer_std(r0, r1), value_3 = layer_std(rp, r1), passed = NA)
}
lic <- read_tsv(file.path(agg, "mandate_licences.tsv"))
for (c in c("M", "C", "P")) {
  s <- lic[lic$cohort == c, ]
  m_checks[[length(m_checks) + 1L]] <- data.frame(check = paste("S3 CC BY share among licensed-open, cohort", c),
    value_1 = sum(s$n[s$licence == "cc-by"]) / s$n_open_at_publisher[1], value_2 = s$n_open_at_publisher[1],
    value_3 = NA, passed = NA)
}
m_checks <- do.call(rbind, m_checks)
s2 <- m_checks[grepl("^S2", m_checks$check), ]
m_checks$passed[grepl("^S2", m_checks$check)] <- (s2$value_1[1] - s2$value_2[1]) > (s2$value_1[2] - s2$value_2[2])
s3 <- m_checks[grepl("^S3", m_checks$check), ]
m_checks$passed[grepl("^S3", m_checks$check)] <- s3$value_1[1] > s3$value_1[2]
write_tsv(m_checks, "mandate_checks.tsv")

# ======================================================= combined outputs
bounds_all <- do.call(rbind, c(bounds, m_bounds)); write_tsv(bounds_all, "attribution_bounds.tsv")
flows_all <- do.call(rbind, c(flows, m_flows)); write_tsv(flows_all, "flow_bounds.tsv")
COMMON <- c("study", "group", "year", "method", "binary_compatible", "PNS", "PNS_upper", "PN", "PN_upper",
            "categories_compatible")
mono_common <- rbind(do.call(rbind, mono)[, COMMON], do.call(rbind, m_mono)[, COMMON])
write_tsv(do.call(rbind, c(mono)), "agreement_monotone.tsv")
write_tsv(do.call(rbind, c(m_mono)), "mandate_monotone.tsv")
regions_all <- do.call(rbind, c(regions, m_regions)); write_tsv(regions_all, "attribution_regions.tsv")
breakdown_all <- do.call(rbind, c(breakdown, m_breakdown)); write_tsv(breakdown_all, "breakdown.tsv")
# Pooled sampling projection: pooled exposed share (Clopper-Pearson) and the
# article-weighted counterfactual share (Hoeffding over all comparison articles).
pooled_sampling <- function(study, group, method, n, x, qs, vs, ds, fs) {
  N <- sum(n); w <- n / N; p_ci <- clopper(sum(x), N, 1 - ALPHA / 2)
  q <- sum(w * qs); eps <- hoeffding(sum(w^2 * vs), 1 - ALPHA / 2)
  bands <- list(sampling = c(q - eps, q + eps), sampling_margin = c(q - eps - sum(w * ds), q + eps + sum(w * ds)))
  if (!all(is.na(fs))) bands$sampling_margin_choice <-
    c(sum(w * (1 - fs) * (qs - ds)) - eps, sum(w * (qs + ds + fs * (1 - qs))) + eps)
  do.call(rbind, lapply(names(bands), function(b) {
    r <- region(p_ci[["lo"]], p_ci[["hi"]], bands[[b]][1], bands[[b]][2])
    data.frame(study = study, group = group, method = method, region = b, articles = N,
               years = "", PNS_lower = r[["PNS_lower"]], PNS_upper = r[["PNS_upper"]], PN_lower = r[["PN_lower"]],
               PN_upper = r[["PN_upper"]], opened_lower = r[["PNS_lower"]] * N, opened_upper = r[["PNS_upper"]] * N,
               stringsAsFactors = FALSE)
  }))
}
vrows <- do.call(rbind, breakdown)
for (g in names(GROUPS)) for (m in c("a", "b")) {
  s <- agreement_inputs[agreement_inputs$group == g, ]
  v <- vrows[vrows$group == g & vrows$method == m, ]
  vs <- (v$q_halfwidth^2) * 2 / log(2 / (ALPHA / 2))
  ok <- if (g == "2019") s$year >= 2020 else rep(TRUE, nrow(s))
  pr <- pooled_sampling("agreements", g, m, s$n, s$licensed_count, s[[paste0("q", m, "_licensed")]], vs,
                        s[[paste0("margin_", m)]], rep(NA, nrow(s)))
  pc <- pooled_sampling("agreements", g, m, s$n[ok], s$licensed_count[ok], s[[paste0("q", m, "_licensed")]][ok], vs[ok],
                        s[[paste0("margin_", m)]][ok], s$attracted[ok])
  pr$years <- paste(range(s$year), collapse = "-"); pc$years <- paste(range(s$year[ok]), collapse = "-")
  pooled_rows[[length(pooled_rows) + 1L]] <- pr
  pooled_rows[[length(pooled_rows) + 1L]] <- pc[pc$region == "sampling_margin_choice", ]
}
mv <- do.call(rbind, m_breakdown)
for (m in c("a", "b")) {
  v <- mv[mv$method == m, ]; vs <- (v$q_halfwidth^2) * 2 / log(2 / (ALPHA / 2))
  pr <- pooled_sampling("mandate", "NCN", m, mandate_inputs$n, mandate_inputs$ccby_count,
                        mandate_inputs[[paste0("q", m, "_ccby")]], vs, rep(m_margin, nrow(mandate_inputs)), NA)
  pr$years <- "2022-2025"; pooled_rows[[length(pooled_rows) + 1L]] <- pr
}
pooled_all <- do.call(rbind, pooled_rows); write_tsv(pooled_all, "pooled.tsv")

# ======================================================= verification
record("Li-Pearl and exact LP agree for PNS and PN",
       max(abs(c(bounds_all$exact_lower - bounds_all$lipearl_lower, bounds_all$exact_upper - bounds_all$lipearl_upper))) < 1e-8)
record("Shu-Wang-Li comparator agrees where it applies",
       all(is.na(bounds_all$shu_lower) | abs(bounds_all$shu_lower - bounds_all$exact_lower) < 1e-8 &
             abs(bounds_all$shu_upper - bounds_all$exact_upper) < 1e-8))
pns <- bounds_all[bounds_all$estimand == "PNS", ]
key_inputs <- rbind(data.frame(study = "agreements", group = agreement_inputs$group, year = agreement_inputs$year,
                               p = agreement_inputs$p_licensed, qa = agreement_inputs$qa_licensed, qb = agreement_inputs$qb_licensed),
                    data.frame(study = "mandate", group = "NCN", year = mandate_inputs$year,
                               p = mandate_inputs$p_ccby, qa = mandate_inputs$qa_ccby, qb = mandate_inputs$qb_ccby))
merged <- merge(pns, key_inputs, by = c("study", "group", "year"))
qq <- ifelse(merged$method == "a", merged$qa, merged$qb)
record("Exact PNS equals max(0, p - q) and min(p, 1 - q)",
       max(abs(merged$exact_lower - pmax(0, merged$p - qq)), abs(merged$exact_upper - pmin(merged$p, 1 - qq))) < 1e-8)
pt <- regions_all[regions_all$region == "point", ]
record("Point regions reproduce the complete-margin intervals",
       max(abs(pt$PNS_lower - pns$exact_lower[match(paste(pt$study, pt$group, pt$year, pt$method),
           paste(pns$study, pns$group, pns$year, pns$method))])) < 1e-8)
record("Regions are nested: point within margin within margin-and-choice", {
  ok <- TRUE
  for (key in unique(paste(regions_all$study, regions_all$group, regions_all$year, regions_all$method))) {
    s <- regions_all[paste(regions_all$study, regions_all$group, regions_all$year, regions_all$method) == key, ]
    rr <- function(n) unlist(s[s$region == n, c("PNS_lower", "PNS_upper", "PN_lower", "PN_upper")])
    inside <- function(a, b) all(is.na(b)) || (a[1] >= b[1] - TOL && a[2] <= b[2] + TOL && a[3] >= b[3] - TOL && a[4] <= b[4] + TOL)
    ok <- ok && inside(rr("point"), rr("margin")) && inside(rr("margin"), rr("margin_choice")) &&
      inside(rr("point"), rr("sampling")) && inside(rr("sampling"), rr("sampling_margin_choice"))
  }
  ok })
record("Monotone PNS equals p - q wherever the binary monotone model is compatible", {
  mm <- merge(mono_common, key_inputs, by = c("study", "group", "year"))
  qv <- ifelse(mm$method == "a", mm$qa, mm$qb); ok <- mm$binary_compatible
  all(ok == (mm$p >= qv - TOL)) && max(abs(mm$PNS[ok] - (mm$p - qv)[ok]), abs(mm$PNS_upper[ok] - mm$PNS[ok])) < 1e-8 })
record("Pooled strata point results equal article-weighted sums of yearly bounds", {
  ok <- TRUE
  for (i in which(pooled_all$region == "point")) {
    s <- pns[pns$study == pooled_all$study[i] & pns$group == pooled_all$group[i] & pns$method == pooled_all$method[i], ]
    ns <- if (pooled_all$study[i] == "agreements") agreement_inputs$n[agreement_inputs$group == pooled_all$group[i]] else mandate_inputs$n
    ok <- ok && abs(sum(ns * s$exact_lower) / sum(ns) - pooled_all$PNS_lower[i]) < 1e-8 &&
      abs(sum(ns * s$exact_upper) / sum(ns) - pooled_all$PNS_upper[i]) < 1e-8
  }
  ok })

write_tsv(checks, "checks.tsv")

# ======================================================= figures
# Legend whose columns are as wide as the longest label in the current font.
legend_text <- function(..., legend, cex) {
  graphics::legend(..., legend = legend, cex = cex, text.width = max(graphics::strwidth(legend, cex = cex)))
}
fig_data <- list()
palette <- c(a = "black", b = "grey45")
year_axis <- function(years) graphics::axis(1, at = years, labels = years, cex.axis = .8, mgp = c(2, .4, 0), tcl = -.25)
draw_intervals <- function(study, group, estimand, years, title) {
  graphics::plot(NA, xlim = range(years) + c(-.5, .5), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                 main = title, cex.main = .85, font.main = 1)
  year_axis(years)
  graphics::axis(2, at = seq(0, 1, .25), labels = c("0", ".25", ".50", ".75", "1"), las = 1, cex.axis = .8,
                 mgp = c(2, .5, 0), tcl = -.25)
  graphics::abline(h = seq(0, 1, .25), col = "grey90", lwd = .6)
  for (m in c("a", "b")) {
    off <- if (m == "a") -.14 else .14
    for (t in years) {
      s <- regions_all[regions_all$study == study & regions_all$group == group & regions_all$year == t &
                       regions_all$method == m, ]
      inner <- unlist(s[s$region == "point", paste0(estimand, c("_lower", "_upper"))])
      outer_name <- if (study == "agreements" && !is.na(s$PNS_lower[s$region == "margin_choice"])) "margin_choice" else "margin"
      outer <- unlist(s[s$region == outer_name, paste0(estimand, c("_lower", "_upper"))])
      graphics::segments(t + off, outer[1], t + off, outer[2], col = palette[[m]], lwd = .9)
      graphics::segments(t + off, inner[1], t + off, inner[2], col = palette[[m]], lwd = 4.5, lend = 1)
      v <- mono_common[mono_common$study == study & mono_common$group == group & mono_common$year == t &
                       mono_common$method == m, ]
      if (nrow(v) && isTRUE(v$binary_compatible))
        graphics::points(t + off, v[[estimand]], pch = 21, bg = "white", col = palette[[m]], cex = .75)
      fig_data[[length(fig_data) + 1L]] <<- data.frame(study = study, group = group, year = t, method = m,
        estimand = estimand, point_lower = inner[1], point_upper = inner[2], outer_region = outer_name,
        outer_lower = outer[1], outer_upper = outer[2],
        monotone = if (nrow(v) && isTRUE(v$binary_compatible)) v[[estimand]] else NA)
    }
  }
}
render_intervals <- function(file, device) {
  if (device == "pdf") grDevices::pdf(file, width = 7, height = 5.6, family = "Helvetica") else
    grDevices::png(file, width = 7, height = 5.6, units = "in", res = 300, family = "Helvetica")
  graphics::par(mfrow = c(2, 3), mar = c(2.2, 2.6, 1.6, .4), oma = c(4.4, 1.4, 0, 0))
  for (est in c("PNS", "PN")) {
    draw_intervals("agreements", "2020", est, 2020:2023, paste0(est, ": agreements from 2020"))
    draw_intervals("agreements", "2019", est, 2019:2023, paste0(est, ": agreements from 2019"))
    draw_intervals("mandate", "NCN", est, 2022:2025, paste0(est, ": NCN mandate"))
  }
  graphics::mtext("Publication year", side = 1, outer = TRUE, line = .5, cex = .75)
  graphics::mtext("Probability", side = 2, outer = TRUE, line = .2, cex = .75)
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot.new()
  legend_text("bottom", ncol = 2, bty = "n", cex = .72, inset = c(0, .004), x.intersp = .8,
    legend = c("Method (a): level before the policy", "Method (b): comparison trend",
               "Point value under monotonicity", "Identification interval at point inputs",
               "Assumption region (margin; publisher choice)", ""),
    col = c("black", "grey45", "black", "black", "black", NA), lwd = c(4.5, 4.5, NA, 4.5, .9, NA),
    pch = c(NA, NA, 21, NA, NA, NA), pt.bg = "white")
  grDevices::dev.off()
}
fig_data <- list()
render_intervals(file.path(out, "fig-attribution.pdf"), "pdf")
write_tsv(do.call(rbind, fig_data), "fig-attribution.data.tsv")
render_intervals(file.path(out, "fig-attribution.png"), "png")

# Design figure: observed licensed shares and counterfactual estimates.
series <- list()
for (g in names(GROUPS)) {
  pubs <- GROUPS[[g]]$publishers
  for (corr in c("PL", "foreign")) for (y in 2016:2023) {
    if (corr == "PL" && y %in% GROUPS[[g]]$years) next
    v <- Reduce(`+`, Filter(Negate(is.null), lapply(pubs, function(pub) counts(pub, corr, y))))
    series[[length(series) + 1L]] <- data.frame(panel = g, series = paste0(corr, "_observed"), year = y,
                                                value = v[["licensed"]] / sum(v), n = sum(v))
  }
  s <- agreement_inputs[agreement_inputs$group == g, ]
  series[[length(series) + 1L]] <- data.frame(panel = g, series = "PL_covered", year = s$year, value = s$p_licensed, n = s$n)
  series[[length(series) + 1L]] <- data.frame(panel = g, series = "counterfactual_a", year = s$year, value = s$qa_licensed, n = NA)
  series[[length(series) + 1L]] <- data.frame(panel = g, series = "counterfactual_b", year = s$year, value = s$qb_licensed, n = NA)
}
series[[length(series) + 1L]] <- data.frame(panel = "NCN", series = "PL_covered", year = mandate_inputs$year, value = mandate_inputs$p_ccby, n = mandate_inputs$n)
series[[length(series) + 1L]] <- data.frame(panel = "NCN", series = "counterfactual_a", year = mandate_inputs$year, value = mandate_inputs$qa_ccby, n = NA)
series[[length(series) + 1L]] <- data.frame(panel = "NCN", series = "counterfactual_b", year = mandate_inputs$year, value = mandate_inputs$qb_ccby, n = NA)
series[[length(series) + 1L]] <- data.frame(panel = "NCN", series = "placebo", year = mandate_inputs$year, value = mandate_inputs$placebo_ccby, n = NA)
series <- do.call(rbind, series)
write_tsv(series, "fig-design.data.tsv")
render_design <- function(file, device) {
  if (device == "pdf") grDevices::pdf(file, width = 7, height = 3.2, family = "Helvetica") else
    grDevices::png(file, width = 7, height = 3.2, units = "in", res = 300, family = "Helvetica")
  graphics::par(mfrow = c(1, 3), mar = c(2.2, 2.6, 1.6, .4), oma = c(4, 1.4, 0, 0))
  panel <- function(key, years, title, start) {
    graphics::plot(NA, xlim = range(years) + c(-.3, .3), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                   main = title, cex.main = .85, font.main = 1)
    year_axis(years)
    graphics::axis(2, at = seq(0, 1, .25), labels = c("0", ".25", ".50", ".75", "1"), las = 1, cex.axis = .8,
                   mgp = c(2, .5, 0), tcl = -.25)
    graphics::abline(h = seq(0, 1, .25), col = "grey90", lwd = .6)
    if (!is.na(start)) graphics::abline(v = start - .5, lty = 3, col = "grey40")
    s <- series[series$panel == key, ]
    line <- function(name, ...) { z <- s[s$series == name, ]; z <- z[order(z$year), ]
      if (nrow(z)) graphics::lines(z$year, z$value, ...) }
    if (key != "NCN") {
      line("PL_observed", type = "o", pch = 19, cex = .6, lwd = 1.3)
      line("PL_covered", type = "o", pch = 19, cex = .6, lwd = 1.3)
      line("foreign_observed", type = "o", pch = 1, cex = .6, lty = 2, col = "grey30")
    } else line("PL_covered", type = "o", pch = 19, cex = .6, lwd = 1.3)
    line("counterfactual_a", type = "o", pch = 22, cex = .7, lty = 1, col = "grey35", bg = "white")
    line("counterfactual_b", type = "o", pch = 24, cex = .7, lty = 1, col = "grey55", bg = "white")
    if (key == "NCN") line("placebo", type = "o", pch = 4, cex = .7, lty = 3, col = "grey40")
  }
  panel("2020", 2016:2023, "Agreements from 2020", 2020)
  panel("2019", 2016:2023, "Agreements from 2019", 2019)
  panel("NCN", 2022:2025, "NCN mandate (CC BY)", NA)
  graphics::mtext("Publication year", side = 1, outer = TRUE, line = .4, cex = .75)
  graphics::mtext("Share of articles", side = 2, outer = TRUE, line = .2, cex = .75)
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot.new()
  legend_text("bottom", ncol = 3, bty = "n", cex = .7, inset = c(0, .004),
    legend = c("Polish corresponding author", "Foreign corresponding author",
               "Counterfactual (a)", "Counterfactual (b)", "Earlier grants, standardized"),
    lty = c(1, 2, 1, 1, 3), pch = c(19, 1, 22, 24, 4), col = c("black", "grey30", "grey35", "grey55", "grey40"), pt.bg = "white")
  grDevices::dev.off()
}
render_design(file.path(out, "fig-design.pdf"), "pdf"); render_design(file.path(out, "fig-design.png"), "png")

# Flow figure: decomposition of the licensed share under monotonicity, method (a).
flow_rows <- list()
for (g in c("2020", "2019")) {
  s <- agreement_inputs[agreement_inputs$group == g, ]
  for (i in seq_len(nrow(s))) {
    mo <- do.call(rbind, mono); v <- mo[mo$group == g & mo$year == s$year[i] & mo$method == "a", ]
    flow_rows[[length(flow_rows) + 1L]] <- data.frame(panel = g, year = s$year[i], compatible = v$categories_compatible,
      retained = s$qa_licensed[i], closed = s$qa_closed[i] - s$p_closed[i], green = s$qa_green[i] - s$p_green[i],
      bronze = s$qa_bronze[i] - s$p_bronze[i], other_licence = NA, total = s$p_licensed[i])
  }
}
for (i in seq_len(nrow(mandate_inputs))) {
  v <- do.call(rbind, m_mono); v <- v[v$year == mandate_inputs$year[i] & v$method == "a", ]
  mi <- mandate_inputs[i, ]
  flow_rows[[length(flow_rows) + 1L]] <- data.frame(panel = "NCN", year = mi$year, compatible = v$categories_compatible,
    retained = mi$qa_ccby, closed = mi$qa_closed - mi$p_closed, green = mi$qa_green - mi$p_green,
    bronze = mi$qa_bronze - mi$p_bronze, other_licence = mi$qa_other_licence - mi$p_other_licence, total = mi$p_ccby)
}
flow_rows <- do.call(rbind, flow_rows)
write_tsv(flow_rows, "fig-flows.data.tsv")
render_flows <- function(file, device) {
  if (device == "pdf") grDevices::pdf(file, width = 7, height = 3.2, family = "Helvetica") else
    grDevices::png(file, width = 7, height = 3.2, units = "in", res = 300, family = "Helvetica")
  graphics::par(mfrow = c(1, 3), mar = c(2.2, 2.6, 1.6, .4), oma = c(4, 1.4, 0, 0))
  fills <- c(retained = "grey15", closed = "grey45", green = "grey72", bronze = "grey90", other_licence = "white")
  for (key in c("2020", "2019", "NCN")) {
    s <- flow_rows[flow_rows$panel == key, ]
    graphics::plot(NA, xlim = range(s$year) + c(-.6, .6), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                   main = c("2020" = "Agreements from 2020", "2019" = "Agreements from 2019", NCN = "NCN mandate (CC BY)")[[key]],
                   cex.main = .85, font.main = 1)
    year_axis(s$year)
    graphics::axis(2, at = seq(0, 1, .25), labels = c("0", ".25", ".50", ".75", "1"), las = 1, cex.axis = .8,
                   mgp = c(2, .5, 0), tcl = -.25)
    for (i in seq_len(nrow(s))) {
      parts <- c(retained = s$retained[i], closed = s$closed[i], green = s$green[i], bronze = s$bronze[i],
                 other_licence = if (is.na(s$other_licence[i])) 0 else s$other_licence[i])
      low <- 0
      for (k in names(parts)) {
        h <- max(0, parts[[k]])
        graphics::rect(s$year[i] - .3, low, s$year[i] + .3, low + h, col = fills[[k]], border = "grey20", lwd = .5)
        low <- low + h
      }
      if (!isTRUE(s$compatible[i])) graphics::text(s$year[i], min(1, low + .05), "*", cex = 1.1)
    }
  }
  graphics::mtext("Publication year", side = 1, outer = TRUE, line = .4, cex = .75)
  graphics::mtext("Share of exposed articles", side = 2, outer = TRUE, line = .2, cex = .75)
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot.new()
  legend_text("bottom", ncol = 3, bty = "n", cex = .7, inset = c(0, .004),
    legend = c("Event also without the policy", "Opened from closed", "Opened from green",
               "Opened from bronze", "From another license (mandate)"),
    fill = fills, border = "grey20")
  grDevices::dev.off()
}
render_flows(file.path(out, "fig-flows.pdf"), "pdf"); render_flows(file.path(out, "fig-flows.png"), "png")

# ======================================================= manifest
writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
files <- sort(list.files(out))
manifest <- list(script = "analysis.R", pnsbib = as.character(packageVersion("pnsbib")), backend = BACKEND,
                 alpha = ALPHA, checks_passed = all(checks$passed), n_checks = nrow(checks),
                 outputs_sha256 = setNames(as.list(unname(vapply(file.path(out, files), function(f)
                   digest::digest(file = f, algo = "sha256"), ""))), files),
                 inputs_sha256 = setNames(as.list(unname(vapply(list.files(agg, full.names = TRUE), function(f)
                   digest::digest(file = f, algo = "sha256"), ""))), file.path("aggregates", list.files(agg))))
writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, digits = NA), file.path(out, "manifest.json"))
print(checks)
if (!all(checks$passed)) quit(status = 1L)
