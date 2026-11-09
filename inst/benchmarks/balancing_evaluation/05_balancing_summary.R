# Targeted balancing evaluation: declared summaries of the five outcomes.
# Design: 0_materials/validation_protocol/balancing-evaluation-v1.md, section 5.
#
#   Rscript 2_syntax/05_balancing_summary.R --self-test
#   Rscript 2_syntax/05_balancing_summary.R --in <folder with shard files> --out <folder for tables>
#
# Every table is per cell; nothing is aggregated across cells. Base R only.

TOLERANCES <- c(0.10, 0.05, 0.15)
RMSE_MARGIN <- 0.02
BOOT <- 2000L
OUTLIER_SD <- 0.25
CONTAMINATION <- 0.10

slope_population <- function(p_b, b_a, b_b, q, s, sd_out = OUTLIER_SD) {
  w <- c(1 - p_b, p_b * (1 - q), p_b * q)
  mx <- c(0, 0, s)
  my <- c(0, 0, -s)
  vx <- c(1, 1, sd_out^2)
  cxy <- c(b_a, b_b, 0)
  (sum(w * (cxy + mx * my)) - sum(w * mx) * sum(w * my)) / (sum(w * (vx + mx^2)) - sum(w * mx)^2)
}

se <- function(z) sqrt(stats::var(z, na.rm = TRUE) / sum(!is.na(z)))

read_run <- function(dir) {
  cond <- utils::read.csv(file.path(dir, "conditions.csv"))
  read_all <- function(pattern) {
    f <- sort(list.files(dir, pattern, full.names = TRUE))
    if (!length(f)) stop("no files matching ", pattern, " in ", dir)
    do.call(rbind, lapply(f, utils::read.csv, stringsAsFactors = FALSE))
  }
  list(cells = cond, datasets = read_all("-datasets\\.csv\\.gz$"), repeats = read_all("-repeats\\.csv\\.gz$"))
}

# One row per dataset and method: mean change over successful repeats, plus M1 from the weighted fit.
dataset_means <- function(repeats, datasets) {
  ok <- repeats[repeats$status == "ok", ]
  key <- paste(ok$cell, ok$rep, ok$method)
  out <- data.frame(key = sort(unique(key)), stringsAsFactors = FALSE)
  out$dbar <- as.numeric(tapply(ok$change, key, mean)[out$key])
  out$sd_repeats <- as.numeric(tapply(ok$change, key, stats::sd)[out$key])
  out$n_ok <- as.integer(tapply(ok$change, key, length)[out$key])
  out$slope_gen_b <- as.numeric(tapply(ok$slope_generated_b, key, mean)[out$key])
  parts <- strsplit(out$key, " ", fixed = TRUE)
  out$cell <- as.integer(vapply(parts, `[`, "", 1L))
  out$rep <- as.integer(vapply(parts, `[`, "", 2L))
  out$method <- vapply(parts, `[`, "", 3L)
  out$key <- NULL
  w <- data.frame(dbar = datasets$change_weighted, sd_repeats = NA_real_, n_ok = 1L, slope_gen_b = NA_real_,
    cell = datasets$cell, rep = datasets$rep, method = "M1", stringsAsFactors = FALSE)
  rbind(out, w)
}

references <- function(cells) {
  clean <- cells[cells$severity == 0, ]
  a <- cells[cells$severity > 0, ]
  rbind(
    data.frame(cell = clean$cell, reference = "clean", delta_star = clean$delta_artifact, beta_b_star = clean$b_b),
    data.frame(cell = a$cell, reference = "artifact", delta_star = a$delta_artifact, beta_b_star = a$b_b),
    data.frame(cell = a$cell, reference = "genuine", delta_star = a$delta_genuine,
      beta_b_star = mapply(function(b, s) slope_population(1, 0, b, CONTAMINATION, s), a$b_b, a$severity)))
}

by_group <- function(d, groups, fun) {
  key <- do.call(paste, c(d[groups], sep = "\r"))
  pieces <- split(seq_len(nrow(d)), key)
  res <- lapply(pieces, function(i) cbind(d[i[1L], groups, drop = FALSE], fun(d[i, , drop = FALSE])))
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}

summarize_run <- function(run, boot = BOOT, tolerances = TOLERANCES) {
  cells <- run$cells
  means <- dataset_means(run$repeats, run$datasets)
  refs <- references(cells)
  scored <- merge(means, refs, by = "cell")
  scored$error <- scored$dbar - scored$delta_star

  # Completion, counting every expected repeat and every failure.
  completion <- by_group(run$repeats, c("cell", "method"), function(x) data.frame(
    expected = nrow(x), ok = sum(x$status == "ok"), failed = sum(x$status != "ok"),
    with_warning = sum(x$status == "ok" & !is.na(x$message) & nzchar(x$message)),
    wrong_row_count = sum(x$status == "ok" & x$n_returned != x$n_returned_b)))

  # Outcome 1: tracking.
  tracking <- by_group(scored, c("cell", "method", "reference"), function(x) data.frame(
    n_datasets = sum(!is.na(x$error)), delta_star = x$delta_star[1L], mean_change = mean(x$dbar, na.rm = TRUE),
    signed_error = mean(x$error, na.rm = TRUE), se_signed_error = se(x$error),
    rmse = sqrt(mean(x$error^2, na.rm = TRUE))))

  # Outcome 2: classification against each tolerance.
  classification <- do.call(rbind, lapply(tolerances, function(tol) {
    tab <- by_group(scored, c("cell", "method", "reference"), function(x) {
      material <- abs(x$delta_star[1L]) > tol
      flag <- mean(abs(x$dbar) > tol, na.rm = TRUE)
      data.frame(delta_star = x$delta_star[1L], truth_material = material, flag_rate = flag,
        false_alarm_rate = if (material) NA_real_ else flag, miss_rate = if (material) 1 - flag else NA_real_)
    })
    cbind(tolerance = tol, tab)
  }))

  # Outcome 3: influence, contaminated minus its paired clean dataset.
  info <- cells[, c("cell", "base_cell", "severity")]
  with_info <- merge(means, info, by = "cell")
  clean <- with_info[with_info$severity == 0, c("base_cell", "rep", "method", "dbar")]
  names(clean)[4L] <- "dbar_clean"
  paired <- merge(with_info[with_info$severity > 0, ], clean, by = c("base_cell", "rep", "method"))
  paired$influence <- paired$dbar - paired$dbar_clean
  influence <- by_group(paired, c("cell", "method"), function(x) data.frame(
    n_pairs = sum(!is.na(x$influence)), mean_influence = mean(x$influence, na.rm = TRUE),
    sd_influence = stats::sd(x$influence, na.rm = TRUE), se_influence = se(x$influence)))

  # Outcome 4: slope within generated group B rows against the reference group B slope.
  gen <- scored[scored$method != "M1", ]
  observed <- merge(run$datasets[, c("cell", "rep", "slope_b")], refs, by = "cell")
  retention <- rbind(
    by_group(gen, c("cell", "method", "reference"), function(x) data.frame(
      beta_b_star = x$beta_b_star[1L], mean_slope_generated_b = mean(x$slope_gen_b, na.rm = TRUE),
      signed_difference = mean(x$slope_gen_b - x$beta_b_star, na.rm = TRUE),
      se_difference = se(x$slope_gen_b - x$beta_b_star))),
    by_group(transform(observed, method = "observed"), c("cell", "method", "reference"), function(x) data.frame(
      beta_b_star = x$beta_b_star[1L], mean_slope_generated_b = mean(x$slope_b, na.rm = TRUE),
      signed_difference = mean(x$slope_b - x$beta_b_star, na.rm = TRUE),
      se_difference = se(x$slope_b - x$beta_b_star))))

  # Outcome 5: generator variability.
  ok <- run$repeats[run$repeats$status == "ok", ]
  joined <- merge(ok[, c("cell", "rep", "method", "change")], means[, c("cell", "rep", "method", "dbar")],
    by = c("cell", "rep", "method"))
  variability <- by_group(means[means$method != "M1", ], c("cell", "method"), function(x) {
    within <- mean(x$sd_repeats, na.rm = TRUE)
    between <- stats::sd(x$dbar, na.rm = TRUE)
    data.frame(sd_within_dataset = within, sd_between_datasets = between, ratio = within / between)
  })
  for (tol in tolerances) {
    name <- sprintf("single_run_disagrees_at_%.2f", tol)
    d <- by_group(joined, c("cell", "method"), function(x) {
      v <- mean((abs(x$change) > tol) != (abs(x$dbar) > tol), na.rm = TRUE)
      setNames(data.frame(v), name)
    })
    variability <- merge(variability, d, by = c("cell", "method"), all.x = TRUE)
  }

  # RMSE against weighting, paired over datasets; the accepted margin is 0.02.
  set.seed(20261109L)
  w <- scored[scored$method == "M1", c("cell", "rep", "reference", "error")]
  names(w)[4L] <- "error_w"
  vs <- merge(scored[scored$method != "M1", c("cell", "rep", "reference", "method", "error")], w,
    by = c("cell", "rep", "reference"))
  rmse_vs <- by_group(vs, c("cell", "method", "reference"), function(x) {
    keep <- !is.na(x$error) & !is.na(x$error_w)
    a <- x$error[keep]^2
    b <- x$error_w[keep]^2
    n <- length(a)
    diff <- sqrt(mean(a)) - sqrt(mean(b))
    idx <- matrix(sample.int(n, n * boot, replace = TRUE), nrow = n)
    boots <- sqrt(colMeans(matrix(a[idx], nrow = n))) - sqrt(colMeans(matrix(b[idx], nrow = n)))
    q <- stats::quantile(boots, c(0.025, 0.975), names = FALSE)
    data.frame(n_datasets = n, rmse_minus_weighting = diff, lower = q[1L], upper = q[2L],
      within_margin = diff <= RMSE_MARGIN)
  })

  list(completion = completion, tracking = tracking, classification = classification, influence = influence,
    retention = retention, variability = variability, rmse_vs_weighting = rmse_vs)
}

# ---- self-test on invented numbers with hand-calculated answers ----------------------------

self_test <- function() {
  ok <- TRUE
  check <- function(name, pass) {
    cat(sprintf("%-4s %s\n", if (isTRUE(pass)) "PASS" else "FAIL", name))
    if (!isTRUE(pass)) ok <<- FALSE
  }
  near <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-9))
  cells <- data.frame(cell = 1:2, base_cell = c(1L, 1L), severity = c(0, 2), share = 0.2, b_a = 0, b_b = 0.4,
    delta_artifact = 0.12, delta_genuine = c(0.12, -0.003))
  datasets <- data.frame(cell = rep(1:2, each = 2), base_cell = 1L, rep = rep(1:2, 2), t0 = 0,
    change_weighted = c(0.10, 0.14, 0.20, 0.24), sd_x = 1, sd_y = 1, slope_a = 0, slope_b = c(0.3, 0.5, 0.3, 0.5))
  grid <- expand.grid(cell = 1:2, rep = 1:2, repeat_id = 1:2)
  repeats <- data.frame(grid, base_cell = 1L, method = "M2", seed = 1L, status = "ok", n_returned = 60L,
    n_returned_b = 60L, elapsed = 0, message = "", stringsAsFactors = FALSE)
  repeats$change <- ifelse(repeats$repeat_id == 1, 0.1, 0.3) + ifelse(repeats$cell == 2, 0.3, 0)
  repeats$slope_generated_b <- 0.2
  failed <- repeats[1L, ]
  failed$repeat_id <- 3L
  failed$status <- "error"
  failed$change <- NA_real_
  failed$slope_generated_b <- NA_real_
  failed$n_returned <- NA_integer_
  failed$n_returned_b <- NA_integer_
  res <- summarize_run(list(cells = cells, datasets = datasets, repeats = rbind(repeats, failed)), boot = 50L)

  p <- function(t, ...) {
    f <- list(...)
    t[Reduce(`&`, Map(function(n, v) t[[n]] == v, names(f), f)), ]
  }
  x <- p(res$completion, cell = 1L, method = "M2")
  check("completion counts the failed repeat", x$expected == 5L && x$ok == 4L && x$failed == 1L)
  x <- p(res$tracking, cell = 1L, method = "M2", reference = "clean")
  check("tracking mean change 0.2", near(x$mean_change, 0.2))
  check("tracking signed error 0.08", near(x$signed_error, 0.08))
  check("tracking RMSE 0.08", near(x$rmse, 0.08))
  x <- p(res$tracking, cell = 2L, method = "M1", reference = "genuine")
  check("weighted change scored against the genuine reference", near(x$signed_error, 0.22 + 0.003))
  check("clean cell has no artifact row", nrow(p(res$tracking, cell = 1L, method = "M2", reference = "artifact")) == 0L)
  x <- p(res$classification, tolerance = 0.1, cell = 1L, method = "M2", reference = "clean")
  check("material truth gives a miss rate of 0", isTRUE(x$truth_material) && near(x$miss_rate, 0) && is.na(x$false_alarm_rate))
  x <- p(res$classification, tolerance = 0.1, cell = 2L, method = "M2", reference = "genuine")
  check("non-material truth gives a false alarm rate of 1", !x$truth_material && near(x$false_alarm_rate, 1) && is.na(x$miss_rate))
  x <- p(res$influence, cell = 2L, method = "M2")
  check("influence is the paired difference 0.3", near(x$mean_influence, 0.3) && x$n_pairs == 2L)
  x <- p(res$influence, cell = 2L, method = "M1")
  check("influence includes weighting at 0.1", near(x$mean_influence, 0.1))
  x <- p(res$retention, cell = 1L, method = "M2", reference = "clean")
  check("retention difference is -0.2", near(x$signed_difference, -0.2))
  x <- p(res$retention, cell = 2L, method = "M2", reference = "genuine")
  check("genuine group B slope uses the mixture", near(x$beta_b_star, slope_population(1, 0, 0.4, 0.1, 2)))
  x <- p(res$retention, cell = 1L, method = "observed", reference = "clean")
  check("observed group B slope is reported", near(x$mean_slope_generated_b, 0.4))
  x <- p(res$variability, cell = 1L, method = "M2")
  check("within-dataset SD", near(x$sd_within_dataset, sqrt(0.02)))
  check("one single run in two disagrees at 0.10", near(x$single_run_disagrees_at_0.10, 0.5))
  x <- p(res$rmse_vs_weighting, cell = 1L, method = "M2", reference = "clean")
  w_rmse <- sqrt(mean((c(0.10, 0.14) - 0.12)^2))
  check("RMSE minus weighting", near(x$rmse_minus_weighting, 0.08 - w_rmse))
  check("margin flag applied", identical(x$within_margin, 0.08 - w_rmse <= RMSE_MARGIN))
  if (!ok) quit(save = "no", status = 1)
  cat("All self-tests passed.\n")
}

# ---- main --------------------------------------------------------------------------------

argv <- commandArgs(trailingOnly = TRUE)
if ("--self-test" %in% argv) {
  self_test()
} else {
  get <- function(flag) {
    i <- match(flag, argv)
    if (is.na(i) || i == length(argv)) stop("missing ", flag)
    argv[[i + 1L]]
  }
  out <- get("--out")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  res <- summarize_run(read_run(get("--in")))
  for (n in names(res)) utils::write.csv(res[[n]], file.path(out, paste0(n, ".csv")), row.names = FALSE)
  cat("Wrote", length(res), "tables to", out, "\n")
}
