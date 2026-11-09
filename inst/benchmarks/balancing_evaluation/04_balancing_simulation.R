# Targeted balancing evaluation: simulation runner.
# Design: 0_materials/validation_protocol/balancing-evaluation-v1.md (settings accepted, Relay decision 017).
#
# Run from the project root:
#   Rscript 2_syntax/04_balancing_simulation.R --mode check
#   Rscript 2_syntax/04_balancing_simulation.R --mode pilot --reps 2 --repeats 2
#   Rscript 2_syntax/04_balancing_simulation.R --mode full --shard 1 --shards 20 --confirm-authorized-run
#
# check  software assertions on separate seeds; prints PASS or FAIL and timings, never an outcome.
# pilot  timing only; writes elapsed seconds per method and discards every estimate.
# full   the declared study. Refuses to run without --confirm-authorized-run and never overwrites a shard.

suppressPackageStartupMessages(library(sydni))

DESIGN <- list(
  n_total = c(100L, 200L, 400L),
  share = c(0.10, 0.20, 0.30),
  slopes = list(c(0.2, 0.2), c(0, 0.4)),
  severity = c(0, 2, 3),
  contamination_rate = 0.10,
  outlier_sd = 0.25,
  reps = 200L,
  repeats = 10L
)

METHODS <- data.frame(
  label = paste0("M", 2:9),
  index = 2:9,
  kind = c("oversample", "oversample_noise", rep("generator", 6)),
  engine = c(NA, NA, "graphical", "graphical", "graphical", "graphical", "copula", "bart"),
  noise = c("none", "default", "none", "default", "none", "default", "default", "default"),
  k = c(NA, NA, 0, 0, 2, 2, 2, 2),
  stringsAsFactors = FALSE
)

# ---- arguments ---------------------------------------------------------------

parse_args <- function(argv) {
  out <- list(mode = "check", shard = 1L, shards = 1L, reps = NA_integer_, repeats = NA_integer_,
    out = "3_results/04_balancing_simulation", confirm = FALSE)
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (a == "--confirm-authorized-run") {
      out$confirm <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(argv)) stop("missing value for ", a)
    v <- argv[[i + 1L]]
    switch(a,
      "--mode" = out$mode <- v,
      "--shard" = out$shard <- as.integer(v),
      "--shards" = out$shards <- as.integer(v),
      "--reps" = out$reps <- as.integer(v),
      "--repeats" = out$repeats <- as.integer(v),
      "--out" = out$out <- v,
      stop("unknown argument ", a))
    i <- i + 2L
  }
  if (!out$mode %in% c("check", "pilot", "full")) stop("--mode must be check, pilot or full")
  if (out$shard < 1L || out$shard > out$shards) stop("--shard must lie in 1..shards")
  out
}

# ---- conditions and population references ---------------------------------------

slope_population <- function(p_b, b_a, b_b, q, s, sd_out) {
  w <- c(1 - p_b, p_b * (1 - q), p_b * q)
  mx <- c(0, 0, s)
  my <- c(0, 0, -s)
  vx <- c(1, 1, sd_out^2)
  cxy <- c(b_a, b_b, 0)
  (sum(w * (cxy + mx * my)) - sum(w * mx) * sum(w * my)) / (sum(w * (vx + mx^2)) - sum(w * mx)^2)
}

build_conditions <- function(design) {
  base <- expand.grid(
    n_total = design$n_total, share = design$share, slope_set = seq_along(design$slopes),
    KEEP.OUT.ATTRS = FALSE)
  base$base_cell <- seq_len(nrow(base))
  cells <- merge(base, data.frame(severity = design$severity), all = TRUE)
  cells <- cells[order(cells$base_cell, cells$severity), ]
  cells$cell <- seq_len(nrow(cells))
  cells$b_a <- vapply(design$slopes[cells$slope_set], `[`, numeric(1), 1)
  cells$b_b <- vapply(design$slopes[cells$slope_set], `[`, numeric(1), 2)
  cells$n_b <- as.integer(round(cells$share * cells$n_total))
  cells$n_a <- cells$n_total - cells$n_b
  cells$k_out <- ifelse(cells$severity == 0, 0L,
    pmax(1L, as.integer(round(design$contamination_rate * cells$n_b))))
  ref <- function(p, q, s, a, b) slope_population(p, a, b, q, s, design$outlier_sd)
  cells$delta_artifact <- mapply(function(p, a, b) ref(0.5, 0, 0, a, b) - ref(p, 0, 0, a, b),
    cells$share, cells$b_a, cells$b_b)
  cells$delta_genuine <- mapply(function(p, a, b, s) {
    q <- if (s == 0) 0 else design$contamination_rate
    ref(0.5, q, s, a, b) - ref(p, q, s, a, b)
  }, cells$share, cells$b_a, cells$b_b, cells$severity)
  rownames(cells) <- NULL
  cells
}

# ---- random streams ---------------------------------------------------------------
# The base sample and outlier draws depend only on base_cell and rep, so clean and contaminated
# versions, and both severities, are paired. Generator streams are shared across severities too.

stream <- function(kind, base_cell, rep, method = 0L, repeat_id = 0L, offset = 0) {
  unit <- (base_cell - 1L) * 1000L + rep
  as.integer(switch(kind,
    data = 1e8 + unit,
    outlier = 2e8 + unit,
    generator = 3e8 + unit * 2000 + method * 100 + repeat_id) + offset)
}

# ---- data ----------------------------------------------------------------------------

make_data <- function(cell, rep, offset, design) {
  set.seed(stream("data", cell$base_cell, rep, offset = offset))
  g <- c(rep("A", cell$n_a), rep("B", cell$n_b))
  b <- ifelse(g == "A", cell$b_a, cell$b_b)
  x <- rnorm(cell$n_total)
  y <- b * x + sqrt(1 - b^2) * rnorm(cell$n_total)
  if (cell$k_out > 0L) {
    set.seed(stream("outlier", cell$base_cell, rep, offset = offset))
    idx <- seq.int(cell$n_a + 1L, length.out = cell$k_out)
    z_x <- rnorm(cell$k_out)
    z_y <- rnorm(cell$k_out)
    x[idx] <- cell$severity + design$outlier_sd * z_x
    y[idx] <- -cell$severity + design$outlier_sd * z_y
  }
  data.frame(group = factor(g, levels = c("A", "B")), x = x, y = y)
}

# ---- analysis ------------------------------------------------------------------------

slope <- function(x, y, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  mx <- sum(w * x) / sum(w)
  my <- sum(w * y) / sum(w)
  sum(w * (x - mx) * (y - my)) / sum(w * (x - mx)^2)
}

balance_weights <- function(group) {
  p <- as.numeric(table(group)[as.character(group)]) / length(group)
  0.5 / p
}

# ---- constructions -------------------------------------------------------------------

construct <- function(d, m, seed) {
  n_aug <- as.integer(sum(d$group == "A") - sum(d$group == "B"))
  if (m$kind %in% c("oversample", "oversample_noise")) {
    set.seed(seed)
    b_rows <- d[d$group == "B", , drop = FALSE]
    out <- b_rows[sample.int(nrow(b_rows), n_aug, replace = TRUE), , drop = FALSE]
    if (m$kind == "oversample_noise") out <- as.data.frame(add_noise(out, variables = c("x", "y")))
    rownames(out) <- NULL
    return(out)
  }
  set.seed(seed)
  args <- list(data = d, variables = c("x", "y"), covariates = "group", engine = m$engine,
    balance_by = "group", balance_value = "B", k = m$k, seed = seed)
  if (m$noise == "none") {
    args$noise_method <- "discrete"
    args$c <- 0
  }
  do.call(sydni, args)$data[, c("group", "x", "y")]
}

run_repeat <- function(d, t0, m, seed, keep_estimates = TRUE) {
  warned <- character()
  started <- Sys.time()
  syn <- tryCatch(
    withCallingHandlers(construct(d, m, seed), warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) e)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  failed <- inherits(syn, "error")
  res <- list(status = if (failed) "error" else "ok",
    message = if (failed) conditionMessage(syn) else paste(unique(warned), collapse = " | "),
    n_returned = if (failed) NA_integer_ else nrow(syn),
    n_returned_b = if (failed) NA_integer_ else sum(syn$group == "B"),
    elapsed = elapsed, change = NA_real_, slope_generated_b = NA_real_)
  if (!failed && keep_estimates) {
    pooled <- rbind(d[, c("group", "x", "y")], syn)
    res$change <- slope(pooled$x, pooled$y) - t0
    b <- syn[syn$group == "B", , drop = FALSE]
    res$slope_generated_b <- if (nrow(b) > 2L) slope(b$x, b$y) else NA_real_
  }
  res
}

run_unit <- function(cell, rep, repeats, offset, design, methods, keep_estimates = TRUE) {
  d <- make_data(cell, rep, offset, design)
  t0 <- slope(d$x, d$y)
  t_w <- slope(d$x, d$y, balance_weights(d$group))
  dataset <- data.frame(cell = cell$cell, base_cell = cell$base_cell, rep = rep, t0 = t0,
    change_weighted = t_w - t0, sd_x = sd(d$x), sd_y = sd(d$y),
    slope_a = slope(d$x[d$group == "A"], d$y[d$group == "A"]),
    slope_b = slope(d$x[d$group == "B"], d$y[d$group == "B"]))
  rows <- vector("list", nrow(methods) * repeats)
  j <- 0L
  for (mi in seq_len(nrow(methods))) {
    m <- methods[mi, ]
    for (g in seq_len(repeats)) {
      j <- j + 1L
      s <- stream("generator", cell$base_cell, rep, m$index, g, offset)
      r <- run_repeat(d, t0, m, s, keep_estimates)
      rows[[j]] <- data.frame(cell = cell$cell, base_cell = cell$base_cell, rep = rep, method = m$label,
        repeat_id = g, seed = s, status = r$status, n_returned = r$n_returned,
        n_returned_b = r$n_returned_b, change = r$change, slope_generated_b = r$slope_generated_b,
        elapsed = r$elapsed, message = r$message, stringsAsFactors = FALSE)
    }
  }
  list(dataset = dataset, repeats = do.call(rbind, rows))
}

# ---- provenance ----------------------------------------------------------------------

script_path <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) normalizePath(sub("^--file=", "", a[[1]])) else NA_character_
}

write_manifest <- function(path, args, design, reps, repeats, n_units) {
  sp <- script_path()
  lines <- c(
    paste("mode:", args$mode),
    paste("started_utc:", format(Sys.time(), tz = "UTC")),
    paste("r_version:", R.version.string),
    paste("sydni_version:", as.character(packageVersion("sydni"))),
    paste("script_md5:", if (is.na(sp)) "unknown" else unname(tools::md5sum(sp))),
    paste("rng_kind:", paste(RNGkind(), collapse = " / ")),
    paste("shard:", args$shard, "of", args$shards),
    paste("reps:", reps),
    paste("repeats:", repeats),
    paste("units_in_shard:", n_units))
  writeLines(lines, path)
}

# ---- software checks: never print an outcome ---------------------------------------------

run_checks <- function(design, methods) {
  offset <- 1e9
  ok <- TRUE
  report <- function(name, pass) {
    cat(sprintf("%-4s %s\n", if (isTRUE(pass)) "PASS" else "FAIL", name))
    if (!isTRUE(pass)) ok <<- FALSE
  }
  cells <- build_conditions(design)
  report("54 cells and 18 base cells", nrow(cells) == 54L && max(cells$base_cell) == 18L)
  report("clean references match the closed form",
    isTRUE(all.equal(cells$delta_artifact, (0.5 - cells$share) * (cells$b_b - cells$b_a))))
  report("no contamination makes both references equal",
    isTRUE(all.equal(cells$delta_artifact[cells$severity == 0], cells$delta_genuine[cells$severity == 0])))
  report("group sizes match the design for every cell",
    all(cells$n_a + cells$n_b == cells$n_total) && all(cells$n_b > cells$k_out))

  ids <- unlist(lapply(1:18, function(b) lapply(1:3, function(r) c(stream("data", b, r),
    stream("outlier", b, r), unlist(lapply(2:9, function(m) stream("generator", b, r, m, 1:10)))))))
  report("no stream is reused across a sample of units", !anyDuplicated(ids))
  report("largest stream fits an R integer",
    !is.na(stream("generator", 18L, design$reps, 9L, design$repeats, offset)))

  cl <- cells[cells$base_cell == 8L & cells$severity == 0, ]
  c2 <- cells[cells$base_cell == 8L & cells$severity == 2, ]
  c3 <- cells[cells$base_cell == 8L & cells$severity == 3, ]
  d0 <- make_data(cl, 1L, offset, design)
  d2 <- make_data(c2, 1L, offset, design)
  d3 <- make_data(c3, 1L, offset, design)
  changed <- which(d0$x != d2$x | d0$y != d2$y)
  report("contamination changes only the first k group B rows",
    identical(changed, seq.int(cl$n_a + 1L, length.out = c2$k_out)))
  report("the two severities use the same draws",
    isTRUE(all.equal(d3$x[changed] - 3, d2$x[changed] - 2)))
  report("group sizes in generated data", sum(d0$group == "B") == cl$n_b)

  w <- balance_weights(d0$group)
  report("weighted slope equals lm with weights",
    isTRUE(all.equal(slope(d0$x, d0$y, w), unname(coef(lm(y ~ x, d0, weights = w))[["x"]]))))
  report("weights give the groups equal total weight",
    isTRUE(all.equal(sum(w[d0$group == "A"]), sum(w[d0$group == "B"]))))

  t0 <- proc.time()[["elapsed"]]
  u1 <- run_unit(cl, 1L, 2L, offset, design, methods)
  cat(sprintf("     two repeats of all 8 stochastic methods at n = %d took %.1f seconds\n",
    cl$n_total, proc.time()[["elapsed"]] - t0))
  u2 <- run_unit(cl, 1L, 2L, offset, design, methods)
  report("one row per method and repeat", nrow(u1$repeats) == nrow(methods) * 2L)
  report("every construction ran", all(u1$repeats$status == "ok"))
  need <- cl$n_a - cl$n_b
  report("every construction returned exactly the rows needed for parity",
    all(u1$repeats$n_returned == need) && all(u1$repeats$n_returned_b == need))
  report("rerunning a unit reproduces it exactly",
    isTRUE(all.equal(u1$repeats[, c("change", "slope_generated_b")], u2$repeats[, c("change", "slope_generated_b")])))
  by_method <- tapply(u1$repeats$elapsed, u1$repeats$method, mean)
  cat("     mean seconds per repeat:", paste(names(by_method), round(by_method, 2), collapse = ", "), "\n")

  d <- d0
  none <- methods[methods$label == "M2", ]
  plain <- construct(d, none, 77L)
  rows <- paste(round(d$x, 10), round(d$y, 10))
  report("random oversampling returns only observed group B rows",
    all(paste(round(plain$x, 10), round(plain$y, 10)) %in% rows[d$group == "B"]))
  noisy <- construct(d, methods[methods$label == "M3", ], 77L)
  report("noise-perturbed oversampling changes the rows",
    mean(paste(round(noisy$x, 10), round(noisy$y, 10)) %in% rows) < 0.5)

  guard <- tryCatch({ guard_full(list(mode = "full", confirm = FALSE)); FALSE }, error = function(e) TRUE)
  report("full mode refuses without the authorization flag", guard)
  if (!ok) quit(save = "no", status = 1)
  cat("All checks passed.\n")
}

guard_full <- function(args) {
  if (args$mode == "full" && !isTRUE(args$confirm)) {
    stop("Full mode runs the declared study once. Pass --confirm-authorized-run only after the owner has authorized it.")
  }
  invisible(TRUE)
}

# ---- main --------------------------------------------------------------------------------

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (args$mode == "check") return(run_checks(DESIGN, METHODS))
  guard_full(args)

  cells <- build_conditions(DESIGN)
  reps <- if (is.na(args$reps)) DESIGN$reps else args$reps
  repeats <- if (is.na(args$repeats)) DESIGN$repeats else args$repeats
  units <- merge(data.frame(cell = cells$cell), data.frame(rep = seq_len(reps)), all = TRUE)
  units <- units[order(units$cell, units$rep), ]
  units <- units[(seq_len(nrow(units)) - 1L) %% args$shards == args$shard - 1L, ]

  pilot <- args$mode == "pilot"
  offset <- if (pilot) 5e8 else 0
  run_dir <- file.path(args$out, args$mode)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- sprintf("shard-%03d-of-%03d", args$shard, args$shards)
  files <- file.path(run_dir, paste0(stem, c("-datasets.csv.gz", "-repeats.csv.gz", "-manifest.txt", "-timing.csv")))
  if (!pilot && any(file.exists(files[1:3]))) stop("This shard already has output. The study runs once; nothing was overwritten.")
  if (!pilot) utils::write.csv(cells, file.path(run_dir, "conditions.csv"), row.names = FALSE)
  write_manifest(files[[3]], args, DESIGN, reps, repeats, nrow(units))

  datasets <- vector("list", nrow(units))
  results <- vector("list", nrow(units))
  for (i in seq_len(nrow(units))) {
    u <- run_unit(cells[cells$cell == units$cell[i], ], units$rep[i], repeats, offset, DESIGN, METHODS,
      keep_estimates = !pilot)
    datasets[[i]] <- u$dataset
    results[[i]] <- u$repeats
  }
  all_rows <- do.call(rbind, results)
  if (pilot) {
    timing <- aggregate(elapsed ~ method + cell, all_rows, function(z) c(mean = mean(z), max = max(z), n = length(z)))
    timing <- do.call(data.frame, timing)
    utils::write.csv(timing, files[[4]], row.names = FALSE)
    cat("Pilot wrote timings only:", files[[4]], "\n")
    return(invisible(NULL))
  }
  gz <- function(x, p) {
    con <- gzfile(p, "w")
    on.exit(close(con))
    utils::write.csv(x, con, row.names = FALSE)
  }
  gz(do.call(rbind, datasets), files[[1]])
  gz(all_rows, files[[2]])
  cat("Wrote", nrow(all_rows), "repeat rows for", nrow(units), "units to", run_dir, "\n")
}

main()
