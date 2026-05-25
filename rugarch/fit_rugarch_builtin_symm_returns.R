#!/usr/bin/env Rscript

# Fit the symmetric-news-impact models from xfit_symm_gen_garch_returns.f90
# that are directly built into rugarch.

suppressPackageStartupMessages(library(rugarch))

`%||%` <- function(x, y) if (is.null(x)) y else x

trading_days <- 252.0
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "."
root_dir <- normalizePath(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), ".."), mustWork = FALSE)
default_prices <- file.path(root_dir, "spy_efa_eem_tlt_lqd.csv")
default_output <- file.path(root_dir, "rugarch", "rugarch_builtin_symm_results.csv")

model_names <- c(
  "SYMM_GARCH", "SYMM_GARCH_2_1", "SYMM_GARCH_1_2", "SYMM_GARCH_2_2",
  "FIGARCH", "CSGARCH", "TGARCH", "AVGARCH"
)

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    prices = default_prices,
    output = default_output,
    assets = NULL,
    models = model_names,
    scale = 100.0
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--prices", "--output", "--scale")) {
      if (i == length(args)) stop("Missing value after ", key)
      val <- args[[i + 1]]
      if (key == "--prices") out$prices <- val
      if (key == "--output") out$output <- val
      if (key == "--scale") out$scale <- as.numeric(val)
      i <- i + 2
    } else if (key %in% c("--assets", "--models")) {
      vals <- character()
      i <- i + 1
      while (i <= length(args) && !startsWith(args[[i]], "--")) {
        vals <- c(vals, args[[i]])
        i <- i + 1
      }
      if (key == "--assets") out$assets <- vals
      if (key == "--models") out$models <- toupper(vals)
    } else {
      stop("Unknown argument: ", key)
    }
  }
  out
}

read_prices <- function(path) {
  x <- read.csv(path, check.names = FALSE)
  dates <- as.character(x[[1]])
  prices <- as.data.frame(lapply(x[-1], as.numeric), check.names = FALSE)
  names(prices) <- names(x)[-1]
  list(dates = dates, prices = prices)
}

demeaned_log_returns <- function(x) {
  r <- diff(log(as.numeric(x)))
  r - mean(r)
}

normal_moments <- function(z) {
  z <- z - mean(z)
  v <- mean(z^2)
  if (!is.finite(v) || v <= 0.0) return(c(skew = 0.0, ekurt = 0.0))
  c(skew = mean(z^3) / v^1.5, ekurt = mean(z^4) / v^2 - 3.0)
}

format_date <- function(x) {
  x <- as.character(x)
  if (grepl("^\\d{8}$", x)) paste0(substr(x, 1, 4), "-", substr(x, 5, 6), "-", substr(x, 7, 8)) else x
}

rugarch_model_spec <- function(model) {
  switch(model,
    SYMM_GARCH = list(kind = "sGARCH", order = c(1, 1)),
    SYMM_GARCH_2_1 = list(kind = "sGARCH", order = c(2, 1)),
    SYMM_GARCH_1_2 = list(kind = "sGARCH", order = c(1, 2)),
    SYMM_GARCH_2_2 = list(kind = "sGARCH", order = c(2, 2)),
    FIGARCH = list(kind = "fiGARCH", order = c(1, 1)),
    CSGARCH = list(kind = "csGARCH", order = c(1, 1)),
    TGARCH = list(kind = "fGARCH", submodel = "TGARCH", order = c(1, 1)),
    AVGARCH = list(kind = "fGARCH", submodel = "AVGARCH", order = c(1, 1)),
    NULL
  )
}

figarch_persist <- function(phi, d, beta, truncation = 1000) {
  lambda <- numeric(truncation)
  lambda[[1]] <- phi - beta + d
  delta_prev <- d
  if (truncation > 1) {
    for (i in 2:truncation) {
      delta_cur <- ((i - 1) - d) / i * delta_prev
      lambda[[i]] <- beta * lambda[[i - 1]] + delta_cur - phi * delta_prev
      delta_prev <- delta_cur
    }
  }
  sum(lambda)
}

fit_one <- function(y, model, scale) {
  spec_info <- rugarch_model_spec(model)
  if (is.null(spec_info)) stop("No rugarch built-in mapping for model ", model)

  ys <- y * scale
  spec <- ugarchspec(
    variance.model = if (spec_info$kind == "fGARCH") {
      list(model = spec_info$kind, submodel = spec_info$submodel, garchOrder = spec_info$order)
    } else {
      list(model = spec_info$kind, garchOrder = spec_info$order)
    },
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"
  )
  t0 <- proc.time()[["elapsed"]]
  fit <- ugarchfit(
    spec, ys, solver = "hybrid",
    solver.control = list(trace = 0),
    fit.control = list(scale = 0),
    out.sample = 0
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  cf <- coef(fit)
  get <- function(nm) if (nm %in% names(cf)) unname(cf[[nm]]) else 0.0
  alpha_names <- grep("^alpha", names(cf), value = TRUE)
  beta_names <- grep("^beta", names(cf), value = TRUE)
  alpha <- sum(cf[alpha_names])
  beta <- sum(cf[beta_names])
  gamma <- 0.0
  theta <- 0.0
  persist <- alpha + beta

  if (model == "FIGARCH") {
    alpha <- get("alpha1")
    beta <- get("beta1")
    theta <- get("delta")
    persist <- figarch_persist(alpha, theta, beta)
  } else if (model == "CSGARCH") {
    alpha <- get("alpha1")
    beta <- get("beta1")
    gamma <- get("eta21")
    persist <- get("eta11")
  } else if (model == "TGARCH") {
    alpha <- get("alpha1")
    beta <- get("beta1")
    gamma <- get("eta11")
    persist <- persistence(fit)
  } else if (model == "AVGARCH") {
    alpha <- get("alpha1")
    beta <- get("beta1")
    gamma <- get("eta11")
    theta <- get("eta21")
    persist <- persistence(fit)
  }

  h <- pmax(as.numeric(sigma(fit))^2 / (scale * scale), 1.0e-12)
  omega <- if (model %in% c("TGARCH", "AVGARCH")) get("omega") / scale else get("omega") / (scale * scale)
  h_unc <- mean(h)
  if (model == "SYMM_GARCH") h_unc <- omega / max(1.0 - persist, 1.0e-8)
  logl <- likelihood(fit) + length(y) * log(scale)
  z <- y / sqrt(h)
  mom <- normal_moments(z)

  nparam <- length(cf)
  list(
    model = model, omega = omega, alpha = alpha, gamma = gamma, beta = beta,
    theta = theta, twist = 0.0, persist = persist,
    vol_ann = sqrt(trading_days * h_unc) * 100.0,
    logl = logl, aic = 2.0 * nparam - 2.0 * logl,
    bic = log(length(y)) * nparam - 2.0 * logl,
    nparam = nparam, niter = NA_integer_, conv = convergence(fit) == 0,
    skew = mom[["skew"]], ekurt = mom[["ekurt"]], sec = elapsed
  )
}

rank_rows <- function(rows, asset) {
  idx <- which(vapply(rows, function(x) x$asset == asset, logical(1)))
  aic <- vapply(rows[idx], function(x) x$aic, numeric(1))
  bic <- vapply(rows[idx], function(x) x$bic, numeric(1))
  for (k in seq_along(idx)) {
    rows[[idx[[k]]]]$aic_rank <- 1L + sum(aic < aic[[k]])
    rows[[idx[[k]]]]$bic_rank <- 1L + sum(bic < bic[[k]])
  }
  rows
}

print_main_table <- function(rows) {
  cat("Model            Asset        omega   alpha   gamma    beta   theta   twist  persist  vol_ann%        logL         AIC         BIC #param iter conv    skew   ekurt AIC_rank BIC_rank\n")
  cat(paste0(strrep("-", 183), "\n"))
  for (row in rows) {
    cat(sprintf(
      "%16s %9s%12.3E%8.4f%8.4f%8.4f%8.4f%8.4f%9.4f%10.2f%12.2f%12.2f%12.2f%7d%5s %1s%8.3f%8.3f%9d%9d\n",
      row$model, row$asset, row$omega, row$alpha, row$gamma, row$beta,
      row$theta, row$twist, row$persist, row$vol_ann, row$logl, row$aic,
      row$bic, row$nparam, ifelse(is.na(row$niter), "NA", as.character(row$niter)),
      substr(as.character(row$conv), 1, 1), row$skew, row$ekurt,
      row$aic_rank, row$bic_rank
    ))
  }
}

print_selection_counts <- function(rows, models) {
  cat("\nModel selection counts:\n")
  cat("Model            #param  AIC_wins  BIC_wins AIC_avg_rank BIC_avg_rank  AIC_symbols\n")
  cat(paste0(strrep("-", 100), "\n"))
  df <- do.call(rbind, lapply(models, function(m) {
    r <- rows[vapply(rows, function(x) x$model == m, logical(1))]
    aic_winners <- vapply(r, function(x) x$aic_rank == 1, logical(1))
    bic_winners <- vapply(r, function(x) x$bic_rank == 1, logical(1))
    data.frame(
      model = m, nparam = r[[length(r)]]$nparam,
      aic_wins = sum(aic_winners), bic_wins = sum(bic_winners),
      aic_avg_rank = mean(vapply(r, function(x) x$aic_rank, numeric(1))),
      bic_avg_rank = mean(vapply(r, function(x) x$bic_rank, numeric(1))),
      aic_symbols = paste(vapply(r[aic_winners], function(x) x$asset, character(1)), collapse = " ")
    )
  }))
  df <- df[order(-df$aic_wins, -df$bic_wins, df$model), ]
  for (i in seq_len(nrow(df))) {
    cat(sprintf("%16s%8d%10d%10d%13.2f%13.2f  %s\n",
                df$model[[i]], df$nparam[[i]], df$aic_wins[[i]], df$bic_wins[[i]],
                df$aic_avg_rank[[i]], df$bic_avg_rank[[i]], df$aic_symbols[[i]]))
  }
}

print_fit_times <- function(rows, models) {
  cat("\nTotal fitting time by model:\n")
  cat("Model            #param  fit_seconds  cumul_time  cumul_frac\n")
  cat(paste0(strrep("-", 62), "\n"))
  df <- do.call(rbind, lapply(models, function(m) {
    r <- rows[vapply(rows, function(x) x$model == m, logical(1))]
    data.frame(model = m, nparam = r[[length(r)]]$nparam,
               sec = sum(vapply(r, function(x) x$sec, numeric(1))))
  }))
  df <- df[order(-df$sec), ]
  total <- sum(df$sec)
  cumul <- 0.0
  for (i in seq_len(nrow(df))) {
    cumul <- cumul + df$sec[[i]]
    frac <- if (total > 0.0) cumul / total else 0.0
    cat(sprintf("%16s%8d%13.4f%13.4f%13.4f\n", df$model[[i]], df$nparam[[i]], df$sec[[i]], cumul, frac))
  }
}

write_results <- function(rows, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(do.call(rbind, lapply(rows, as.data.frame)), path, row.names = FALSE)
}

args <- parse_args()
unknown <- setdiff(args$models, model_names)
if (length(unknown)) stop("Unknown rugarch built-in model(s): ", paste(unknown, collapse = ", "))

t0 <- proc.time()[["elapsed"]]
data <- read_prices(args$prices)
assets <- if (is.null(args$assets)) names(data$prices) else args$assets
missing_assets <- setdiff(assets, names(data$prices))
if (length(missing_assets)) stop("Missing asset columns: ", paste(missing_assets, collapse = ", "))
nobs <- nrow(data$prices) - 1
start_date <- format_date(data$dates[[2]])
end_date <- format_date(data$dates[[length(data$dates)]])

cat(sprintf("Prices file: %s\n", basename(args$prices)))
cat(sprintf("Using %d demeaned log returns for %d assets from %s to %s\n",
            nobs, length(assets), start_date, end_date))
cat("Package: rugarch built-in GARCH-family models\n")
cat(sprintf("rugarch fits use returns multiplied by %.6g; printed omega/logL/AIC/BIC are on raw-return scale.\n",
            args$scale))

rows <- list()
for (asset in assets) {
  y <- demeaned_log_returns(data$prices[[asset]])
  for (model in args$models) {
    fit <- fit_one(y, model, args$scale)
    fit$asset <- asset
    fit$aic_rank <- NA_integer_
    fit$bic_rank <- NA_integer_
    rows[[length(rows) + 1L]] <- fit
  }
  rows <- rank_rows(rows, asset)
}

print_main_table(rows)
print_selection_counts(rows, args$models)
print_fit_times(rows, args$models)
write_results(rows, args$output)

cat(sprintf("\nWrote CSV results: %s\n", args$output))
cat(sprintf("Elapsed wall time: %10.3f seconds\n", proc.time()[["elapsed"]] - t0))
