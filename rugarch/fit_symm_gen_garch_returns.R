#!/usr/bin/env Rscript

# Fit the symmetric-news-impact model set used by xfit_symm_gen_garch_returns.f90.
# rugarch is used where it has a matching model. Small local likelihoods are used
# for EWMA, HARCH, RiskMetrics 2006, and the Fortran MIDAS-hyperbolic row.

suppressPackageStartupMessages(library(rugarch))

`%||%` <- function(x, y) if (is.null(x)) y else x

trading_days <- 252.0
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "."
root_dir <- normalizePath(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), ".."), mustWork = FALSE)
default_prices <- file.path(root_dir, "spy_efa_eem_tlt_lqd.csv")
default_output <- file.path(root_dir, "rugarch", "symm_gen_garch_results.csv")
model_names <- c(
  "SYMM_GARCH", "SYMM_GARCH_2_1", "SYMM_GARCH_1_2", "SYMM_GARCH_2_2",
  "FIGARCH", "CSGARCH", "HARCH", "RM2006", "MIDASHYP", "EWMA"
)

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    prices = default_prices,
    output = default_output,
    assets = NULL,
    models = model_names,
    scale = 100.0,
    maxit = 1000
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--prices", "--output", "--scale", "--maxit")) {
      if (i == length(args)) stop("Missing value after ", key)
      val <- args[[i + 1]]
      if (key == "--prices") out$prices <- val
      if (key == "--output") out$output <- val
      if (key == "--scale") out$scale <- as.numeric(val)
      if (key == "--maxit") out$maxit <- as.integer(val)
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

nll_from_variance <- function(y, h) {
  h <- pmax(h, 1.0e-12)
  mean(0.5 * log(2.0 * pi) + 0.5 * (log(h) + y^2 / h))
}

sigmoid <- function(x) 1.0 / (1.0 + exp(-x))
logit <- function(x) {
  x <- min(max(x, 1.0e-8), 1.0 - 1.0e-8)
  log(x / (1.0 - x))
}

simplex_transform <- function(p) {
  e <- exp(p)
  e / (1.0 + sum(e))
}

simplex_inv <- function(x) {
  slack <- max(1.0 - sum(x), 1.0e-8)
  log(pmax(x, 1.0e-12) / slack)
}

garch_variance <- function(y, omega, alpha, beta) {
  n <- length(y)
  p <- length(alpha)
  q <- length(beta)
  h <- numeric(n)
  backcast <- max(mean(y^2), 1.0e-12)
  for (t in seq_len(n)) {
    ht <- omega
    for (i in seq_len(p)) {
      idx <- t - i
      ht <- ht + alpha[[i]] * if (idx >= 1) y[[idx]]^2 else backcast
    }
    for (j in seq_len(q)) {
      idx <- t - j
      ht <- ht + beta[[j]] * if (idx >= 1) h[[idx]] else backcast
    }
    h[[t]] <- max(ht, 1.0e-12)
  }
  h
}

ewma_variance <- function(y, lambda) {
  n <- length(y)
  h <- numeric(n)
  ht <- max(mean(y^2), 1.0e-12)
  for (t in seq_len(n)) {
    h[[t]] <- ht
    ht <- lambda * ht + (1.0 - lambda) * y[[t]]^2
    ht <- max(ht, 1.0e-12)
  }
  h
}

fit_ewma_local <- function(y, maxit) {
  starts <- c(0.90, 0.94, 0.97, 0.99)
  best <- NULL
  for (s in starts) {
    opt <- optim(logit(s), function(p) nll_from_variance(y, ewma_variance(y, sigmoid(p[[1]]))),
                 method = "BFGS", control = list(maxit = maxit))
    if (is.null(best) || opt$value < best$value) best <- opt
  }
  lambda <- sigmoid(best$par[[1]])
  h <- ewma_variance(y, lambda)
  list(
    omega = 0.0, alpha = 1.0 - lambda, gamma = 0.0, beta = lambda, theta = 0.0,
    twist = 0.0, persist = 1.0, nparam = 1, niter = unname(best$counts[["function"]]),
    conv = best$convergence == 0, f = best$value, variance = h
  )
}

harch_lag_matrix <- function(y) {
  n <- length(y)
  backcast <- max(mean(y^2), 1.0e-12)
  x <- matrix(0.0, nrow = n, ncol = 3)
  lags <- c(1L, 5L, 22L)
  for (t in seq_len(n)) {
    for (k in seq_along(lags)) {
      lag <- lags[[k]]
      total <- 0.0
      for (j in seq_len(lag)) {
        idx <- t - j
        total <- total + if (idx >= 1) y[[idx]]^2 else backcast
      }
      x[t, k] <- total / lag
    }
  }
  x
}

harch_variance_from_matrix <- function(x, omega, alpha1, alpha5, alpha22) {
  as.numeric(pmax(omega + x %*% c(alpha1, alpha5, alpha22), 1.0e-12))
}

fit_harch_local <- function(y, maxit) {
  starts <- rbind(
    c(1.0e-6, 0.10, 0.40, 0.40),
    c(1.0e-6, 0.00, 0.45, 0.40),
    c(1.0e-6, 0.10, 0.20, 0.60),
    c(5.0e-6, 0.05, 0.30, 0.50)
  )
  x <- harch_lag_matrix(y)
  best <- NULL
  for (i in seq_len(nrow(starts))) {
    omega0 <- max((1.0 - sum(starts[i, 2:4])) * mean(y^2), 1.0e-12)
    p0 <- c(log(omega0), simplex_inv(starts[i, 2:4]))
    opt <- optim(p0, function(p) {
      a <- simplex_transform(p[2:4])
      nll_from_variance(y, harch_variance_from_matrix(x, exp(p[[1]]), a[[1]], a[[2]], a[[3]]))
    }, method = "BFGS", control = list(maxit = maxit))
    if (is.null(best) || opt$value < best$value) best <- opt
  }
  a <- simplex_transform(best$par[2:4])
  h <- harch_variance_from_matrix(x, exp(best$par[[1]]), a[[1]], a[[2]], a[[3]])
  list(
    omega = exp(best$par[[1]]), alpha = a[[1]], gamma = a[[2]], beta = a[[3]],
    theta = 0.0, twist = 0.0, persist = sum(a), nparam = 4,
    niter = unname(best$counts[["function"]]), conv = best$convergence == 0,
    f = best$value, variance = h, h_unc = exp(best$par[[1]]) / max(1.0 - sum(a), 1.0e-8)
  )
}

rm2006_variance <- function(y) {
  kmax <- 14
  tau0 <- 1560.0
  tau1 <- 4.0
  rho <- sqrt(2.0)
  tau <- tau1 * rho^(0:(kmax - 1))
  weights <- 1.0 - log(tau) / log(tau0)
  weights <- weights / sum(weights)
  mu <- exp(-1.0 / tau)
  n <- length(y)
  backcast <- numeric(kmax)
  for (k in seq_len(kmax)) {
    endpoint <- floor(log(0.01) / log(mu[[k]]))
    endpoint <- as.integer(max(min(endpoint, n), k - 1, 1))
    w <- mu[[k]]^(0:(endpoint - 1))
    backcast[[k]] <- sum(w * y[seq_len(endpoint)]^2) / sum(w)
  }
  comp <- backcast
  h <- numeric(n)
  for (t in seq_len(n)) {
    h[[t]] <- sum(weights * comp)
    comp <- mu * comp + (1.0 - mu) * y[[t]]^2
  }
  pmax(h, 1.0e-12)
}

fit_rm2006_local <- function(y) {
  h <- rm2006_variance(y)
  list(
    omega = 0.0, alpha = 0.0, gamma = 0.0, beta = 0.0, theta = 0.0,
    twist = 0.0, persist = 1.0, nparam = 0, niter = 0, conv = TRUE,
    f = nll_from_variance(y, h), variance = h, h_unc = mean(h)
  )
}

midas_weights <- function(theta) {
  raw <- numeric(22)
  raw[[1]] <- theta
  for (i in 2:22) raw[[i]] <- raw[[i - 1]] * ((i - 1) + theta) / i
  raw / sum(raw)
}

midas_backcast <- function(y) {
  tau <- min(75, length(y))
  w <- 0.94^(0:(tau - 1))
  max(sum(w * y[seq_len(tau)]^2) / sum(w), 1.0e-12)
}

midas_lag_matrix <- function(y) {
  n <- length(y)
  backcast <- midas_backcast(y)
  x <- matrix(0.0, nrow = n, ncol = 22)
  for (t in seq_len(n)) {
    for (i in seq_len(22)) {
      idx <- t - i
      x[t, i] <- if (idx >= 1) y[[idx]]^2 else backcast
    }
  }
  x
}

midas_variance_from_matrix <- function(x, omega, alpha, theta) {
  weights <- midas_weights(theta)
  as.numeric(pmax(omega + alpha * (x %*% weights), 1.0e-12))
}

fit_midas_local <- function(y, maxit) {
  starts_alpha <- c(0.80, 0.90, 0.95, 0.98)
  starts_theta <- c(0.10, 0.50, 0.80, 0.90)
  x <- midas_lag_matrix(y)
  best <- NULL
  for (i in seq_along(starts_alpha)) {
    omega0 <- max((1.0 - min(starts_alpha[[i]], 0.99)) * mean(y^2), 1.0e-12)
    p0 <- c(log(omega0), logit(starts_alpha[[i]]), logit(starts_theta[[i]]))
    opt <- optim(p0, function(p) {
      alpha <- sigmoid(p[[2]])
      theta <- sigmoid(p[[3]])
      nll_from_variance(y, midas_variance_from_matrix(x, exp(p[[1]]), alpha, theta))
    }, method = "BFGS", control = list(maxit = maxit))
    if (is.null(best) || opt$value < best$value) best <- opt
  }
  alpha <- sigmoid(best$par[[2]])
  theta <- sigmoid(best$par[[3]])
  omega <- exp(best$par[[1]])
  h <- midas_variance_from_matrix(x, omega, alpha, theta)
  list(
    omega = omega, alpha = alpha, gamma = 0.0, beta = 0.0, theta = theta,
    twist = 0.0, persist = alpha, nparam = 3,
    niter = unname(best$counts[["function"]]), conv = best$convergence == 0,
    f = best$value, variance = h, h_unc = omega / max(1.0 - alpha, 1.0e-8)
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

rugarch_model_spec <- function(model) {
  switch(model,
    SYMM_GARCH = list(kind = "sGARCH", order = c(1, 1)),
    SYMM_GARCH_2_1 = list(kind = "sGARCH", order = c(2, 1)),
    SYMM_GARCH_1_2 = list(kind = "sGARCH", order = c(1, 2)),
    SYMM_GARCH_2_2 = list(kind = "sGARCH", order = c(2, 2)),
    FIGARCH = list(kind = "fiGARCH", order = c(1, 1)),
    CSGARCH = list(kind = "csGARCH", order = c(1, 1)),
    NULL
  )
}

fit_rugarch_one <- function(y, model, scale, maxit) {
  spec_info <- rugarch_model_spec(model)
  ys <- y * scale
  spec <- ugarchspec(
    variance.model = list(model = spec_info$kind, garchOrder = spec_info$order),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"
  )
  fit <- ugarchfit(
    spec, ys, solver = "hybrid",
    solver.control = list(trace = 0),
    fit.control = list(scale = 0),
    numderiv.control = list(),
    out.sample = 0
  )
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
  }
  h <- as.numeric(sigma(fit))^2 / (scale * scale)
  h_unc <- mean(h)
  if (model == "SYMM_GARCH") h_unc <- get("omega") / (scale * scale) / max(1.0 - persist, 1.0e-8)
  list(
    omega = get("omega") / (scale * scale), alpha = alpha, gamma = gamma,
    beta = beta, theta = theta, twist = 0.0, persist = persist,
    nparam = length(cf), niter = if (!is.null(fit@fit$solver$sol$iterations)) fit@fit$solver$sol$iterations else NA_integer_,
    conv = convergence(fit) == 0,
    f = -(likelihood(fit) + length(y) * log(scale)) / length(y),
    variance = pmax(h, 1.0e-12), h_unc = h_unc
  )
}

fit_one <- function(y, model, scale, maxit) {
  if (!is.null(rugarch_model_spec(model))) return(fit_rugarch_one(y, model, scale, maxit))
  if (model == "EWMA") return(fit_ewma_local(y, maxit))
  if (model == "HARCH") return(fit_harch_local(y, maxit))
  if (model == "RM2006") return(fit_rm2006_local(y))
  if (model == "MIDASHYP") return(fit_midas_local(y, maxit))
  stop("No fit implementation for model ", model)
}

format_date <- function(x) {
  x <- as.character(x)
  if (grepl("^\\d{8}$", x)) paste0(substr(x, 1, 4), "-", substr(x, 5, 6), "-", substr(x, 7, 8)) else x
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

make_row <- function(asset, model, fit, nobs) {
  logl <- -nobs * fit$f
  z <- demeaned_z <- NULL
  z <- current_y / sqrt(pmax(fit$variance, 1.0e-12))
  mom <- normal_moments(z)
  list(
    model = model, asset = asset, omega = fit$omega, alpha = fit$alpha,
    gamma = fit$gamma, beta = fit$beta, theta = fit$theta, twist = fit$twist,
    persist = fit$persist,
    vol_ann = sqrt(trading_days * (fit$h_unc %||% mean(fit$variance))) * 100.0,
    logl = logl, aic = 2.0 * fit$nparam - 2.0 * logl,
    bic = log(nobs) * fit$nparam - 2.0 * logl,
    nparam = fit$nparam, niter = fit$niter, conv = fit$conv,
    skew = mom[["skew"]], ekurt = mom[["ekurt"]], sec = fit$sec %||% NA_real_,
    aic_rank = NA_integer_, bic_rank = NA_integer_
  )
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
  summary <- lapply(models, function(m) {
    r <- rows[vapply(rows, function(x) x$model == m, logical(1))]
    if (!length(r)) {
      return(data.frame(model = m, nparam = 0, aic_wins = 0, bic_wins = 0,
                        aic_avg_rank = 0, bic_avg_rank = 0, aic_symbols = ""))
    }
    aic_winners <- vapply(r, function(x) x$aic_rank == 1, logical(1))
    bic_winners <- vapply(r, function(x) x$bic_rank == 1, logical(1))
    data.frame(
      model = m,
      nparam = r[[length(r)]]$nparam,
      aic_wins = sum(aic_winners),
      bic_wins = sum(bic_winners),
      aic_avg_rank = mean(vapply(r, function(x) x$aic_rank, numeric(1))),
      bic_avg_rank = mean(vapply(r, function(x) x$bic_rank, numeric(1))),
      aic_symbols = paste(vapply(r[aic_winners], function(x) x$asset, character(1)), collapse = " ")
    )
  })
  df <- do.call(rbind, summary)
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
    data.frame(model = m,
               nparam = if (length(r)) r[[length(r)]]$nparam else 0,
               sec = sum(vapply(r, function(x) x$sec, numeric(1)), na.rm = TRUE))
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
  df <- do.call(rbind, lapply(rows, as.data.frame))
  write.csv(df, path, row.names = FALSE)
}

args <- parse_args()
unknown <- setdiff(args$models, model_names)
if (length(unknown)) stop("Unknown model(s): ", paste(unknown, collapse = ", "))

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
cat("Package: rugarch plus local likelihoods for EWMA/HARCH/RM2006/MIDASHYP\n")
cat(sprintf("rugarch fits use returns multiplied by %.6g; printed omega/logL/AIC/BIC are on raw-return scale.\n",
            args$scale))

rows <- list()
for (asset in assets) {
  current_y <- demeaned_log_returns(data$prices[[asset]])
  for (model in args$models) {
    t_fit <- proc.time()[["elapsed"]]
    fit <- fit_one(current_y, model, args$scale, args$maxit)
    fit$sec <- proc.time()[["elapsed"]] - t_fit
    rows[[length(rows) + 1L]] <- make_row(asset, model, fit, length(current_y))
  }
  rows <- rank_rows(rows, asset)
}

print_main_table(rows)
print_selection_counts(rows, args$models)
print_fit_times(rows, args$models)
write_results(rows, args$output)

cat(sprintf("\nWrote CSV results: %s\n", args$output))
cat(sprintf("Elapsed wall time: %10.3f seconds\n", proc.time()[["elapsed"]] - t0))
