"""Fit Python arch models that overlap with the Fortran GARCH code.

The script reads the same daily price CSV files used by the Fortran programs,
computes demeaned log returns, and fits common Normal-innovation volatility
models using the installed ``arch`` package.

Default scale is 100.0 since ``arch`` optimizes more reliably on percent returns.
Omega and log-likelihood are converted back to the raw-return scale in the
printed table.
"""

from __future__ import annotations

import argparse
import csv
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from arch.univariate import ZeroMean, arch_model
from arch.univariate.distribution import Normal
from arch.univariate.volatility import FIGARCH, MIDASHyperbolic, RiskMetrics2006


TRADING_DAYS = 252.0
ROOT = Path(__file__).resolve().parents[1]

MODEL_SPECS = {
    "ARCH1": {"vol": "ARCH", "p": 1, "o": 0, "q": 0, "power": 2.0},
    "SYMM_GARCH": {"vol": "GARCH", "p": 1, "o": 0, "q": 1, "power": 2.0},
    "FIGARCH": {"vol": "FIGARCH", "p": 1, "q": 1, "power": 2.0},
    "GJR_GARCH": {"vol": "GARCH", "p": 1, "o": 1, "q": 1, "power": 2.0},
    "EGARCH": {"vol": "EGARCH", "p": 1, "o": 1, "q": 1},
    "APARCH": {"vol": "APARCH", "p": 1, "o": 1, "q": 1},
    "HARCH": {"vol": "HARCH", "p": [1, 5, 22]},
    "RM2006": {"vol": "RiskMetrics2006"},
    "MIDASHYP": {"vol": "MIDASHyperbolic"},
    "MIDASHYP_ASYM": {"vol": "MIDASHyperbolic", "asym": True},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prices", default=str(ROOT / "spy_efa_eem_tlt_lqd.csv"))
    parser.add_argument("--assets", nargs="*", default=None)
    parser.add_argument("--models", nargs="*", default=list(MODEL_SPECS))
    parser.add_argument("--scale", type=float, default=100.0)
    parser.add_argument("--output", default=str(ROOT / "arch" / "arch_common_results.csv"))
    parser.add_argument("--maxiter", type=int, default=1000)
    return parser.parse_args()


def read_prices(path: Path) -> pd.DataFrame:
    prices = pd.read_csv(path, index_col=0)
    prices.index.name = "date"
    return prices.apply(pd.to_numeric, errors="coerce")


def demeaned_log_returns(prices: pd.Series) -> np.ndarray:
    values = prices.to_numpy(dtype=float)
    ret = np.diff(np.log(values))
    return ret - ret.mean()


def normal_moments(z: np.ndarray) -> tuple[float, float]:
    zc = z - z.mean()
    var = float(np.mean(zc * zc))
    if var <= 0.0:
        return 0.0, 0.0
    skew = float(np.mean(zc**3) / var**1.5)
    ekurt = float(np.mean(zc**4) / var**2 - 3.0)
    return skew, ekurt


def param_value(params: pd.Series, *names: str) -> float:
    for name in names:
        if name in params:
            return float(params[name])
    return 0.0


def figarch_persistence(phi: float, d: float, beta: float, truncation: int = 1000) -> float:
    lam = np.empty(truncation)
    lam[0] = phi - beta + d
    delta_prev = d
    for i in range(1, truncation):
        delta_cur = (i - d) / (i + 1.0) * delta_prev
        lam[i] = beta * lam[i - 1] + delta_cur - phi * delta_prev
        delta_prev = delta_cur
    return float(np.sum(lam))


def fit_one(asset: str, y_raw: np.ndarray, model_name: str, scale: float, maxiter: int) -> dict[str, object]:
    spec = MODEL_SPECS[model_name]
    y = y_raw * scale
    started = time.perf_counter()
    if model_name == "RM2006":
        model = ZeroMean(y, volatility=RiskMetrics2006(), distribution=Normal(), rescale=False)
    elif model_name == "FIGARCH":
        model = ZeroMean(y, volatility=FIGARCH(p=1, q=1, power=2.0), distribution=Normal(), rescale=False)
    elif model_name == "MIDASHYP":
        model = ZeroMean(y, volatility=MIDASHyperbolic(), distribution=Normal(), rescale=False)
    elif model_name == "MIDASHYP_ASYM":
        model = ZeroMean(y, volatility=MIDASHyperbolic(asym=True), distribution=Normal(), rescale=False)
    else:
        model = arch_model(y, mean="Zero", dist="normal", rescale=False, **spec)
    result = model.fit(disp="off", options={"maxiter": maxiter})
    elapsed = time.perf_counter() - started

    params = result.params
    omega_scaled = param_value(params, "omega")
    alpha = param_value(params, "alpha[1]", "alpha")
    gamma = param_value(params, "gamma[1]", "gamma")
    beta = param_value(params, "beta[1]")
    delta = param_value(params, "delta")
    if model_name == "FIGARCH":
        alpha = param_value(params, "phi")
        delta = param_value(params, "d")
        beta = param_value(params, "beta")

    if model_name == "EGARCH":
        omega_raw = omega_scaled - (1.0 - beta) * 2.0 * math.log(scale)
        persist = beta
        h_unc_raw = math.exp(omega_scaled / max(1.0 - beta, 1.0e-8)) / (scale * scale)
    elif model_name == "MIDASHYP":
        delta = param_value(params, "theta")
        omega_raw = omega_scaled / (scale * scale)
        persist = alpha
        h_unc_raw = omega_raw / max(1.0 - persist, 1.0e-8)
    elif model_name == "MIDASHYP_ASYM":
        delta = param_value(params, "theta")
        omega_raw = omega_scaled / (scale * scale)
        persist = alpha + 0.5 * gamma
        h_unc_raw = omega_raw / max(1.0 - persist, 1.0e-8)
    elif model_name == "RM2006":
        omega_raw = 0.0
        persist = 1.0
        cond_vol_raw = np.asarray(result.conditional_volatility, dtype=float) / scale
        h_unc_raw = float(np.mean(cond_vol_raw * cond_vol_raw))
    elif model_name == "FIGARCH":
        delta = param_value(params, "d")
        omega_raw = omega_scaled / (scale * scale)
        cond_vol_raw = np.asarray(result.conditional_volatility, dtype=float) / scale
        h_unc_raw = float(np.mean(cond_vol_raw * cond_vol_raw))
        persist = figarch_persistence(alpha, delta, beta)
    elif model_name == "APARCH":
        omega_raw = omega_scaled / (scale**delta)
        persist = alpha + beta
        cond_vol_raw = np.asarray(result.conditional_volatility, dtype=float) / scale
        h_unc_raw = float(np.mean(cond_vol_raw * cond_vol_raw))
    elif model_name == "HARCH":
        gamma = param_value(params, "alpha[5]")
        beta = param_value(params, "alpha[22]")
        omega_raw = omega_scaled / (scale * scale)
        persist = alpha + gamma + beta
        h_unc_raw = omega_raw / max(1.0 - persist, 1.0e-8)
    elif model_name == "GJR_GARCH":
        omega_raw = omega_scaled / (scale * scale)
        persist = alpha + 0.5 * gamma + beta
        h_unc_raw = omega_raw / max(1.0 - persist, 1.0e-8)
    else:
        omega_raw = omega_scaled / (scale * scale)
        persist = alpha + beta
        h_unc_raw = omega_raw / max(1.0 - persist, 1.0e-8)

    logl_raw = float(result.loglikelihood) + len(y_raw) * math.log(scale)
    nparam = int(result.num_params)
    aic_raw = 2.0 * nparam - 2.0 * logl_raw
    bic_raw = math.log(len(y_raw)) * nparam - 2.0 * logl_raw

    cond_vol_raw = np.asarray(result.conditional_volatility, dtype=float) / scale
    z = y_raw / np.maximum(cond_vol_raw, 1.0e-12)
    skew, ekurt = normal_moments(z)

    return {
        "asset": asset,
        "model": model_name,
        "omega": omega_raw,
        "alpha": alpha,
        "gamma": gamma,
        "beta": beta,
        "delta": delta,
        "persist": persist,
        "vol_ann_pct": math.sqrt(max(h_unc_raw, 0.0) * TRADING_DAYS) * 100.0,
        "logL": logl_raw,
        "AIC": aic_raw,
        "BIC": bic_raw,
        "nparam": nparam,
        "iter": int(getattr(result.optimization_result, "nit", -1)),
        "conv": result.convergence_flag == 0,
        "skew": skew,
        "ekurt": ekurt,
        "sec": elapsed,
    }


def print_table(rows: list[dict[str, object]], prices_file: Path, nobs: int, nassets: int, start: str, end: str) -> None:
    print(f"Prices file: {prices_file}")
    print(f"Using {nobs} demeaned log returns for {nassets} assets from {start} to {end}")
    print("Package: arch")
    print("Input returns are multiplied by --scale for fitting; printed logL and omega are converted to raw-return scale.")
    print(
        "Model            Asset        omega   alpha   gamma    beta   delta  persist  "
        "vol_ann%        logL         AIC         BIC #param iter conv    skew   ekurt      sec"
    )
    print("-" * 174)
    for row in rows:
        print(
            f"{row['model']:>16} {row['asset']:>9}"
            f"{row['omega']:12.3E}{row['alpha']:8.4f}{row['gamma']:8.4f}{row['beta']:8.4f}"
            f"{row['delta']:8.4f}"
            f"{row['persist']:9.4f}{row['vol_ann_pct']:10.2f}"
            f"{row['logL']:12.2f}{row['AIC']:12.2f}{row['BIC']:12.2f}"
            f"{row['nparam']:7d}{row['iter']:5d} {str(row['conv'])[0]}"
            f"{row['skew']:9.3f}{row['ekurt']:8.3f}{row['sec']:9.3f}"
        )


def write_csv(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    t_total0 = time.perf_counter()
    args = parse_args()
    prices_file = Path(args.prices)
    t_read0 = time.perf_counter()
    prices = read_prices(prices_file)
    assets = list(args.assets) if args.assets else list(prices.columns)
    models = [m.upper() for m in args.models]
    unknown = sorted(set(models) - set(MODEL_SPECS))
    if unknown:
        raise SystemExit(f"Unknown model(s): {', '.join(unknown)}")

    rows: list[dict[str, object]] = []
    nobs = prices.shape[0] - 1
    start = str(prices.index[1])
    end = str(prices.index[-1])
    t_read1 = time.perf_counter()

    t_fit0 = time.perf_counter()
    for asset in assets:
        y = demeaned_log_returns(prices[asset])
        for model_name in models:
            rows.append(fit_one(asset, y, model_name, args.scale, args.maxiter))
    t_fit1 = time.perf_counter()

    t_print0 = time.perf_counter()
    print_table(rows, prices_file, nobs, len(assets), start, end)
    t_print1 = time.perf_counter()

    t_write0 = time.perf_counter()
    write_csv(rows, Path(args.output))
    t_write1 = time.perf_counter()
    print(f"\nWrote {args.output}")
    t_total1 = time.perf_counter()
    print("\nTiming:")
    print(f"  read/prep:       {t_read1 - t_read0:10.3f} seconds")
    print(f"  fit models:      {t_fit1 - t_fit0:10.3f} seconds")
    print(f"  print output:    {t_print1 - t_print0:10.3f} seconds")
    print(f"  write CSV:       {t_write1 - t_write0:10.3f} seconds")
    print(f"  elapsed wall:    {t_total1 - t_total0:10.3f} seconds")


if __name__ == "__main__":
    main()
