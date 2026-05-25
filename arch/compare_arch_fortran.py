"""Run Python arch and Fortran common-model fitters and compare results."""

from __future__ import annotations

import argparse
import subprocess
import time
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PY_SCRIPT = ROOT / "arch" / "fit_arch_common.py"
PY_CSV = ROOT / "arch" / "arch_common_results.csv"
FTN_CSV = ROOT / "arch" / "fortran_arch_common_results.csv"
COMPARE_CSV = ROOT / "arch" / "arch_fortran_comparison.csv"

PARAM_COLS = ["omega", "alpha", "gamma", "beta", "delta", "persist", "vol_ann_pct", "logL"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python-cmd", default="python")
    parser.add_argument("--make-cmd", default="make")
    parser.add_argument("--prices", default="spy_efa_eem_tlt_lqd.csv")
    parser.add_argument("--assets", nargs="*", default=None)
    parser.add_argument("--models", nargs="*", default=None)
    parser.add_argument("--scale", type=float, default=100.0)
    parser.add_argument("--rel-tol", type=float, default=0.03)
    parser.add_argument("--abs-tol", type=float, default=1.0e-5)
    parser.add_argument("--logl-tol", type=float, default=1.0)
    parser.add_argument("--comparison-output", type=Path, default=COMPARE_CSV)
    parser.add_argument("--keep-going", action="store_true")
    return parser.parse_args()


def run_command(cmd: list[str], label: str) -> tuple[float, str, str, int]:
    started = time.perf_counter()
    proc = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True)
    elapsed = time.perf_counter() - started
    print(f"{label} command: {' '.join(cmd)}")
    print(f"{label} wall time: {elapsed:.3f} seconds")
    if proc.stdout:
        print(f"\n--- {label} stdout ---")
        print(proc.stdout.rstrip())
    if proc.stderr:
        print(f"\n--- {label} stderr ---")
        print(proc.stderr.rstrip())
    print()
    return elapsed, proc.stdout, proc.stderr, proc.returncode


def python_command(args: argparse.Namespace) -> list[str]:
    cmd = [
        args.python_cmd,
        str(PY_SCRIPT),
        "--prices",
        args.prices,
        "--scale",
        str(args.scale),
    ]
    if args.assets:
        cmd += ["--assets", *args.assets]
    if args.models:
        cmd += ["--models", *args.models]
    return cmd


def read_results(path: Path, source: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["asset"] = df["asset"].astype(str).str.strip()
    df["model"] = df["model"].astype(str).str.strip()
    df["source"] = source
    return df


def significant_diff(value_py: float, value_ftn: float, col: str, args: argparse.Namespace) -> bool:
    diff = abs(value_py - value_ftn)
    if col == "logL":
        return diff > args.logl_tol
    scale = max(abs(value_py), abs(value_ftn), args.abs_tol)
    return diff > args.abs_tol and diff / scale > args.rel_tol


def compare(args: argparse.Namespace) -> pd.DataFrame:
    py = read_results(PY_CSV, "python")
    ftn = read_results(FTN_CSV, "fortran")
    merged = py.merge(ftn, on=["asset", "model"], suffixes=("_py", "_ftn"))
    rows: list[dict[str, object]] = []
    for _, row in merged.iterrows():
        for col in PARAM_COLS:
            py_val = float(row[f"{col}_py"])
            ftn_val = float(row[f"{col}_ftn"])
            diff = py_val - ftn_val
            denom = max(abs(py_val), abs(ftn_val), args.abs_tol)
            rows.append(
                {
                    "asset": row["asset"],
                    "model": row["model"],
                    "field": col,
                    "python": py_val,
                    "fortran": ftn_val,
                    "diff": diff,
                    "abs_diff": abs(diff),
                    "rel_diff": abs(diff) / denom,
                    "flag": significant_diff(py_val, ftn_val, col, args),
                }
            )
    return pd.DataFrame(rows)


def print_comparison(comp: pd.DataFrame, args: argparse.Namespace) -> None:
    flags = comp[comp["flag"]].copy()
    print("Comparison thresholds:")
    print(f"  parameters: abs_diff > {args.abs_tol:g} and rel_diff > {args.rel_tol:g}")
    print(f"  logL:       abs_diff > {args.logl_tol:g}")
    print()

    if flags.empty:
        print("No significant differences found under the configured thresholds.")
    else:
        print("Significant differences:")
        print(
            flags.sort_values(["asset", "model", "field"]).to_string(
                index=False,
                formatters={
                    "python": "{:.8g}".format,
                    "fortran": "{:.8g}".format,
                    "diff": "{:.8g}".format,
                    "abs_diff": "{:.8g}".format,
                    "rel_diff": "{:.4%}".format,
                },
            )
        )

    print("\nLargest absolute differences by field:")
    idx = comp.groupby("field")["abs_diff"].idxmax()
    print(
        comp.loc[idx].sort_values("field").to_string(
            index=False,
            formatters={
                "python": "{:.8g}".format,
                "fortran": "{:.8g}".format,
                "diff": "{:.8g}".format,
                "abs_diff": "{:.8g}".format,
                "rel_diff": "{:.4%}".format,
            },
        )
    )


def main() -> None:
    args = parse_args()
    py_elapsed, _, _, py_code = run_command(python_command(args), "Python arch")
    if py_code != 0 and not args.keep_going:
        raise SystemExit(py_code)

    ftn_elapsed, _, _, ftn_code = run_command([args.make_cmd, "run_fit_arch_common"], "Fortran")
    if ftn_code != 0 and not args.keep_going:
        raise SystemExit(ftn_code)

    comp = compare(args)
    args.comparison_output.parent.mkdir(parents=True, exist_ok=True)
    comp.to_csv(args.comparison_output, index=False)
    print("Runner wall times:")
    print(f"  Python arch: {py_elapsed:.3f} seconds")
    print(f"  Fortran:     {ftn_elapsed:.3f} seconds")
    print()
    print(f"Wrote {args.comparison_output}")
    print()
    print_comparison(comp, args)


if __name__ == "__main__":
    main()
