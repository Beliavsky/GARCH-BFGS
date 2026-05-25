import argparse

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def fit_asymmetric_u_shape(x, y):
    best = None
    best_sse = np.inf

    for center in np.linspace(x[1], x[-2], 200):
        left = np.where(x <= center, (x - center) ** 2, 0.0)
        right = np.where(x > center, (x - center) ** 2, 0.0)
        design = np.column_stack([np.ones_like(x), left, right])
        coef, *_ = np.linalg.lstsq(design, y, rcond=None)
        if coef[0] <= 0.0 or coef[1] < 0.0 or coef[2] < 0.0:
            continue
        fitted = design @ coef
        sse = np.sum((y - fitted) ** 2)
        if sse < best_sse:
            best_sse = sse
            best = (center, coef, fitted)

    if best is None:
        center = x[len(x) // 2]
        left = np.where(x <= center, (x - center) ** 2, 0.0)
        right = np.where(x > center, (x - center) ** 2, 0.0)
        design = np.column_stack([np.ones_like(x), left, right])
        coef, *_ = np.linalg.lstsq(design, y, rcond=None)
        best = (center, coef, design @ coef)

    return best


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", default="mcsgarch_diurnal_curve.csv")
    parser.add_argument("--fit-u-shape", action="store_true")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    x = np.arange(len(df))
    y = df["vol_mult"].to_numpy()

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(x, y, marker="o", linewidth=1.5, markersize=3, label="estimated")
    if args.fit_u_shape:
        center, coef, fitted = fit_asymmetric_u_shape(x.astype(float), y)
        ax.plot(x, fitted, linewidth=2.0, label="asymmetric U fit")
        center_idx = int(round(center))
        center_time = df.loc[min(max(center_idx, 0), len(df) - 1), "time"][:5]
        ax.set_title(f"MCS-GARCH Intraday Diurnal Volatility Curve (U minimum near {center_time})")
    else:
        ax.set_title("MCS-GARCH Intraday Diurnal Volatility Curve")
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)

    ax.set_xlabel("Time of day")
    ax.set_ylabel("Volatility multiplier")
    ax.grid(True, alpha=0.3)
    ax.legend()

    tick_step = max(1, len(df) // 10)
    tick_idx = np.arange(0, len(df), tick_step)
    if tick_idx[-1] != len(df) - 1:
        tick_idx = np.append(tick_idx, len(df) - 1)
    ax.set_xticks(tick_idx)
    tick_labels = df.loc[tick_idx, "time"].str.slice(0, 5)
    ax.set_xticklabels(tick_labels, rotation=45, ha="right")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
