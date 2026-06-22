# GARCH-BFGS

Fortran implementations of GARCH-family volatility models, distributional likelihoods, simulation programs, and command-line examples for fitting volatility models to close-to-close, OHLC, split overnight/intraday, and intraday return data.

The code is research-oriented. Most programs are small executable examples built from reusable modules, with model comparison based on log likelihood, AIC, BIC, and diagnostic summaries.

## Features

- Symmetric GARCH, NAGARCH, GJR-GARCH, EGARCH, and related flexible GARCH variants.
- Innovation distributions including normal, Student t, GED, Laplace, logistic, hyperbolic secant, NIG, and Fernandez-Steel skewed t in the newer distribution-aware fitters.
- Simulation routines for stationary close-to-close GARCH models.
- Close-to-close return fitting with normal and non-normal innovations.
- Two-step fitting that fits the volatility model under normal noise and then fits distributions to standardized residuals.
- Split close-open/open-close models for daily OHLC data.
- Intraday OHLCV readers and intraday MCS-GARCH style models with diurnal volatility curves.
- Simple deterministic intraday diurnal variance baselines independent of GARCH estimation.
- Compact unformatted stream storage for intraday OHLCV tick data, useful when CSV parsing dominates runtime.
- Intraday models with overnight effects, prior-day range predictors, and prior-day open-to-close NAGARCH news impact predictors.
- Intraday ACF diagnostics for signed price changes, absolute price changes, and high-low ranges.
- Cross-asset intraday covariance/correlation matrices at multiple aggregation scales, with diagnostics for pairs whose correlations change unusually across frequencies.
- Daily realized-volatility forecast comparisons using close-to-close, OHLC, realized measures, HAR-style models, HEAVY/realized-GARCH-style models, and implied-volatility correlations.
- Calendar-day annual seasonal volatility tests.
- Iid distribution fitting/simulation and normal-mixture EM examples.
- Univariate value-at-risk estimators including empirical, Harrell-Davis, kernel, parametric, EVT, and Monte Carlo variants.
- Basic DCC/ADCC, GAS, and stochastic volatility examples retained from earlier experiments.

## Requirements

- `gfortran`
- `make`

The Makefile currently uses:

```make
FC      = gfortran
FFLAGS  = -Wall -Wextra -Werror -Wno-compare-reals -fbounds-check -O2
```

On Windows, the examples have been run from PowerShell/Git Bash style shells using MinGW gfortran.

## Build

Build a specific executable:

```sh
make xfit_gen_garch_dist_returns.exe
```

Build all executable targets defined in the Makefile:

```sh
make all
```

Run a configured target:

```sh
make run_fit_gen_garch_dist_returns
```

Clean generated objects, modules, and executables:

```sh
make clean
```

Building `all` compiles all `.exe` targets currently defined in the Makefile and may take longer than building a specific executable.

## Data

Some example programs use hard-coded default input files. Update the file name in the program or pass a command-line argument where supported.

Common expected inputs:

- Daily adjusted close price CSVs, such as `spy_efa_eem_tlt_lqd.csv`.
- Daily OHLC price CSVs for split close-open/open-close or range models.
- Intraday OHLCV CSVs, such as 1-second or 5-minute bars with timestamp, open, high, low, close, and volume.

The example files `prices_ohlc.csv`, `spy_efa_eem_tlt_lqd.csv`, and `vix_spy.csv` were obtained from Yahoo Finance. Intraday files referenced in source-code defaults, such as `c:\python\databento\spy_1s_databento.csv` and `c:\python\intraday_prices\spy_5min_databento.csv`, are Databento data files and are not included. Databento currently offers free signup credits for historical market data.

Large market data files are not necessarily included in this repository.

## Main Modules

- `garch.f90`, `nagarch.f90`, `gjr.f90`, `egarch.f90`, `fgarch.f90`: model likelihood/filter logic.
- `garch_types.f90`: shared GARCH parameter/result types.
- `garch_fit.f90`: normal-noise fitting interface across GARCH model families.
- `garch_fit_dist.f90`: distribution-aware GARCH fitting.
- `garch_sim.f90`: simulation routines for stationary GARCH models.
- `garch_split_fit_dist.f90`: split close-open/open-close distribution-aware fitting.
- `garch_mcsgarch.f90`: reusable intraday MCS-GARCH filters and helpers.
- `intraday_vol_baseline.f90`: lagged/EWMA daily variance, deterministic diurnal multipliers, and simple intraday EWMA baseline forecasts.
- `intraday_summary.f90`: summary statistics and time-gap diagnostics for intraday OHLCV files.
- `intraday_realized_measures.f90`, `realized_vol_forecast.f90`, `realized_garch.f90`: daily realized-measure construction and realized-volatility forecast models.
- `distributions.f90`, `special.f90`, `random.f90`: distribution densities, special functions, and random variates.
- `market_data.f90`: intraday OHLCV containers, CSV and stream readers/writers, resampling, session filtering, and intraday transformations.
- `intraday_returns.f90`, `intraday_correlation_report.f90`, `matrix_print.f90`: intraday return alignment, cross-asset correlation reporting, and matrix printing.
- `csv.f90`, `date.f90`, `strings.f90`, `path_utils.f90`, `glob.f90`: data input and utility modules.
- `stats.f90`: sample statistics, correlations/covariances, autocorrelations, sorting, and ACF table printing.
- `seasonal_vol.f90`: calendar seasonal volatility regression helpers.
- `normal_mixture_em.f90`: normal-mixture EM code shared by the mixture examples.
- `var_univariate.f90`: univariate VaR estimators.

## Selected Programs

Close-to-close examples:

- `xfit_gen_garch_returns.f90`: fit several GARCH models with normal noise to returns.
- `xfit_gen_garch_dist_returns.f90`: fit several GARCH models with several innovation distributions.
- `xfit_garch_twostep_returns.f90`: fit GARCH under normal noise, then fit distributions to standardized residuals.
- `xsim_garch_fit.f90`: simulate GARCH processes, fit them back, and compare true versus fitted parameters.
- `xsim_garch_fit_dist.f90`: simulate and fit GARCH models with non-normal innovations.

OHLC and split-return examples:

- `xfit_split_garch_dist_returns.f90`: fit GARCH models to close-open and open-close returns.
- `xfit_split_range_garch_dist_returns.f90`: fit split models using close-open, open-close, and high-low information.
- `xfit_garch_ohlc_iv_returns.f90`: fit OHLC-based models and compare volatility forecasts.

Intraday examples:

- `xread_intraday_prices.f90`: read and summarize intraday OHLCV data.
- `xsummary_intraday.f90`: print summary statistics for one or more intraday OHLCV CSV or `.bin` files, with wildcard expansion.
- `xcsv_to_intraday_tick_stream.f90`: convert one or more intraday OHLCV CSV files, or all CSV files in a directory, to compact `.bin` stream files.
- `xroundtrip_intraday_tick_stream.f90`: test CSV to stream to stream-read round trips for intraday tick data.
- `xacf_intraday_measures.f90`: compute ACFs for signed close changes, absolute close changes, and high-low ranges; reads `.bin` files directly when given one.
- `xdiurnal_variance_baseline.f90`: estimate a deterministic time-of-day variance multiplier from intraday returns using `estimate_diurnal_variance_baseline`.
- `xfit_mcsgarch_intraday.f90`: fit intraday MCS-GARCH models with diurnal volatility curves.
- `xfit_mcsgarch_intraday_batch.f90`: run intraday MCS-GARCH fits over multiple files.
- `xfit_mcsgarch_on_intraday.f90`: fit joint overnight/intraday MCS-GARCH models.
- `xfit_mcsgarch_on_range_intraday.f90`: fit overnight/intraday models with prior-day range and open-to-close predictors.
- `xfit_mcsegarch_on_intraday.f90`: fit an EGARCH-style intraday model with overnight effects.
- `xcompare_intraday_ewma_ohlc.f90`: compare simple intraday EWMA baselines using close-close, Parkinson, and Garman-Klass proxies.
- `xcompare_intraday_ewma_freq.f90`: compare EWMA predictors from higher-frequency bars for lower-frequency target volatility.
- `xcorrel_intraday_assets.f90`: compute realized volatilities plus covariance/correlation matrices for multiple assets at several intraday aggregation scales; accepts individual files or `dir=...`.

Daily realized-volatility comparisons:

- `xcompare_daily_realized_vol_forecasts.f90`: compare daily volatility forecasts from realized measures, HAR-family models, realized-GARCH-style models, close-to-close GARCH/EWMA baselines, and optional implied-volatility correlations.
- `xcompare_daily_intraday_garch.f90`: compare daily forecasts from daily GARCH-style models and intraday MCS-GARCH-style models.
- `xcompare_regular_allhours_rv.f90`: compare regular-session and all-hours realized volatility as predictors of close-to-close volatility.

Diagnostics and utilities:

- `xseasonal_vol_calendar.f90`: test for annual calendar seasonal volatility.
- `xcompare_nagarch_news.f90`: compare NAGARCH news impact forms.
- `xfit_dist_returns.f90`: distribution fitting/testing utility for returns.
- `xfit_dist.f90`: fit iid distributions to numeric columns read from a CSV file.
- `xsim_dist.f90`: simulate iid distribution samples in the CSV format read by `xfit_dist.f90`.
- `xcalibrate_dist_warm_starts.f90`: calibrate distribution warm-start shape parameters.
- `xmix.f90`, `xmix_ic.f90`: normal-mixture EM simulation and information-criterion examples using `normal_mixture_em.f90`.
- `xvar_univariate.f90`: exercise the univariate VaR estimators on simulated returns.
- `xf90_make_deps.py`: report Makefile target dependencies by source file and audit `x*.f90` executable coverage.
- `xcompare_dirs.py`: compare selected source file types across two directories.

## Intraday Binary Stream Workflow

CSV parsing is slow for large intraday files. Convert a CSV once:

```sh
make xcsv_to_intraday_tick_stream.exe
./xcsv_to_intraday_tick_stream.exe c:\python\databento\spy_1s_databento.csv
```

This writes `spy_1s_databento.bin` in the current directory. Prices are stored as integer multiples of `tick_size`, which defaults to `0.001` dollars. Pass a second argument to use a different tick size, and a third argument to cap the number of data rows read:

```sh
./xcsv_to_intraday_tick_stream.exe c:\python\databento\spy_1s_databento.csv 0.005 100000
```

Named arguments are also supported and are clearer for batch conversion:

```sh
./xcsv_to_intraday_tick_stream.exe c:\python\databento\spy_1s_databento.csv tick_size=0.005 max_obs=100000
```

Convert all price CSV files in a directory and write the `.bin` files to a separate output directory:

```sh
./xcsv_to_intraday_tick_stream.exe dir=c:\python\intraday_prices\continuous out_dir=bin_dir tick_size=0.0000001
```

Programs can infer stream input from the `.bin` extension. For example:

```sh
make xacf_intraday_measures.exe
./xacf_intraday_measures.exe spy_1s_databento.bin
```

If `spy_1s_databento.bin` exists, `xacf_intraday_measures.exe` uses it by default; otherwise it falls back to the configured CSV path.

Summarize one or more intraday files:

```sh
make xsummary_intraday.exe
./xsummary_intraday.exe spy_1s_databento.bin
./xsummary_intraday.exe "c:\python\intraday_prices\continuous\*.csv"
```

The summary program reads CSV or `.bin` files, reports per-file OHLCV/return summaries, and can print time-gap diagnostics.

For cross-asset intraday correlations, pass files explicitly:

```sh
make xcorrel_intraday_assets.exe
./xcorrel_intraday_assets.exe ES.bin JY.bin TY.bin
```

or run on every `.bin` file in a directory:

```sh
./xcorrel_intraday_assets.exe dir=bin_dir
```

When a directory is supplied, `xcorrel_intraday_assets.exe` uses all `.bin` files if any are present; otherwise it uses all price CSV files in the directory. The output includes realized volatilities, covariance matrices, correlation matrices, and a final correlation-change diagnostic comparing lower-frequency correlations with the highest-frequency matrix using Fisher z differences.

Estimate a standalone deterministic intraday diurnal variance baseline:

```sh
make xdiurnal_variance_baseline.exe
./xdiurnal_variance_baseline.exe c:\python\intraday_prices\spy_5min_databento.csv diurnal_variance_baseline.csv
```

The program filters to regular-session bars, forms within-day log close-to-close returns, estimates lag-1 daily realized variance forecasts, calls `estimate_diurnal_variance_baseline`, prints the populated time-of-day bins, and writes the curve to CSV.

## Source and Makefile Utilities

List `.f90` files by how many Makefile targets depend on them:

```sh
python xf90_make_deps.py
python xf90_make_deps.py --show-targets
```

Audit `x*.f90` program files against matching `.exe` targets and `make all` coverage:

```sh
python xf90_make_deps.py --audit-x
```

Compare selected source files in two directories:

```sh
python xcompare_dirs.py dir1 dir2
```

By default it compares `*.f90`, `*.py`, `*.c`, `*.cpp`, `*.r`, and `*make*`, grouped by pattern.

## Intraday Range and OC Predictor Configuration

`xfit_mcsgarch_on_range_intraday.f90` has explicit configuration arrays near the top of the file. For example, active innovation distributions are controlled by:

```fortran
character(len=8), parameter :: fit_dist_names(*) = [character(len=8) :: &
    "NORMAL", "T"]
```

Add `"FS_SKEWT"` to restore Fernandez-Steel skewed-t fits:

```fortran
character(len=8), parameter :: fit_dist_names(*) = [character(len=8) :: &
    "NORMAL", "T", "FS_SKEWT"]
```

The same program compares:

- `R0`: no prior-day range predictor.
- `R5M`: prior-day sum of 5-minute Parkinson ranges.
- `RDAY`: prior-day full regular-session Parkinson range.
- `RBOTH`: both range predictors.
- `OC0`: no prior-day open-to-close news impact.
- `OCNAG`: prior-day open-to-close NAGARCH-style news impact.

## Current Notes

- The codebase contains many executable experiments. Prefer building one target at a time while developing.
- Some intraday model grids can be slow. Use normal-only distributions while developing model structure, then rerun a shortlist with Student t or skewed t.
- Several programs have default file paths from the local research environment. Treat these as examples and adjust them for your data layout.
- AIC and BIC comparisons are only directly meaningful when models are fit to the same observations and likelihood target.

