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
- Compact unformatted stream storage for intraday OHLCV tick data, useful when CSV parsing dominates runtime.
- Intraday models with overnight effects, prior-day range predictors, and prior-day open-to-close NAGARCH news impact predictors.
- Intraday ACF diagnostics for signed price changes, absolute price changes, and high-low ranges.
- Calendar-day annual seasonal volatility tests.
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

Run a configured target:

```sh
make run_fit_gen_garch_dist_returns
```

Clean generated objects, modules, and executables:

```sh
make clean
```

Building `all` will compile many experimental programs and may take longer than building a specific executable.

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
- `distributions.f90`, `special.f90`, `random.f90`: distribution densities, special functions, and random variates.
- `market_data.f90`: intraday OHLCV containers, CSV and stream readers/writers, resampling, session filtering, and intraday transformations.
- `csv.f90`, `date.f90`, `strings.f90`, `path_utils.f90`: data input and utility modules.
- `stats.f90`: sample statistics, autocorrelations, sorting, and ACF table printing.
- `seasonal_vol.f90`: calendar seasonal volatility regression helpers.

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
- `xcsv_to_intraday_tick_stream.f90`: convert an intraday OHLCV CSV to a compact `.bin` stream file in the current directory.
- `xroundtrip_intraday_tick_stream.f90`: test CSV to stream to stream-read round trips for intraday tick data.
- `xacf_intraday_measures.f90`: compute ACFs for signed close changes, absolute close changes, and high-low ranges; reads `.bin` files directly when given one.
- `xfit_mcsgarch_intraday.f90`: fit intraday MCS-GARCH models with diurnal volatility curves.
- `xfit_mcsgarch_intraday_batch.f90`: run intraday MCS-GARCH fits over multiple files.
- `xfit_mcsgarch_on_intraday.f90`: fit joint overnight/intraday MCS-GARCH models.
- `xfit_mcsgarch_on_range_intraday.f90`: fit overnight/intraday models with prior-day range and open-to-close predictors.
- `xfit_mcsegarch_on_intraday.f90`: fit an EGARCH-style intraday model with overnight effects.

Diagnostics and utilities:

- `xseasonal_vol_calendar.f90`: test for annual calendar seasonal volatility.
- `xcompare_nagarch_news.f90`: compare NAGARCH news impact forms.
- `xdist.f90`: distribution fitting/testing utility.

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

Programs can infer stream input from the `.bin` extension. For example:

```sh
make xacf_intraday_measures.exe
./xacf_intraday_measures.exe spy_1s_databento.bin
```

If `spy_1s_databento.bin` exists, `xacf_intraday_measures.exe` uses it by default; otherwise it falls back to the configured CSV path.

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

