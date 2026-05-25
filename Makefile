FC      = gfortran
FFLAGS  = -Wall -Wextra -Werror -Wno-compare-reals -fbounds-check -O2
EXE     = garch.exe
EXE2    = garch_scaling.exe
EXE3    = xgarch_t.exe
EXE4    = xgarch_dist.exe
EXE5    = xgarch_choose_dist.exe
EXE6    = xnagarch.exe
EXE7    = xgjr.exe
EXE8    = xegarch.exe
EXE9    = xgarch_flex.exe
EXE10   = xread_csv.exe
EXE11   = xfit_spy.exe
EXE12   = xfit_garch_returns.exe
EXE13   = xfit_garch_dist_returns.exe
EXE14   = xfit_gas_returns.exe
EXE15   = xgarch_gas.exe
EXE16   = xgas_scaling.exe
EXE17   = xgas_garch_scaling.exe
EXE18   = xsv_scaling.exe
EXE19   = xsv_lev_scaling.exe
EXE20   = xfit_sv_returns.exe
EXE21   = xsv_t_scaling.exe
EXE22   = xfit_sv_garch_returns.exe
EXE23   = xnagarch_mix.exe
EXE24   = xnagarch_mix_t.exe
EXE25   = xstgarch.exe
EXE26   = xarch_ew.exe
EXE27   = xarch_lw.exe
EXE28   = xarch.exe
EXE29   = xdcc.exe
EXE30   = xadcc.exe
EXE31   = xadcc_t.exe
EXE32   = xdcc_vt.exe
EXE33   = xadcc_vt.exe
EXE34   = xadcc_t_vt.exe
EXE35   = xdist.exe
EXE36   = xcompare_nagarch_news.exe
EXE37   = xfgarch.exe
EXE38   = xfit_nagarch_returns.exe
EXE39   = xfit_gjr_returns.exe
EXE40   = xfit_symm_garch_returns.exe
EXE41   = xfit_gen_garch_returns.exe
EXE42   = xfit_egarch_returns.exe
EXE43   = xfit_gen_egarch_returns.exe
EXE44   = xfit_gen_garch_iv_returns.exe
EXE45   = xfit_garch_ohlc_iv_returns.exe
EXE46   = xfit_nagarch_ohlc_returns.exe
EXE47   = xfit_split_ohlc_iv_returns.exe
EXE48   = xfit_split_range_ohlc_iv_returns.exe
EXE49   = xfit_arch_common.exe
EXE50   = xfit_symm_gen_garch_returns.exe
EXE51   = xfit_rugarch_builtin_symm_returns.exe
SHARED  = kind.o math_const.o garch.o bfgs.o
OBJS    = $(SHARED) garch_main.o
OBJS2   = $(SHARED) garch_scaling.o
OBJS3   = $(SHARED) special.o distributions.o random.o garch_t.o xgarch_t.o
OBJS4   = $(SHARED) special.o distributions.o random.o garch_t.o garch_sech.o garch_ged.o garch_laplace.o garch_logistic.o garch_nig.o xgarch_dist.o
OBJS5   = $(SHARED) special.o distributions.o random.o garch_t.o garch_sech.o garch_ged.o garch_laplace.o garch_logistic.o garch_nig.o garch_choose.o xgarch_choose_dist.o
OBJS6   = $(SHARED) nagarch.o xnagarch.o
OBJS7   = $(SHARED) gjr.o xgjr.o
OBJS8   = $(SHARED) egarch.o xegarch.o
OBJS9   = $(SHARED) special.o distributions.o random.o garch_t.o nagarch.o gjr.o egarch.o \
          garch_flex.o xgarch_flex.o
OBJS10  = kind.o date.o strings.o csv.o xread_csv.o
OBJS11  = kind.o math_const.o special.o distributions.o garch.o bfgs.o nagarch.o \
          gjr.o egarch.o garch_flex.o strings.o csv.o rank.o xfit_spy.o
OBJS12  = kind.o math_const.o special.o distributions.o garch.o bfgs.o nagarch.o \
          gjr.o egarch.o garch_flex.o strings.o csv.o stats.o rank.o xfit_garch_returns.o
OBJS13  = kind.o math_const.o special.o distributions.o garch.o bfgs.o nagarch.o \
          gjr.o egarch.o garch_flex.o strings.o csv.o stats.o rank.o xfit_garch_dist_returns.o
OBJS14  = kind.o math_const.o special.o distributions.o random.o bfgs.o gas.o strings.o csv.o stats.o rank.o xfit_gas_returns.o
OBJS15  = kind.o math_const.o special.o garch.o bfgs.o nagarch.o \
          gjr.o egarch.o garch_flex.o distributions.o random.o gas.o strings.o csv.o stats.o rank.o xgarch_gas.o
OBJS16  = kind.o math_const.o special.o distributions.o random.o bfgs.o gas.o xgas_scaling.o
OBJS17  = kind.o math_const.o special.o distributions.o random.o bfgs.o gas.o \
          garch.o nagarch.o gjr.o egarch.o garch_flex.o \
          xgas_garch_scaling.o
OBJS18  = kind.o math_const.o special.o distributions.o random.o bfgs.o sv.o xsv_scaling.o
OBJS19  = kind.o math_const.o special.o distributions.o random.o bfgs.o sv.o xsv_lev_scaling.o
OBJS20  = kind.o math_const.o special.o distributions.o random.o bfgs.o sv.o strings.o csv.o stats.o rank.o xfit_sv_returns.o
OBJS21  = kind.o math_const.o special.o distributions.o random.o bfgs.o sv.o xsv_t_scaling.o
OBJS22  = kind.o math_const.o special.o garch.o bfgs.o nagarch.o \
          gjr.o egarch.o garch_flex.o fgarch.o distributions.o random.o sv.o \
          strings.o csv.o stats.o rank.o xfit_sv_garch_returns.o
OBJS23  = kind.o math_const.o special.o distributions.o garch.o bfgs.o nagarch.o \
          gjr.o egarch.o garch_flex.o strings.o csv.o stats.o rank.o xnagarch_mix.o
OBJS24  = kind.o math_const.o special.o distributions.o garch.o bfgs.o nagarch.o \
          gjr.o egarch.o garch_flex.o strings.o csv.o stats.o rank.o xnagarch_mix_t.o
OBJS25  = kind.o math_const.o special.o distributions.o garch.o bfgs.o nagarch.o \
          gjr.o egarch.o garch_flex.o strings.o csv.o stats.o rank.o xstgarch.o
OBJS26  = kind.o math_const.o garch.o bfgs.o strings.o csv.o stats.o xarch_ew.o
OBJS27  = kind.o math_const.o garch.o bfgs.o strings.o csv.o stats.o xarch_lw.o
OBJS28  = kind.o math_const.o garch.o bfgs.o strings.o csv.o stats.o xarch.o
OBJS29  = kind.o math_const.o garch.o bfgs.o strings.o csv.o stats.o \
          linalg.o dcc.o xdcc.o
OBJS30  = kind.o math_const.o garch.o nagarch.o bfgs.o \
          strings.o csv.o stats.o linalg.o dcc.o xadcc.o
OBJS31  = kind.o math_const.o garch.o nagarch.o bfgs.o \
          strings.o csv.o stats.o linalg.o dcc.o xadcc_t.o
OBJS32  = kind.o math_const.o garch.o bfgs.o \
          strings.o csv.o stats.o linalg.o dcc.o xdcc_vt.o
OBJS33  = kind.o math_const.o garch.o nagarch.o bfgs.o \
          strings.o csv.o stats.o linalg.o dcc.o xadcc_vt.o
OBJS34  = kind.o math_const.o garch.o nagarch.o bfgs.o \
          strings.o csv.o stats.o linalg.o dcc.o xadcc_t_vt.o
OBJS35  = kind.o math_const.o special.o nagarch.o bfgs.o \
          strings.o csv.o stats.o distributions.o xdist.o
OBJS36  = kind.o math_const.o date.o strings.o csv.o nagarch.o bfgs.o xcompare_nagarch_news.o
OBJS37  = kind.o math_const.o garch_types.o garch.o nagarch.o gjr.o egarch.o fgarch.o garch_fit.o strings.o csv.o stats.o bfgs.o xfgarch.o
OBJS38  = kind.o math_const.o garch_types.o garch.o nagarch.o gjr.o egarch.o fgarch.o garch_fit.o strings.o csv.o stats.o bfgs.o xfit_nagarch_returns.o
OBJS39  = kind.o math_const.o garch_types.o garch.o nagarch.o gjr.o egarch.o fgarch.o garch_fit.o strings.o csv.o stats.o bfgs.o xfit_gjr_returns.o
OBJS40  = kind.o math_const.o garch_types.o garch.o nagarch.o gjr.o egarch.o fgarch.o garch_fit.o strings.o csv.o stats.o bfgs.o xfit_symm_garch_returns.o
OBJS41  = kind.o math_const.o garch_types.o strings.o date.o csv.o stats.o garch.o nagarch.o gjr.o egarch.o fgarch.o garch_fit.o garch_forecast.o model_selection.o bfgs.o xfit_gen_garch_returns.o
OBJS42  = kind.o math_const.o garch_types.o garch.o nagarch.o gjr.o egarch.o fgarch.o garch_fit.o strings.o csv.o stats.o bfgs.o xfit_egarch_returns.o
OBJS43  = kind.o math_const.o garch_types.o strings.o date.o csv.o stats.o garch.o nagarch.o gjr.o egarch.o fgarch.o garch_fit.o model_selection.o bfgs.o xfit_gen_egarch_returns.o
OBJS44  = kind.o math_const.o garch_types.o strings.o date.o csv.o stats.o garch.o nagarch.o gjr.o egarch.o fgarch.o garch_fit.o garch_forecast.o model_selection.o time_series_compare.o vol_forecast_compare.o bfgs.o xfit_gen_garch_iv_returns.o
OBJS45  = kind.o math_const.o garch_types.o strings.o date.o csv.o stats.o garch.o nagarch.o gjr.o egarch.o fgarch.o garch_fit.o garch_forecast.o model_selection.o time_series_compare.o vol_forecast_compare.o bfgs.o xfit_garch_ohlc_iv_returns.o
OBJS46  = kind.o math_const.o garch_types.o garch.o nagarch.o gjr.o egarch.o fgarch.o garch_fit.o strings.o csv.o stats.o bfgs.o xfit_nagarch_ohlc_returns.o
OBJS47  = kind.o math_const.o garch_types.o strings.o date.o csv.o stats.o nagarch.o time_series_compare.o vol_forecast_compare.o bfgs.o xfit_split_ohlc_iv_returns.o
OBJS48  = kind.o math_const.o garch_types.o strings.o date.o csv.o stats.o nagarch.o time_series_compare.o vol_forecast_compare.o bfgs.o xfit_split_range_ohlc_iv_returns.o
OBJS49  = kind.o math_const.o garch_types.o date.o strings.o csv.o stats.o garch.o nagarch.o gjr.o egarch.o fgarch.o garch_fit.o bfgs.o xfit_arch_common.o
OBJS50  = kind.o math_const.o garch_types.o strings.o date.o csv.o stats.o garch.o nagarch.o gjr.o egarch.o fgarch.o garch_fit.o garch_forecast.o model_selection.o bfgs.o xfit_symm_gen_garch_returns.o
OBJS51  = kind.o math_const.o garch_types.o date.o strings.o csv.o stats.o garch.o nagarch.o gjr.o egarch.o fgarch.o garch_fit.o garch_forecast.o model_selection.o bfgs.o xfit_rugarch_builtin_symm_returns.o

.PHONY: all run run_scaling run_t run_dist run_choose_dist run_nagarch run_fit_gas_returns run_gas_scaling run_gas_garch_scaling run_sv_scaling run_sv_lev_scaling run_fit_sv_returns run_sv_t_scaling run_fit_sv_garch_returns run_nagarch_mix run_nagarch_mix_t run_stgarch run_arch_ew run_arch_lw run_arch run_dcc run_adcc run_adcc_t run_dcc_vt run_adcc_vt run_adcc_t_vt run_xdist run_compare_nagarch_news run_fgarch run_fgarch_full run_fit_nagarch_returns run_fit_gjr_returns run_fit_symm_garch_returns run_fit_gen_garch_returns run_fit_egarch_returns run_fit_gen_egarch_returns run_fit_gen_garch_iv_returns run_fit_garch_ohlc_iv_returns run_fit_nagarch_ohlc_returns run_fit_split_ohlc_iv_returns run_fit_split_range_ohlc_iv_returns run_fit_arch_common run_fit_symm_gen_garch_returns run_fit_rugarch_builtin_symm_returns clean
.PHONY: run_xsim_garch_fit run_xsim_garch_fit_dist run_fit_gen_garch_dist_returns run_fit_garch_twostep_returns run_seasonal_vol_calendar run_fit_split_garch_dist_returns run_fit_split_range_garch_dist_returns run_xmcsgarch run_xread_intraday_prices run_csv_to_intraday_tick_stream run_roundtrip_intraday_tick_stream run_fit_mcsgarch_intraday

all: $(EXE) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11) $(EXE12) $(EXE13) $(EXE14) $(EXE15) $(EXE16) $(EXE17) $(EXE18) $(EXE19) $(EXE20) $(EXE21) $(EXE22) $(EXE23) $(EXE24) $(EXE25) $(EXE26) $(EXE27) $(EXE28) $(EXE29) $(EXE30) $(EXE31) $(EXE32) $(EXE33) $(EXE34) $(EXE35) $(EXE36) $(EXE37) $(EXE38) $(EXE39) $(EXE40) $(EXE41) $(EXE42) $(EXE43) $(EXE44) $(EXE45) $(EXE46) $(EXE47) $(EXE48) $(EXE49) $(EXE50) $(EXE51)

$(EXE): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE2): $(OBJS2)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE3): $(OBJS3)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE4): $(OBJS4)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE5): $(OBJS5)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE6): $(OBJS6)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE7): $(OBJS7)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE8): $(OBJS8)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE9): $(OBJS9)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE10): $(OBJS10)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE11): $(OBJS11)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE12): $(OBJS12)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE13): $(OBJS13)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE14): $(OBJS14)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE15): $(OBJS15)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE16): $(OBJS16)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE17): $(OBJS17)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE18): $(OBJS18)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE19): $(OBJS19)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE20): $(OBJS20)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE21): $(OBJS21)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE22): $(OBJS22)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE23): $(OBJS23)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE24): $(OBJS24)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE25): $(OBJS25)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE26): $(OBJS26)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE27): $(OBJS27)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE28): $(OBJS28)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE29): $(OBJS29)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE30): $(OBJS30)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE31): $(OBJS31)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE32): $(OBJS32)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE33): $(OBJS33)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE34): $(OBJS34)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE35): $(OBJS35)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE36): $(OBJS36)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE37): $(OBJS37)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE38): $(OBJS38)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE39): $(OBJS39)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE40): $(OBJS40)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE41): $(OBJS41)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE42): $(OBJS42)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE43): $(OBJS43)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE44): $(OBJS44)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE45): $(OBJS45)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE46): $(OBJS46)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE47): $(OBJS47)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE48): $(OBJS48)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE49): $(OBJS49)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE50): $(OBJS50)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE51): $(OBJS51)
	$(FC) $(FFLAGS) -o $@ $^

xsim_garch_fit.exe: xsim_garch_fit.f90 kind.o math_const.o garch_types.o garch.o nagarch.o gjr.o egarch.o fgarch.o garch_fit.o garch_sim.o bfgs.o
	$(FC) $(FFLAGS) -o $@ $^

xsim_garch_fit_dist.exe: xsim_garch_fit_dist.f90 kind.o math_const.o stats.o special.o distributions.o random.o \
                         garch.o garch_t.o garch_sech.o garch_ged.o garch_laplace.o \
                         garch_logistic.o garch_nig.o garch_types.o nagarch.o gjr.o egarch.o \
                         fgarch.o garch_fit.o garch_fit_dist.o bfgs.o
	$(FC) $(FFLAGS) -o $@ $^

xfit_gen_garch_dist_returns.exe: xfit_gen_garch_dist_returns.f90 kind.o date.o strings.o csv.o stats.o garch_types.o special.o \
                                 distributions.o garch.o nagarch.o gjr.o egarch.o fgarch.o \
                                 garch_fit.o garch_fit_dist.o bfgs.o
	$(FC) $(FFLAGS) -o $@ $^

xfit_garch_twostep_returns.exe: xfit_garch_twostep_returns.f90 kind.o date.o strings.o csv.o stats.o garch_types.o special.o \
                                distributions.o garch.o nagarch.o gjr.o egarch.o fgarch.o \
                                garch_fit.o garch_fit_dist.o bfgs.o
	$(FC) $(FFLAGS) -o $@ $^

xseasonal_vol_calendar.exe: xseasonal_vol_calendar.f90 kind.o date.o strings.o csv.o stats.o linalg.o seasonal_vol.o
	$(FC) $(FFLAGS) -o $@ $^

xfit_split_garch_dist_returns.exe: xfit_split_garch_dist_returns.f90 kind.o date.o strings.o csv.o stats.o math_const.o \
                                   garch_types.o special.o distributions.o garch_fit_dist.o \
                                   garch_split_fit_dist.o garch.o nagarch.o gjr.o egarch.o fgarch.o \
                                   garch_fit.o bfgs.o
	$(FC) $(FFLAGS) -o $@ $^

xfit_split_range_garch_dist_returns.exe: xfit_split_range_garch_dist_returns.f90 kind.o date.o strings.o csv.o stats.o math_const.o \
                                         garch_types.o special.o distributions.o garch_fit_dist.o \
                                         garch_split_fit_dist.o garch.o nagarch.o gjr.o egarch.o fgarch.o \
                                         garch_fit.o bfgs.o
	$(FC) $(FFLAGS) -o $@ $^

xmcsgarch.exe: xmcsgarch.f90 kind.o math_const.o special.o distributions.o stats.o bfgs.o garch_mcsgarch.o
	$(FC) $(FFLAGS) -o $@ $^

xread_intraday_prices.exe: xread_intraday_prices.f90 kind.o date.o strings.o market_data.o
	$(FC) $(FFLAGS) -o $@ $^

xcheck_intraday_penny_grid.exe: xcheck_intraday_penny_grid.f90 kind.o date.o strings.o market_data.o
	$(FC) $(FFLAGS) -o $@ $^

xroundtrip_intraday_tick_stream.exe: xroundtrip_intraday_tick_stream.f90 kind.o date.o strings.o market_data.o
	$(FC) $(FFLAGS) -o $@ $^

xcsv_to_intraday_tick_stream.exe: xcsv_to_intraday_tick_stream.f90 kind.o date.o strings.o market_data.o path_utils.o
	$(FC) $(FFLAGS) -o $@ $^

xfit_intraday_tick_returns.exe: xfit_intraday_tick_returns.f90 kind.o math_const.o date.o strings.o market_data.o
	$(FC) $(FFLAGS) -o $@ $^

xacf_abs_price_changes.exe: xacf_abs_price_changes.f90 kind.o date.o strings.o market_data.o
	$(FC) $(FFLAGS) -o $@ $^

xacf_price_changes.exe: xacf_price_changes.f90 kind.o date.o strings.o market_data.o
	$(FC) $(FFLAGS) -o $@ $^

xacf_price_ranges.exe: xacf_price_ranges.f90 kind.o date.o strings.o market_data.o
	$(FC) $(FFLAGS) -o $@ $^

xacf_intraday_measures.exe: xacf_intraday_measures.f90 kind.o date.o strings.o market_data.o path_utils.o stats.o
	$(FC) $(FFLAGS) -o $@ $^

xcorr_signed_future_abs_price_changes.exe: xcorr_signed_future_abs_price_changes.f90 kind.o date.o strings.o market_data.o
	$(FC) $(FFLAGS) -o $@ $^

xfit_mcsgarch_intraday.exe: xfit_mcsgarch_intraday.f90 kind.o math_const.o special.o distributions.o date.o strings.o market_data.o stats.o bfgs.o garch_mcsgarch.o intraday_vol_baseline.o
	$(FC) $(FFLAGS) -o $@ $^

xfit_mcsgarch_intraday_batch.exe: xfit_mcsgarch_intraday_batch.f90 kind.o math_const.o special.o distributions.o date.o strings.o market_data.o stats.o bfgs.o garch_mcsgarch.o intraday_vol_baseline.o
	$(FC) $(FFLAGS) -o $@ $^

xfit_mcsgarch_on_intraday.exe: xfit_mcsgarch_on_intraday.f90 kind.o math_const.o special.o distributions.o date.o strings.o market_data.o stats.o bfgs.o garch_mcsgarch.o intraday_vol_baseline.o
	$(FC) $(FFLAGS) -o $@ $^

xfit_mcsgarch_on_range_intraday.exe: xfit_mcsgarch_on_range_intraday.f90 kind.o math_const.o special.o distributions.o date.o strings.o market_data.o stats.o bfgs.o garch_mcsgarch.o intraday_vol_baseline.o
	$(FC) $(FFLAGS) -o $@ $^

xfit_mcsegarch_on_intraday.exe: xfit_mcsegarch_on_intraday.f90 kind.o math_const.o special.o distributions.o date.o strings.o market_data.o stats.o bfgs.o garch_mcsgarch.o
	$(FC) $(FFLAGS) -o $@ $^

xcompare_intraday_ewma_ohlc.exe: xcompare_intraday_ewma_ohlc.f90 kind.o math_const.o date.o strings.o market_data.o intraday_vol_baseline.o
	$(FC) $(FFLAGS) -o $@ $^

xcompare_intraday_ewma_freq.exe: xcompare_intraday_ewma_freq.f90 kind.o math_const.o date.o strings.o market_data.o intraday_vol_baseline.o
	$(FC) $(FFLAGS) -o $@ $^

kind.o: kind.f90
	$(FC) $(FFLAGS) -c $<

strings.o: strings.f90
	$(FC) $(FFLAGS) -c $<

garch_types.o: garch_types.f90 kind.o
	$(FC) $(FFLAGS) -c $<

math_const.o: math_const.f90 kind.o
	$(FC) $(FFLAGS) -c $<

garch.o: garch.f90 kind.o math_const.o
	$(FC) $(FFLAGS) -c $<

bfgs.o: bfgs.f90 kind.o
	$(FC) $(FFLAGS) -c $<

garch_main.o: garch_main.f90 kind.o garch.o bfgs.o
	$(FC) $(FFLAGS) -c $<

garch_scaling.o: garch_scaling.f90 kind.o garch.o bfgs.o
	$(FC) $(FFLAGS) -c $<

nagarch.o: nagarch.f90 kind.o math_const.o
	$(FC) $(FFLAGS) -c $<

garch_t.o: garch_t.f90 kind.o math_const.o garch.o special.o random.o
	$(FC) $(FFLAGS) -c $<

garch_sech.o: garch_sech.f90 kind.o math_const.o garch.o
	$(FC) $(FFLAGS) -c $<

xgarch_t.o: xgarch_t.f90 kind.o garch.o garch_t.o bfgs.o
	$(FC) $(FFLAGS) -c $<

garch_ged.o: garch_ged.f90 kind.o garch.o special.o random.o distributions.o
	$(FC) $(FFLAGS) -c $<

garch_laplace.o: garch_laplace.f90 kind.o math_const.o garch.o
	$(FC) $(FFLAGS) -c $<

garch_logistic.o: garch_logistic.f90 kind.o math_const.o garch.o
	$(FC) $(FFLAGS) -c $<

garch_nig.o: garch_nig.f90 kind.o math_const.o garch.o special.o random.o
	$(FC) $(FFLAGS) -c $<

xgarch_dist.o: xgarch_dist.f90 kind.o garch.o garch_t.o garch_sech.o garch_ged.o garch_laplace.o garch_logistic.o garch_nig.o bfgs.o
	$(FC) $(FFLAGS) -c $<

garch_choose.o: garch_choose.f90 kind.o
	$(FC) $(FFLAGS) -c $<

xgarch_choose_dist.o: xgarch_choose_dist.f90 kind.o garch_choose.o garch.o garch_t.o garch_sech.o garch_ged.o garch_laplace.o garch_logistic.o garch_nig.o bfgs.o
	$(FC) $(FFLAGS) -c $<

xnagarch.o: xnagarch.f90 kind.o nagarch.o bfgs.o
	$(FC) $(FFLAGS) -c $<

gjr.o: gjr.f90 kind.o math_const.o
	$(FC) $(FFLAGS) -c $<

xgjr.o: xgjr.f90 kind.o gjr.o bfgs.o
	$(FC) $(FFLAGS) -c $<

egarch.o: egarch.f90 kind.o math_const.o
	$(FC) $(FFLAGS) -c $<

xegarch.o: xegarch.f90 kind.o egarch.o bfgs.o
	$(FC) $(FFLAGS) -c $<

garch_sim.o: garch_sim.f90 kind.o garch_types.o garch.o nagarch.o gjr.o egarch.o
	$(FC) $(FFLAGS) -c $<

garch_mcsgarch.o: garch_mcsgarch.f90 kind.o math_const.o distributions.o stats.o bfgs.o
	$(FC) $(FFLAGS) -c $<

intraday_vol_baseline.o: intraday_vol_baseline.f90 kind.o
	$(FC) $(FFLAGS) -c $<

garch_flex.o: garch_flex.f90 kind.o math_const.o garch.o \
                  nagarch.o gjr.o egarch.o special.o distributions.o
	$(FC) $(FFLAGS) -c $<

xgarch_flex.o: xgarch_flex.f90 kind.o garch_flex.o garch_t.o \
               gjr.o egarch.o bfgs.o
	$(FC) $(FFLAGS) -c $<

date.o: date.f90
	$(FC) $(FFLAGS) -c $<

csv.o: csv.f90 kind.o date.o strings.o
	$(FC) $(FFLAGS) -c $<

market_data.o: market_data.f90 kind.o date.o strings.o
	$(FC) $(FFLAGS) -c $<

path_utils.o: path_utils.f90 strings.o
	$(FC) $(FFLAGS) -c $<

gas.o: gas.f90 kind.o math_const.o random.o
	$(FC) $(FFLAGS) -c $<

stats.o: stats.f90 kind.o
	$(FC) $(FFLAGS) -c $<

rank.o: rank.f90 kind.o
	$(FC) $(FFLAGS) -c $<

seasonal_vol.o: seasonal_vol.f90 kind.o stats.o linalg.o
	$(FC) $(FFLAGS) -c $<

xread_csv.o: xread_csv.f90 kind.o date.o strings.o csv.o
	$(FC) $(FFLAGS) -c $<

xfit_spy.o: xfit_spy.f90 kind.o garch_flex.o garch.o nagarch.o \
            gjr.o egarch.o bfgs.o strings.o csv.o rank.o
	$(FC) $(FFLAGS) -c $<

xfit_garch_returns.o: xfit_garch_returns.f90 kind.o garch_flex.o garch.o \
                      nagarch.o gjr.o egarch.o bfgs.o strings.o csv.o rank.o
	$(FC) $(FFLAGS) -c $<

xfit_garch_dist_returns.o: xfit_garch_dist_returns.f90 kind.o garch_flex.o garch.o \
                           nagarch.o gjr.o egarch.o bfgs.o strings.o csv.o rank.o
	$(FC) $(FFLAGS) -c $<

xfit_gas_returns.o: xfit_gas_returns.f90 kind.o gas.o bfgs.o strings.o csv.o stats.o rank.o
	$(FC) $(FFLAGS) -c $<

xgarch_gas.o: xgarch_gas.f90 kind.o garch_flex.o garch.o nagarch.o \
              gjr.o egarch.o gas.o bfgs.o strings.o csv.o stats.o rank.o
	$(FC) $(FFLAGS) -c $<

xgas_scaling.o: xgas_scaling.f90 kind.o gas.o bfgs.o
	$(FC) $(FFLAGS) -c $<

xgas_garch_scaling.o: xgas_garch_scaling.f90 kind.o gas.o garch_flex.o \
                      garch.o bfgs.o
	$(FC) $(FFLAGS) -c $<

sv.o: sv.f90 kind.o math_const.o special.o random.o
	$(FC) $(FFLAGS) -c $<

xsv_scaling.o: xsv_scaling.f90 kind.o sv.o bfgs.o
	$(FC) $(FFLAGS) -c $<

xsv_lev_scaling.o: xsv_lev_scaling.f90 kind.o sv.o bfgs.o
	$(FC) $(FFLAGS) -c $<

xfit_sv_returns.o: xfit_sv_returns.f90 kind.o sv.o bfgs.o strings.o csv.o stats.o rank.o
	$(FC) $(FFLAGS) -c $<

xsv_t_scaling.o: xsv_t_scaling.f90 kind.o sv.o bfgs.o
	$(FC) $(FFLAGS) -c $<

fgarch.o: fgarch.f90 kind.o math_const.o
	$(FC) $(FFLAGS) -c $<

garch_fit.o: garch_fit.f90 kind.o math_const.o garch_types.o garch.o nagarch.o gjr.o egarch.o fgarch.o bfgs.o
	$(FC) $(FFLAGS) -c $<

garch_fit_dist.o: garch_fit_dist.f90 kind.o garch_types.o distributions.o stats.o garch_fit.o bfgs.o
	$(FC) $(FFLAGS) -c $<

garch_split_fit_dist.o: garch_split_fit_dist.f90 kind.o math_const.o garch_types.o distributions.o garch_fit_dist.o bfgs.o
	$(FC) $(FFLAGS) -c $<

garch_forecast.o: garch_forecast.f90 kind.o date.o strings.o csv.o stats.o garch_types.o garch_fit.o
	$(FC) $(FFLAGS) -c $<

model_selection.o: model_selection.f90 kind.o
	$(FC) $(FFLAGS) -c $<

time_series_compare.o: time_series_compare.f90 kind.o strings.o
	$(FC) $(FFLAGS) -c $<

vol_forecast_compare.o: vol_forecast_compare.f90 kind.o strings.o date.o csv.o time_series_compare.o
	$(FC) $(FFLAGS) -c $<

xfit_sv_garch_returns.o: xfit_sv_garch_returns.f90 kind.o sv.o garch_flex.o \
                          fgarch.o garch.o nagarch.o gjr.o \
                          egarch.o bfgs.o strings.o csv.o stats.o rank.o
	$(FC) $(FFLAGS) -c $<

xnagarch_mix.o: xnagarch_mix.f90 kind.o math_const.o garch_flex.o \
                nagarch.o bfgs.o strings.o csv.o stats.o rank.o
	$(FC) $(FFLAGS) -c $<

xnagarch_mix_t.o: xnagarch_mix_t.f90 kind.o math_const.o garch_flex.o \
                  nagarch.o bfgs.o strings.o csv.o stats.o rank.o
	$(FC) $(FFLAGS) -c $<

xstgarch.o: xstgarch.f90 kind.o math_const.o garch_flex.o nagarch.o \
            bfgs.o strings.o csv.o stats.o rank.o
	$(FC) $(FFLAGS) -c $<

xarch_ew.o: xarch_ew.f90 kind.o math_const.o garch.o bfgs.o \
            strings.o csv.o stats.o
	$(FC) $(FFLAGS) -c $<

xarch_lw.o: xarch_lw.f90 kind.o math_const.o garch.o bfgs.o \
            strings.o csv.o stats.o
	$(FC) $(FFLAGS) -c $<

xarch.o: xarch.f90 kind.o math_const.o garch.o bfgs.o \
         strings.o csv.o stats.o
	$(FC) $(FFLAGS) -c $<

linalg.o: linalg.f90 kind.o
	$(FC) $(FFLAGS) -c $<

dcc.o: dcc.f90 kind.o linalg.o
	$(FC) $(FFLAGS) -c $<

xdcc.o: xdcc.f90 kind.o garch.o dcc.o bfgs.o strings.o csv.o stats.o
	$(FC) $(FFLAGS) -c $<

xadcc.o: xadcc.f90 kind.o nagarch.o dcc.o bfgs.o strings.o csv.o stats.o
	$(FC) $(FFLAGS) -c $<

xadcc_t.o: xadcc_t.f90 kind.o nagarch.o dcc.o bfgs.o strings.o csv.o stats.o
	$(FC) $(FFLAGS) -c $<

xdcc_vt.o: xdcc_vt.f90 kind.o garch.o dcc.o bfgs.o strings.o csv.o stats.o
	$(FC) $(FFLAGS) -c $<

xadcc_vt.o: xadcc_vt.f90 kind.o nagarch.o dcc.o bfgs.o strings.o csv.o stats.o
	$(FC) $(FFLAGS) -c $<

xadcc_t_vt.o: xadcc_t_vt.f90 kind.o nagarch.o dcc.o bfgs.o strings.o csv.o stats.o
	$(FC) $(FFLAGS) -c $<

random.o: random.f90 kind.o math_const.o distributions.o
	$(FC) $(FFLAGS) -c $<

special.o: special.f90 kind.o math_const.o
	$(FC) $(FFLAGS) -c $<

distributions.o: distributions.f90 kind.o math_const.o bfgs.o special.o
	$(FC) $(FFLAGS) -c $<

xdist.o: xdist.f90 kind.o nagarch.o distributions.o bfgs.o strings.o csv.o stats.o
	$(FC) $(FFLAGS) -c $<

xcompare_nagarch_news.o: xcompare_nagarch_news.f90 kind.o date.o strings.o csv.o nagarch.o bfgs.o
	$(FC) $(FFLAGS) -c $<

xfgarch.o: xfgarch.f90 kind.o garch_types.o fgarch.o garch_fit.o strings.o csv.o stats.o bfgs.o
	$(FC) $(FFLAGS) -c $<

xfit_nagarch_returns.o: xfit_nagarch_returns.f90 kind.o date.o strings.o csv.o stats.o nagarch.o garch_types.o garch_fit.o
	$(FC) $(FFLAGS) -c $<

xfit_gjr_returns.o: xfit_gjr_returns.f90 kind.o date.o strings.o csv.o stats.o garch_types.o garch_fit.o
	$(FC) $(FFLAGS) -c $<

xfit_symm_garch_returns.o: xfit_symm_garch_returns.f90 kind.o date.o strings.o csv.o stats.o garch_types.o garch_fit.o
	$(FC) $(FFLAGS) -c $<

xfit_gen_garch_returns.o: xfit_gen_garch_returns.f90 kind.o strings.o date.o csv.o stats.o nagarch.o garch_types.o garch_fit.o garch_forecast.o model_selection.o
	$(FC) $(FFLAGS) -c $<

xfit_egarch_returns.o: xfit_egarch_returns.f90 kind.o date.o strings.o csv.o stats.o garch_types.o garch_fit.o
	$(FC) $(FFLAGS) -c $<

xfit_gen_egarch_returns.o: xfit_gen_egarch_returns.f90 kind.o strings.o date.o csv.o stats.o nagarch.o garch_types.o garch_fit.o model_selection.o
	$(FC) $(FFLAGS) -c $<

xfit_gen_garch_iv_returns.o: xfit_gen_garch_iv_returns.f90 kind.o strings.o date.o csv.o stats.o nagarch.o garch_types.o garch_fit.o garch_forecast.o model_selection.o vol_forecast_compare.o
	$(FC) $(FFLAGS) -c $<

xfit_garch_ohlc_iv_returns.o: xfit_garch_ohlc_iv_returns.f90 kind.o strings.o date.o csv.o stats.o nagarch.o garch_types.o garch_fit.o garch_forecast.o model_selection.o vol_forecast_compare.o
	$(FC) $(FFLAGS) -c $<

xfit_nagarch_ohlc_returns.o: xfit_nagarch_ohlc_returns.f90 kind.o date.o strings.o csv.o stats.o nagarch.o garch_types.o garch_fit.o
	$(FC) $(FFLAGS) -c $<

xfit_split_ohlc_iv_returns.o: xfit_split_ohlc_iv_returns.f90 kind.o math_const.o strings.o date.o csv.o stats.o nagarch.o vol_forecast_compare.o bfgs.o
	$(FC) $(FFLAGS) -c $<

xfit_split_range_ohlc_iv_returns.o: xfit_split_range_ohlc_iv_returns.f90 kind.o math_const.o strings.o date.o csv.o stats.o nagarch.o vol_forecast_compare.o bfgs.o
	$(FC) $(FFLAGS) -c $<

xfit_arch_common.o: xfit_arch_common.f90 kind.o math_const.o date.o strings.o csv.o stats.o garch_types.o garch_fit.o bfgs.o
	$(FC) $(FFLAGS) -c $<

xfit_symm_gen_garch_returns.o: xfit_symm_gen_garch_returns.f90 kind.o strings.o date.o csv.o stats.o garch_types.o garch_fit.o garch_forecast.o model_selection.o
	$(FC) $(FFLAGS) -c $<

xfit_rugarch_builtin_symm_returns.o: xfit_rugarch_builtin_symm_returns.f90 kind.o date.o strings.o csv.o stats.o garch_types.o garch_fit.o garch_forecast.o model_selection.o
	$(FC) $(FFLAGS) -c $<

run: $(EXE)
	./$(EXE)

run_scaling: $(EXE2)
	./$(EXE2)

run_t: $(EXE3)
	./$(EXE3)

run_dist: $(EXE4)
	./$(EXE4)

run_choose_dist: $(EXE5)
	./$(EXE5)

run_nagarch: $(EXE6)
	./$(EXE6)

run_fit_gas_returns: $(EXE14)
	./$(EXE14)

run_gas_scaling: $(EXE16)
	./$(EXE16)

run_gas_garch_scaling: $(EXE17)
	./$(EXE17)

run_sv_scaling: $(EXE18)
	./$(EXE18)

run_sv_lev_scaling: $(EXE19)
	./$(EXE19)

run_fit_sv_returns: $(EXE20)
	./$(EXE20)

run_sv_t_scaling: $(EXE21)
	./$(EXE21)

run_fit_sv_garch_returns: $(EXE22)
	./$(EXE22)

run_nagarch_mix: $(EXE23)
	./$(EXE23)

run_nagarch_mix_t: $(EXE24)
	./$(EXE24)

run_stgarch: $(EXE25)
	./$(EXE25)

run_arch_ew: $(EXE26)
	./$(EXE26)

run_arch_lw: $(EXE27)
	./$(EXE27)

run_arch: $(EXE28)
	./$(EXE28)

run_dcc: $(EXE29)
	./$(EXE29)

run_adcc: $(EXE30)
	./$(EXE30)

run_adcc_t: $(EXE31)
	./$(EXE31)

run_dcc_vt: $(EXE32)
	./$(EXE32)

run_adcc_vt: $(EXE33)
	./$(EXE33)

run_adcc_t_vt: $(EXE34)
	./$(EXE34)

run_xdist: $(EXE35)
	./$(EXE35)

run_compare_nagarch_news: $(EXE36)
	./$(EXE36)

run_fgarch: $(EXE37)
	./$(EXE37)

run_fgarch_full: $(EXE37)
	./$(EXE37) --full

run_fit_nagarch_returns: $(EXE38)
	./$(EXE38)

run_fit_gjr_returns: $(EXE39)
	./$(EXE39)

run_fit_symm_garch_returns: $(EXE40)
	./$(EXE40)

run_fit_gen_garch_returns: $(EXE41)
	./$(EXE41)

run_fit_egarch_returns: $(EXE42)
	./$(EXE42)

run_fit_gen_egarch_returns: $(EXE43)
	./$(EXE43)

run_fit_gen_garch_iv_returns: $(EXE44)
	./$(EXE44)

run_fit_garch_ohlc_iv_returns: $(EXE45)
	./$(EXE45)

run_fit_nagarch_ohlc_returns: $(EXE46)
	./$(EXE46)

run_fit_split_ohlc_iv_returns: $(EXE47)
	./$(EXE47)

run_fit_split_range_ohlc_iv_returns: $(EXE48)
	./$(EXE48)

run_fit_arch_common: $(EXE49)
	./$(EXE49)

run_fit_symm_gen_garch_returns: $(EXE50)
	./$(EXE50)

run_fit_rugarch_builtin_symm_returns: $(EXE51)
	./$(EXE51)

run_xsim_garch_fit: xsim_garch_fit.exe
	./xsim_garch_fit.exe

run_xsim_garch_fit_dist: xsim_garch_fit_dist.exe
	./xsim_garch_fit_dist.exe

run_fit_gen_garch_dist_returns: xfit_gen_garch_dist_returns.exe
	./xfit_gen_garch_dist_returns.exe

run_fit_garch_twostep_returns: xfit_garch_twostep_returns.exe
	./xfit_garch_twostep_returns.exe

run_seasonal_vol_calendar: xseasonal_vol_calendar.exe
	./xseasonal_vol_calendar.exe

run_fit_split_garch_dist_returns: xfit_split_garch_dist_returns.exe
	./xfit_split_garch_dist_returns.exe

run_fit_split_range_garch_dist_returns: xfit_split_range_garch_dist_returns.exe
	./xfit_split_range_garch_dist_returns.exe

run_xmcsgarch: xmcsgarch.exe
	./xmcsgarch.exe

run_xread_intraday_prices: xread_intraday_prices.exe
	./xread_intraday_prices.exe

run_check_intraday_penny_grid: xcheck_intraday_penny_grid.exe
	./xcheck_intraday_penny_grid.exe

run_roundtrip_intraday_tick_stream: xroundtrip_intraday_tick_stream.exe
	./xroundtrip_intraday_tick_stream.exe

run_csv_to_intraday_tick_stream: xcsv_to_intraday_tick_stream.exe
	./xcsv_to_intraday_tick_stream.exe

run_fit_intraday_tick_returns: xfit_intraday_tick_returns.exe
	./xfit_intraday_tick_returns.exe

run_acf_abs_price_changes: xacf_abs_price_changes.exe
	./xacf_abs_price_changes.exe

run_acf_price_changes: xacf_price_changes.exe
	./xacf_price_changes.exe

run_acf_price_ranges: xacf_price_ranges.exe
	./xacf_price_ranges.exe

run_acf_intraday_measures: xacf_intraday_measures.exe
	./xacf_intraday_measures.exe

run_corr_signed_future_abs_price_changes: xcorr_signed_future_abs_price_changes.exe
	./xcorr_signed_future_abs_price_changes.exe

run_fit_mcsgarch_intraday: xfit_mcsgarch_intraday.exe
	./xfit_mcsgarch_intraday.exe

run_fit_mcsgarch_intraday_batch: xfit_mcsgarch_intraday_batch.exe
	./xfit_mcsgarch_intraday_batch.exe

run_fit_mcsgarch_on_intraday: xfit_mcsgarch_on_intraday.exe
	./xfit_mcsgarch_on_intraday.exe

run_fit_mcsgarch_on_range_intraday: xfit_mcsgarch_on_range_intraday.exe
	./xfit_mcsgarch_on_range_intraday.exe

run_fit_mcsegarch_on_intraday: xfit_mcsegarch_on_intraday.exe
	./xfit_mcsegarch_on_intraday.exe

run_compare_intraday_ewma_ohlc: xcompare_intraday_ewma_ohlc.exe
	./xcompare_intraday_ewma_ohlc.exe

run_compare_intraday_ewma_freq: xcompare_intraday_ewma_freq.exe
	./xcompare_intraday_ewma_freq.exe

clean:
	rm -f random.o special.o date.o market_data.o path_utils.o csv_to_intraday_tick_stream.o garch_mcsgarch.o intraday_vol_baseline.o xmcsgarch.exe xread_intraday_prices.exe xcheck_intraday_penny_grid.exe xroundtrip_intraday_tick_stream.exe xcsv_to_intraday_tick_stream.exe xfit_intraday_tick_returns.exe xacf_abs_price_changes.exe xacf_price_changes.exe xacf_price_ranges.exe xacf_intraday_measures.exe xcorr_signed_future_abs_price_changes.exe xfit_mcsgarch_intraday.exe xfit_mcsgarch_intraday_batch.exe xfit_mcsgarch_on_intraday.exe xfit_mcsgarch_on_range_intraday.exe xfit_mcsegarch_on_intraday.exe xcompare_intraday_ewma_ohlc.exe xcompare_intraday_ewma_freq.exe $(OBJS) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJS6) $(OBJS7) $(OBJS8) $(OBJS9) $(OBJS10) $(OBJS11) $(OBJS12) $(OBJS13) $(OBJS14) $(OBJS15) $(OBJS16) $(OBJS17) $(OBJS18) $(OBJS19) $(OBJS20) $(OBJS21) $(OBJS22) $(OBJS23) $(OBJS24) $(OBJS25) $(OBJS26) $(OBJS27) $(OBJS28) $(OBJS29) $(OBJS30) $(OBJS31) $(OBJS32) $(OBJS33) $(OBJS34) $(OBJS35) $(OBJS36) $(OBJS37) $(OBJS38) $(OBJS39) $(OBJS40) $(OBJS41) $(OBJS42) $(OBJS43) $(OBJS44) $(OBJS45) $(OBJS46) $(OBJS47) $(OBJS48) $(OBJS49) $(OBJS50) $(OBJS51) gas.o linalg.o dcc.o stats.o sv.o distributions.o *.mod
