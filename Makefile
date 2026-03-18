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
SHARED  = kind.o math_const_mod.o garch_module.o bfgs_module.o
OBJS    = $(SHARED) garch_main.o
OBJS2   = $(SHARED) garch_scaling.o
OBJS3   = $(SHARED) special.o distributions.o random.o garch_t_module.o xgarch_t.o
OBJS4   = $(SHARED) special.o distributions.o random.o garch_t_module.o garch_sech_module.o garch_ged_module.o garch_laplace_module.o garch_logistic_module.o garch_nig_module.o xgarch_dist.o
OBJS5   = $(SHARED) special.o distributions.o random.o garch_t_module.o garch_sech_module.o garch_ged_module.o garch_laplace_module.o garch_logistic_module.o garch_nig_module.o garch_choose_mod.o xgarch_choose_dist.o
OBJS6   = $(SHARED) nagarch_module.o xnagarch.o
OBJS7   = $(SHARED) gjr_module.o xgjr.o
OBJS8   = $(SHARED) egarch_module.o xegarch.o
OBJS9   = $(SHARED) special.o distributions.o random.o garch_t_module.o nagarch_module.o gjr_module.o egarch_module.o \
          garch_flex_mod.o xgarch_flex.o
OBJS10  = kind.o csv_mod.o xread_csv.o
OBJS11  = kind.o math_const_mod.o special.o distributions.o garch_module.o bfgs_module.o nagarch_module.o \
          gjr_module.o egarch_module.o garch_flex_mod.o csv_mod.o rank_mod.o xfit_spy.o
OBJS12  = kind.o math_const_mod.o special.o distributions.o garch_module.o bfgs_module.o nagarch_module.o \
          gjr_module.o egarch_module.o garch_flex_mod.o csv_mod.o stats_mod.o rank_mod.o xfit_garch_returns.o
OBJS13  = kind.o math_const_mod.o special.o distributions.o garch_module.o bfgs_module.o nagarch_module.o \
          gjr_module.o egarch_module.o garch_flex_mod.o csv_mod.o stats_mod.o rank_mod.o xfit_garch_dist_returns.o
OBJS14  = kind.o math_const_mod.o special.o distributions.o random.o bfgs_module.o gas_mod.o csv_mod.o stats_mod.o rank_mod.o xfit_gas_returns.o
OBJS15  = kind.o math_const_mod.o special.o garch_module.o bfgs_module.o nagarch_module.o \
          gjr_module.o egarch_module.o garch_flex_mod.o distributions.o random.o gas_mod.o csv_mod.o stats_mod.o rank_mod.o xgarch_gas.o
OBJS16  = kind.o math_const_mod.o special.o distributions.o random.o bfgs_module.o gas_mod.o xgas_scaling.o
OBJS17  = kind.o math_const_mod.o special.o distributions.o random.o bfgs_module.o gas_mod.o \
          garch_module.o nagarch_module.o gjr_module.o egarch_module.o garch_flex_mod.o \
          xgas_garch_scaling.o
OBJS18  = kind.o math_const_mod.o special.o distributions.o random.o bfgs_module.o sv_mod.o xsv_scaling.o
OBJS19  = kind.o math_const_mod.o special.o distributions.o random.o bfgs_module.o sv_mod.o xsv_lev_scaling.o
OBJS20  = kind.o math_const_mod.o special.o distributions.o random.o bfgs_module.o sv_mod.o csv_mod.o stats_mod.o rank_mod.o xfit_sv_returns.o
OBJS21  = kind.o math_const_mod.o special.o distributions.o random.o bfgs_module.o sv_mod.o xsv_t_scaling.o
OBJS22  = kind.o math_const_mod.o special.o garch_module.o bfgs_module.o nagarch_module.o \
          gjr_module.o egarch_module.o garch_flex_mod.o fgarch_module.o distributions.o random.o sv_mod.o \
          csv_mod.o stats_mod.o rank_mod.o xfit_sv_garch_returns.o
OBJS23  = kind.o math_const_mod.o special.o distributions.o garch_module.o bfgs_module.o nagarch_module.o \
          gjr_module.o egarch_module.o garch_flex_mod.o csv_mod.o stats_mod.o rank_mod.o xnagarch_mix.o
OBJS24  = kind.o math_const_mod.o special.o distributions.o garch_module.o bfgs_module.o nagarch_module.o \
          gjr_module.o egarch_module.o garch_flex_mod.o csv_mod.o stats_mod.o rank_mod.o xnagarch_mix_t.o
OBJS25  = kind.o math_const_mod.o special.o distributions.o garch_module.o bfgs_module.o nagarch_module.o \
          gjr_module.o egarch_module.o garch_flex_mod.o csv_mod.o stats_mod.o rank_mod.o xstgarch.o
OBJS26  = kind.o math_const_mod.o garch_module.o bfgs_module.o csv_mod.o stats_mod.o xarch_ew.o
OBJS27  = kind.o math_const_mod.o garch_module.o bfgs_module.o csv_mod.o stats_mod.o xarch_lw.o
OBJS28  = kind.o math_const_mod.o garch_module.o bfgs_module.o csv_mod.o stats_mod.o xarch.o
OBJS29  = kind.o math_const_mod.o garch_module.o bfgs_module.o csv_mod.o stats_mod.o \
          linalg_mod.o dcc_mod.o xdcc.o
OBJS30  = kind.o math_const_mod.o garch_module.o nagarch_module.o bfgs_module.o \
          csv_mod.o stats_mod.o linalg_mod.o dcc_mod.o xadcc.o
OBJS31  = kind.o math_const_mod.o garch_module.o nagarch_module.o bfgs_module.o \
          csv_mod.o stats_mod.o linalg_mod.o dcc_mod.o xadcc_t.o
OBJS32  = kind.o math_const_mod.o garch_module.o bfgs_module.o \
          csv_mod.o stats_mod.o linalg_mod.o dcc_mod.o xdcc_vt.o
OBJS33  = kind.o math_const_mod.o garch_module.o nagarch_module.o bfgs_module.o \
          csv_mod.o stats_mod.o linalg_mod.o dcc_mod.o xadcc_vt.o
OBJS34  = kind.o math_const_mod.o garch_module.o nagarch_module.o bfgs_module.o \
          csv_mod.o stats_mod.o linalg_mod.o dcc_mod.o xadcc_t_vt.o
OBJS35  = kind.o math_const_mod.o special.o nagarch_module.o bfgs_module.o \
          csv_mod.o stats_mod.o distributions.o xdist.o

.PHONY: all run run_scaling run_t run_dist run_choose_dist run_gas_scaling run_gas_garch_scaling run_sv_scaling run_sv_lev_scaling run_fit_sv_returns run_sv_t_scaling run_fit_sv_garch_returns run_nagarch_mix run_nagarch_mix_t run_stgarch run_arch_ew run_arch_lw run_arch run_dcc run_adcc run_adcc_t run_dcc_vt run_adcc_vt run_adcc_t_vt run_xdist clean

all: $(EXE) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11) $(EXE12) $(EXE13) $(EXE14) $(EXE15) $(EXE16) $(EXE17) $(EXE18) $(EXE19) $(EXE20) $(EXE21) $(EXE22) $(EXE23) $(EXE24) $(EXE25) $(EXE26) $(EXE27) $(EXE28) $(EXE29) $(EXE30) $(EXE31) $(EXE32) $(EXE33) $(EXE34) $(EXE35)

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

kind.o: kind.f90
	$(FC) $(FFLAGS) -c $<

math_const_mod.o: math_const_mod.f90 kind.o
	$(FC) $(FFLAGS) -c $<

garch_module.o: garch_module.f90 kind.o math_const_mod.o
	$(FC) $(FFLAGS) -c $<

bfgs_module.o: bfgs_module.f90 kind.o
	$(FC) $(FFLAGS) -c $<

garch_main.o: garch_main.f90 kind.o garch_module.o bfgs_module.o
	$(FC) $(FFLAGS) -c $<

garch_scaling.o: garch_scaling.f90 kind.o garch_module.o bfgs_module.o
	$(FC) $(FFLAGS) -c $<

nagarch_module.o: nagarch_module.f90 kind.o math_const_mod.o
	$(FC) $(FFLAGS) -c $<

garch_t_module.o: garch_t_module.f90 kind.o math_const_mod.o garch_module.o special.o random.o
	$(FC) $(FFLAGS) -c $<

garch_sech_module.o: garch_sech_module.f90 kind.o math_const_mod.o garch_module.o
	$(FC) $(FFLAGS) -c $<

xgarch_t.o: xgarch_t.f90 kind.o garch_module.o garch_t_module.o bfgs_module.o
	$(FC) $(FFLAGS) -c $<

garch_ged_module.o: garch_ged_module.f90 kind.o garch_module.o special.o random.o distributions.o
	$(FC) $(FFLAGS) -c $<

garch_laplace_module.o: garch_laplace_module.f90 kind.o math_const_mod.o garch_module.o
	$(FC) $(FFLAGS) -c $<

garch_logistic_module.o: garch_logistic_module.f90 kind.o math_const_mod.o garch_module.o
	$(FC) $(FFLAGS) -c $<

garch_nig_module.o: garch_nig_module.f90 kind.o math_const_mod.o garch_module.o special.o random.o
	$(FC) $(FFLAGS) -c $<

xgarch_dist.o: xgarch_dist.f90 kind.o garch_module.o garch_t_module.o garch_sech_module.o garch_ged_module.o garch_laplace_module.o garch_logistic_module.o garch_nig_module.o bfgs_module.o
	$(FC) $(FFLAGS) -c $<

garch_choose_mod.o: garch_choose_mod.f90 kind.o
	$(FC) $(FFLAGS) -c $<

xgarch_choose_dist.o: xgarch_choose_dist.f90 kind.o garch_choose_mod.o garch_module.o garch_t_module.o garch_sech_module.o garch_ged_module.o garch_laplace_module.o garch_logistic_module.o garch_nig_module.o bfgs_module.o
	$(FC) $(FFLAGS) -c $<

xnagarch.o: xnagarch.f90 kind.o nagarch_module.o bfgs_module.o
	$(FC) $(FFLAGS) -c $<

gjr_module.o: gjr_module.f90 kind.o math_const_mod.o
	$(FC) $(FFLAGS) -c $<

xgjr.o: xgjr.f90 kind.o gjr_module.o bfgs_module.o
	$(FC) $(FFLAGS) -c $<

egarch_module.o: egarch_module.f90 kind.o math_const_mod.o
	$(FC) $(FFLAGS) -c $<

xegarch.o: xegarch.f90 kind.o egarch_module.o bfgs_module.o
	$(FC) $(FFLAGS) -c $<

garch_flex_mod.o: garch_flex_mod.f90 kind.o math_const_mod.o garch_module.o \
                  nagarch_module.o gjr_module.o egarch_module.o special.o distributions.o
	$(FC) $(FFLAGS) -c $<

xgarch_flex.o: xgarch_flex.f90 kind.o garch_flex_mod.o garch_t_module.o \
               gjr_module.o egarch_module.o bfgs_module.o
	$(FC) $(FFLAGS) -c $<

csv_mod.o: csv_mod.f90 kind.o
	$(FC) $(FFLAGS) -c $<

gas_mod.o: gas_mod.f90 kind.o math_const_mod.o random.o
	$(FC) $(FFLAGS) -c $<

stats_mod.o: stats_mod.f90 kind.o
	$(FC) $(FFLAGS) -c $<

rank_mod.o: rank_mod.f90 kind.o
	$(FC) $(FFLAGS) -c $<

xread_csv.o: xread_csv.f90 kind.o csv_mod.o
	$(FC) $(FFLAGS) -c $<

xfit_spy.o: xfit_spy.f90 kind.o garch_flex_mod.o garch_module.o nagarch_module.o \
            gjr_module.o egarch_module.o bfgs_module.o csv_mod.o rank_mod.o
	$(FC) $(FFLAGS) -c $<

xfit_garch_returns.o: xfit_garch_returns.f90 kind.o garch_flex_mod.o garch_module.o \
                      nagarch_module.o gjr_module.o egarch_module.o bfgs_module.o csv_mod.o rank_mod.o
	$(FC) $(FFLAGS) -c $<

xfit_garch_dist_returns.o: xfit_garch_dist_returns.f90 kind.o garch_flex_mod.o garch_module.o \
                           nagarch_module.o gjr_module.o egarch_module.o bfgs_module.o csv_mod.o rank_mod.o
	$(FC) $(FFLAGS) -c $<

xfit_gas_returns.o: xfit_gas_returns.f90 kind.o gas_mod.o bfgs_module.o csv_mod.o stats_mod.o rank_mod.o
	$(FC) $(FFLAGS) -c $<

xgarch_gas.o: xgarch_gas.f90 kind.o garch_flex_mod.o garch_module.o nagarch_module.o \
              gjr_module.o egarch_module.o gas_mod.o bfgs_module.o csv_mod.o stats_mod.o rank_mod.o
	$(FC) $(FFLAGS) -c $<

xgas_scaling.o: xgas_scaling.f90 kind.o gas_mod.o bfgs_module.o
	$(FC) $(FFLAGS) -c $<

xgas_garch_scaling.o: xgas_garch_scaling.f90 kind.o gas_mod.o garch_flex_mod.o \
                      garch_module.o bfgs_module.o
	$(FC) $(FFLAGS) -c $<

sv_mod.o: sv_mod.f90 kind.o math_const_mod.o special.o random.o
	$(FC) $(FFLAGS) -c $<

xsv_scaling.o: xsv_scaling.f90 kind.o sv_mod.o bfgs_module.o
	$(FC) $(FFLAGS) -c $<

xsv_lev_scaling.o: xsv_lev_scaling.f90 kind.o sv_mod.o bfgs_module.o
	$(FC) $(FFLAGS) -c $<

xfit_sv_returns.o: xfit_sv_returns.f90 kind.o sv_mod.o bfgs_module.o csv_mod.o stats_mod.o rank_mod.o
	$(FC) $(FFLAGS) -c $<

xsv_t_scaling.o: xsv_t_scaling.f90 kind.o sv_mod.o bfgs_module.o
	$(FC) $(FFLAGS) -c $<

fgarch_module.o: fgarch_module.f90 kind.o
	$(FC) $(FFLAGS) -c $<

xfit_sv_garch_returns.o: xfit_sv_garch_returns.f90 kind.o sv_mod.o garch_flex_mod.o \
                          fgarch_module.o garch_module.o nagarch_module.o gjr_module.o \
                          egarch_module.o bfgs_module.o csv_mod.o stats_mod.o rank_mod.o
	$(FC) $(FFLAGS) -c $<

xnagarch_mix.o: xnagarch_mix.f90 kind.o math_const_mod.o garch_flex_mod.o \
                nagarch_module.o bfgs_module.o csv_mod.o stats_mod.o rank_mod.o
	$(FC) $(FFLAGS) -c $<

xnagarch_mix_t.o: xnagarch_mix_t.f90 kind.o math_const_mod.o garch_flex_mod.o \
                  nagarch_module.o bfgs_module.o csv_mod.o stats_mod.o rank_mod.o
	$(FC) $(FFLAGS) -c $<

xstgarch.o: xstgarch.f90 kind.o math_const_mod.o garch_flex_mod.o nagarch_module.o \
            bfgs_module.o csv_mod.o stats_mod.o rank_mod.o
	$(FC) $(FFLAGS) -c $<

xarch_ew.o: xarch_ew.f90 kind.o math_const_mod.o garch_module.o bfgs_module.o \
            csv_mod.o stats_mod.o
	$(FC) $(FFLAGS) -c $<

xarch_lw.o: xarch_lw.f90 kind.o math_const_mod.o garch_module.o bfgs_module.o \
            csv_mod.o stats_mod.o
	$(FC) $(FFLAGS) -c $<

xarch.o: xarch.f90 kind.o math_const_mod.o garch_module.o bfgs_module.o \
         csv_mod.o stats_mod.o
	$(FC) $(FFLAGS) -c $<

linalg_mod.o: linalg_mod.f90 kind.o
	$(FC) $(FFLAGS) -c $<

dcc_mod.o: dcc_mod.f90 kind.o linalg_mod.o
	$(FC) $(FFLAGS) -c $<

xdcc.o: xdcc.f90 kind.o garch_module.o dcc_mod.o bfgs_module.o csv_mod.o stats_mod.o
	$(FC) $(FFLAGS) -c $<

xadcc.o: xadcc.f90 kind.o nagarch_module.o dcc_mod.o bfgs_module.o csv_mod.o stats_mod.o
	$(FC) $(FFLAGS) -c $<

xadcc_t.o: xadcc_t.f90 kind.o nagarch_module.o dcc_mod.o bfgs_module.o csv_mod.o stats_mod.o
	$(FC) $(FFLAGS) -c $<

xdcc_vt.o: xdcc_vt.f90 kind.o garch_module.o dcc_mod.o bfgs_module.o csv_mod.o stats_mod.o
	$(FC) $(FFLAGS) -c $<

xadcc_vt.o: xadcc_vt.f90 kind.o nagarch_module.o dcc_mod.o bfgs_module.o csv_mod.o stats_mod.o
	$(FC) $(FFLAGS) -c $<

xadcc_t_vt.o: xadcc_t_vt.f90 kind.o nagarch_module.o dcc_mod.o bfgs_module.o csv_mod.o stats_mod.o
	$(FC) $(FFLAGS) -c $<

random.o: random.f90 kind.o math_const_mod.o distributions.o
	$(FC) $(FFLAGS) -c $<

special.o: special.f90 kind.o math_const_mod.o
	$(FC) $(FFLAGS) -c $<

distributions.o: distributions.f90 kind.o math_const_mod.o bfgs_module.o special.o
	$(FC) $(FFLAGS) -c $<

xdist.o: xdist.f90 kind.o nagarch_module.o distributions.o bfgs_module.o csv_mod.o stats_mod.o
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

clean:
	rm -f random.o special.o $(OBJS) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJS6) $(OBJS7) $(OBJS8) $(OBJS9) $(OBJS10) $(OBJS11) $(OBJS12) $(OBJS13) $(OBJS14) $(OBJS15) $(OBJS16) $(OBJS17) $(OBJS18) $(OBJS19) $(OBJS20) $(OBJS21) $(OBJS22) $(OBJS23) $(OBJS24) $(OBJS25) $(OBJS26) $(OBJS27) $(OBJS28) $(OBJS29) $(OBJS30) $(OBJS31) $(OBJS32) $(OBJS33) $(OBJS34) $(OBJS35) $(EXE) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10) $(EXE11) $(EXE12) $(EXE13) $(EXE14) $(EXE15) $(EXE16) $(EXE17) $(EXE18) $(EXE19) $(EXE20) $(EXE21) $(EXE22) $(EXE23) $(EXE24) $(EXE25) $(EXE26) $(EXE27) $(EXE28) $(EXE29) $(EXE30) $(EXE31) $(EXE32) $(EXE33) $(EXE34) $(EXE35) gas_mod.o linalg_mod.o dcc_mod.o stats_mod.o sv_mod.o distributions.o *.mod
