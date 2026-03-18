! Configuration type, derived types, and print subroutines for xgarch_choose_dist.

module garch_choose_mod
    use kind_mod, only: dp
    implicit none
    private

    ! ---- simulation distribution codes ----
    integer, parameter, public :: dist_normal   = 0
    integer, parameter, public :: dist_t        = 1
    integer, parameter, public :: dist_sech     = 2
    integer, parameter, public :: dist_ged      = 3
    integer, parameter, public :: dist_laplace  = 4
    integer, parameter, public :: dist_logistic = 5
    integer, parameter, public :: dist_nig      = 6

    ! ---- index constants for the fitted-distribution array ----
    integer, parameter, public :: i_normal   = 1
    integer, parameter, public :: i_t        = 2
    integer, parameter, public :: i_sech     = 3
    integer, parameter, public :: i_ged      = 4
    integer, parameter, public :: i_laplace  = 5
    integer, parameter, public :: i_logistic = 6
    integer, parameter, public :: i_nig      = 7
    integer, parameter, public :: n_fit      = 7

    character(len=9), parameter, public :: fit_names(n_fit) = &
        ["Normal   ", "Student-t", "Sech     ", "GED      ", "Laplace  ", "Logistic ", "NIG      "]

    ! number of parameters per fitted distribution
    integer, parameter, public :: np_fit(n_fit) = [3, 4, 3, 4, 3, 3, 4]

    ! has_extra(i): .true. if distribution i has a shape parameter
    logical, parameter, public :: has_extra(n_fit) = &
        [.false., .true., .false., .true., .false., .false., .true.]

    ! ---- user-settable configuration ----
    ! Declare a variable of this type in the main program and override any fields before use.
    type, public :: garch_config_t
        integer  :: nobs      = 1000  ! # of observations in each simulation
        integer  :: nsim       = 10   ! # of simulations
        integer  :: seed_val   = 42
        integer  :: max_iter   = 500
        real(dp) :: gtol       = 1.0e-7_dp
        real(dp) :: true_omega = 0.05_dp
        real(dp) :: true_alpha = 0.10_dp
        real(dp) :: true_beta  = 0.85_dp
        logical  :: fit_on(n_fit) = .true.  ! set .false. to skip fitting a distribution
    end type garch_config_t

    ! ---- simulation cases ----
    integer,  parameter, public :: sim_type(*) = [dist_sech,    dist_normal,   &
                                                   dist_t,       dist_t,        &
                                                   dist_ged,     dist_ged,      &
                                                   dist_laplace, dist_logistic, &
                                                   dist_nig,     dist_nig       ]
    real(dp), parameter, public :: sim_par(*)  = [0.0_dp, 0.0_dp,  &
                                                   5.0_dp, 15.0_dp, &
                                                   1.0_dp, 1.5_dp,  &
                                                   0.0_dp, 0.0_dp,  &
                                                   2.0_dp, 5.0_dp   ]

    ! ---- per-distribution fit results (single replication) ----
    type, public :: fit_result_t
        logical  :: fitted     = .false.
        real(dp) :: omega      = 0.0_dp
        real(dp) :: alpha      = 0.0_dp
        real(dp) :: beta       = 0.0_dp
        real(dp) :: extra_par  = 0.0_dp
        real(dp) :: uncond_vol = 0.0_dp
        real(dp) :: nll_n      = 0.0_dp
        real(dp) :: grad_norm  = 0.0_dp
        real(dp) :: time_s     = 0.0_dp
        integer  :: iterations = 0
        integer  :: rank_ll    = 0
        integer  :: rank_aic   = 0
        integer  :: rank_bic   = 0
        logical  :: converged  = .false.
    end type fit_result_t

    ! ---- accumulators over nsim replications ----
    type, public :: acc_t
        real(dp) :: omega(n_fit)      = 0.0_dp
        real(dp) :: alpha(n_fit)      = 0.0_dp
        real(dp) :: beta(n_fit)       = 0.0_dp
        real(dp) :: extra_par(n_fit)  = 0.0_dp
        real(dp) :: uncond_vol(n_fit) = 0.0_dp
        real(dp) :: nll_n(n_fit)      = 0.0_dp
        real(dp) :: time_s(n_fit)     = 0.0_dp
        integer  :: iterations(n_fit) = 0
        integer  :: rank1_ll(n_fit)   = 0
        integer  :: rank1_aic(n_fit)  = 0
        integer  :: rank1_bic(n_fit)  = 0
        integer  :: n_conv(n_fit)     = 0
    end type acc_t

    public :: print_results, print_summary, print_selection_matrices

contains

    subroutine print_results(res, cfg, sim_type_val, par, y_std, i_rep)
        type(fit_result_t),  intent(in) :: res(n_fit)
        type(garch_config_t),intent(in) :: cfg
        integer,             intent(in) :: sim_type_val, i_rep
        real(dp),            intent(in) :: par, y_std

        character(len=8) :: extra_str
        real(dp) :: true_uncvol, aic_i, bic_i, lnn
        integer  :: i, best_ll, best_aic, best_bic

        lnn = log(real(cfg%nobs, dp))

        print '(A)', ""
        select case (sim_type_val)
        case (dist_normal)
            print '(A,I0,A,I0,A)', " Simulated noise: Normal  [rep ", i_rep, "/", cfg%nsim, "]"
        case (dist_t)
            print '(A,F6.1,A,I0,A,I0,A)', " Simulated noise: Student-t, nu =", par, &
                "  [rep ", i_rep, "/", cfg%nsim, "]"
        case (dist_sech)
            print '(A,I0,A,I0,A)', " Simulated noise: Sech  [rep ", i_rep, "/", cfg%nsim, "]"
        case (dist_ged)
            print '(A,F6.1,A,I0,A,I0,A)', " Simulated noise: GED, nu =", par, &
                "  [rep ", i_rep, "/", cfg%nsim, "]"
        case (dist_laplace)
            print '(A,I0,A,I0,A)', " Simulated noise: Laplace  [rep ", i_rep, "/", cfg%nsim, "]"
        case (dist_logistic)
            print '(A,I0,A,I0,A)', " Simulated noise: Logistic  [rep ", i_rep, "/", cfg%nsim, "]"
        case (dist_nig)
            print '(A,F6.1,A,I0,A,I0,A)', " Simulated noise: NIG, alp =", par, &
                "  [rep ", i_rep, "/", cfg%nsim, "]"
        end select
        print '(A,F10.4)', " Simulated data std dev: ", y_std
        print '(A)', ""

        true_uncvol = sqrt(cfg%true_omega / (1.0_dp - cfg%true_alpha - cfg%true_beta))

        ! ---- parameter table ----
        print '(A)', "  Distribution     omega     alpha      beta     extra   uncond_vol"
        print '(A)', "  ------------   -------   -------   -------   -------  ----------"

        if (sim_type_val == dist_t .or. sim_type_val == dist_ged .or. sim_type_val == dist_nig) then
            write(extra_str, '(F8.4)') par
        else
            extra_str = "     ---"
        end if
        print '(2X,A12,3F10.4,2X,A8,F12.4)', &
            "True        ", cfg%true_omega, cfg%true_alpha, cfg%true_beta, extra_str, true_uncvol

        do i = 1, n_fit
            if (.not. res(i)%fitted) cycle
            if (has_extra(i)) then
                write(extra_str, '(F8.4)') res(i)%extra_par
            else
                extra_str = "     ---"
            end if
            print '(2X,A12,3F10.4,2X,A8,F12.4)', &
                fit_names(i)//"   ", res(i)%omega, res(i)%alpha, res(i)%beta, &
                extra_str, res(i)%uncond_vol
        end do

        print '(A)', ""

        ! ---- fit summary table ----
        print '(A)', "  Distribution   -loglik/n       AIC       BIC   rank  Arank  Brank  iter   grad_norm    time(s)"
        print '(A)', "  ------------   ---------   -------   -------  -----  -----  -----  ----  ----------  --------"
        do i = 1, n_fit
            if (.not. res(i)%fitted) cycle
            aic_i = 2*np_fit(i) + 2*cfg%nobs*res(i)%nll_n
            bic_i = np_fit(i)*lnn + 2*cfg%nobs*res(i)%nll_n
            print '(2X,A12,F12.4,2F10.1,3I7,I6,ES12.3,F10.3)', &
                fit_names(i)//"   ", res(i)%nll_n, aic_i, bic_i, &
                res(i)%rank_ll, res(i)%rank_aic, res(i)%rank_bic, &
                res(i)%iterations, res(i)%grad_norm, res(i)%time_s
        end do

        ! best models
        best_ll = 0;  best_aic = 0;  best_bic = 0
        do i = 1, n_fit
            if (.not. res(i)%fitted) cycle
            if (best_ll  == 0) then; best_ll  = i
            else if (res(i)%rank_ll  < res(best_ll)%rank_ll)  then; best_ll  = i; end if
            if (best_aic == 0) then; best_aic = i
            else if (res(i)%rank_aic < res(best_aic)%rank_aic) then; best_aic = i; end if
            if (best_bic == 0) then; best_bic = i
            else if (res(i)%rank_bic < res(best_bic)%rank_bic) then; best_bic = i; end if
        end do
        if (best_ll > 0) then
            print '(A)', ""
            print '(A,A)', "  Best by -loglik: ", trim(fit_names(best_ll))
            print '(A,A)', "  Best by AIC:     ", trim(fit_names(best_aic))
            print '(A,A)', "  Best by BIC:     ", trim(fit_names(best_bic))
        end if

        do i = 1, n_fit
            if (res(i)%fitted .and. .not. res(i)%converged) &
                print '(A,A,A)', "  Note: ", trim(fit_names(i)), " fit did not converge"
        end do

    end subroutine print_results

    subroutine print_summary(acc, cfg, sim_type_val, par)
        type(acc_t),         intent(in) :: acc
        type(garch_config_t),intent(in) :: cfg
        integer,             intent(in) :: sim_type_val
        real(dp),            intent(in) :: par

        character(len=8) :: extra_str
        real(dp) :: true_uncvol, mean_nll, mean_aic, mean_bic, rn, lnn
        integer  :: i

        rn  = real(cfg%nsim, dp)
        lnn = log(real(cfg%nobs, dp))

        print '(A)', ""
        print '(A)', " ---- Summary over replications ----"
        select case (sim_type_val)
        case (dist_normal);   print '(A)', " Simulated noise: Normal"
        case (dist_t);        print '(A,F6.1)', " Simulated noise: Student-t, nu =", par
        case (dist_sech);     print '(A)', " Simulated noise: Sech"
        case (dist_ged);      print '(A,F6.1)', " Simulated noise: GED, nu =", par
        case (dist_laplace);  print '(A)', " Simulated noise: Laplace"
        case (dist_logistic); print '(A)', " Simulated noise: Logistic"
        case (dist_nig);      print '(A,F6.1)', " Simulated noise: NIG, alp =", par
        end select
        print '(A)', ""

        true_uncvol = sqrt(cfg%true_omega / (1.0_dp - cfg%true_alpha - cfg%true_beta))

        ! ---- mean parameter table ----
        print '(A)', "  Distribution  mean_omega  mean_alpha   mean_beta  mean_extra  mean_uncvol"
        print '(A)', "  ------------  ----------  ----------  ----------  ----------  -----------"

        if (sim_type_val == dist_t .or. sim_type_val == dist_ged .or. sim_type_val == dist_nig) then
            write(extra_str, '(F8.4)') par
        else
            extra_str = "     ---"
        end if
        print '(2X,A12,3F12.4,2X,A8,F13.4)', &
            "True        ", cfg%true_omega, cfg%true_alpha, cfg%true_beta, extra_str, true_uncvol

        do i = 1, n_fit
            if (.not. cfg%fit_on(i)) cycle
            if (has_extra(i)) then
                write(extra_str, '(F8.4)') acc%extra_par(i) / rn
            else
                extra_str = "     ---"
            end if
            print '(2X,A12,3F12.4,2X,A8,F13.4)', &
                fit_names(i)//"   ", acc%omega(i)/rn, acc%alpha(i)/rn, acc%beta(i)/rn, &
                extra_str, acc%uncond_vol(i)/rn
        end do

        print '(A)', ""

        ! ---- mean fit summary table ----
        print '(A)', &
            "  Distribution  mean_nll/n   mean_AIC   mean_BIC  %best_ll  %best_AIC  %best_BIC" // &
            "  mean_iter  total_time  %conv"
        print '(A)', &
            "  ------------  ----------   --------   --------  --------  ---------  ---------" // &
            "  ---------  ----------  -----"
        do i = 1, n_fit
            if (.not. cfg%fit_on(i)) cycle
            mean_nll = acc%nll_n(i) / rn
            mean_aic = 2*np_fit(i) + 2*cfg%nobs * mean_nll
            mean_bic = np_fit(i)*lnn + 2*cfg%nobs * mean_nll
            print '(2X,A12,F12.4,2F11.1,3F11.1,F11.1,F12.3,F7.1)', &
                fit_names(i)//"   ", mean_nll, mean_aic, mean_bic, &
                100.0_dp * acc%rank1_ll(i)  / rn, &
                100.0_dp * acc%rank1_aic(i) / rn, &
                100.0_dp * acc%rank1_bic(i) / rn, &
                real(acc%iterations(i), dp) / rn, &
                acc%time_s(i), &
                100.0_dp * acc%n_conv(i) / rn
        end do

    end subroutine print_summary

    subroutine print_selection_matrices(acc_all, n_cases, cfg, sim_types, sim_pars)
        ! Prints LL, AIC, and BIC selection-proportion matrices.
        ! Rows = simulated distributions; columns = fitted distributions.
        type(acc_t),          intent(in) :: acc_all(n_cases)   ! per-case accumulators
        integer,              intent(in) :: n_cases             ! number of simulation cases
        type(garch_config_t), intent(in) :: cfg                 ! run configuration
        integer,              intent(in) :: sim_types(n_cases)  ! dist_* code for each case
        real(dp),             intent(in) :: sim_pars(n_cases)   ! shape parameter for each case

        character(len=16) :: label
        real(dp) :: rn
        integer  :: ic, j

        rn = real(cfg%nsim, dp)

        print '(A)', ""
        print '(A)', " ---- Selection matrices (% of replications ranked #1) ----"

        print '(A)', ""
        print '(A)', " By -log-likelihood:"
        print '(16X,7A10)', (fit_names(j), j=1, n_fit)
        do ic = 1, n_cases
            call sim_label(sim_types(ic), sim_pars(ic), label)
            print '(A16,7F10.1)', label, (100.0_dp * acc_all(ic)%rank1_ll(j) / rn, j=1, n_fit)
        end do

        print '(A)', ""
        print '(A)', " By AIC:"
        print '(16X,7A10)', (fit_names(j), j=1, n_fit)
        do ic = 1, n_cases
            call sim_label(sim_types(ic), sim_pars(ic), label)
            print '(A16,7F10.1)', label, (100.0_dp * acc_all(ic)%rank1_aic(j) / rn, j=1, n_fit)
        end do

        print '(A)', ""
        print '(A)', " By BIC:"
        print '(16X,7A10)', (fit_names(j), j=1, n_fit)
        do ic = 1, n_cases
            call sim_label(sim_types(ic), sim_pars(ic), label)
            print '(A16,7F10.1)', label, (100.0_dp * acc_all(ic)%rank1_bic(j) / rn, j=1, n_fit)
        end do

    end subroutine print_selection_matrices

    subroutine sim_label(sim_type_val, par, label)
        ! Returns a short row label for a simulated distribution.
        integer,           intent(in)  :: sim_type_val  ! dist_* code
        real(dp),          intent(in)  :: par           ! shape parameter (if any)
        character(len=16), intent(out) :: label         ! output label string

        select case (sim_type_val)
        case (dist_normal);   label = "Normal"
        case (dist_t);        write(label, '(A,F5.1)') "t, nu=", par
        case (dist_sech);     label = "Sech"
        case (dist_ged);      write(label, '(A,F5.1)') "GED, nu=", par
        case (dist_laplace);  label = "Laplace"
        case (dist_logistic); label = "Logistic"
        case (dist_nig);      write(label, '(A,F5.1)') "NIG, alp=", par
        case default;         label = "?"
        end select

    end subroutine sim_label

end module garch_choose_mod
