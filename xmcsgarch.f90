! Simulate and fit a multiplicative component sGARCH intraday return series.
!
! This is a smoke-test driver for garch_mcsgarch_mod.  The model is:
!   var(r_t) = daily_var_t * diurnal_var_t * q_t
! with q_t following a GARCH(1,1) recursion on normalized intraday residuals.

module xmcsgarch_mod
    use kind_mod, only: dp
    use garch_mcsgarch_mod, only: mcsgarch_params_t, mcsgarch_fit_result_t, &
                                  mcsgarch_simulate, fit_mcsgarch
    implicit none
    private
    public :: run_xmcsgarch

contains

    subroutine run_xmcsgarch()
        integer, parameter :: ndays = 500
        integer, parameter :: nbins = 78
        integer, parameter :: n = ndays * nbins
        type(mcsgarch_params_t) :: true_params
        type(mcsgarch_fit_result_t) :: fit
        real(dp), allocatable :: y(:), daily_var(:), true_diurnal(:), fit_diurnal(:), q(:)
        integer, allocatable :: bin_id(:)
        integer :: day, bin, t
        real(dp) :: x, pi

        allocate(y(n), daily_var(n), true_diurnal(n), fit_diurnal(n), q(n), bin_id(n))
        true_params = mcsgarch_params_t(0.04_dp, 0.07_dp, 0.90_dp)
        pi = acos(-1.0_dp)

        do day = 1, ndays
            do bin = 1, nbins
                t = (day - 1)*nbins + bin
                x = real(bin - 1, dp) / real(nbins - 1, dp)
                bin_id(t) = bin
                daily_var(t) = 1.0e-4_dp * (1.0_dp + 0.20_dp*sin(2.0_dp*pi*real(day, dp)/63.0_dp))
                true_diurnal(t) = 0.60_dp + 1.80_dp*(x - 0.5_dp)**2
            end do
        end do
        true_diurnal = true_diurnal / (sum(true_diurnal(1:nbins)) / real(nbins, dp))

        call mcsgarch_simulate(true_params, daily_var, true_diurnal, 12345, y, q)
        call fit_mcsgarch(y, daily_var, bin_id, 120, 1.0e-5_dp, fit, fit_diurnal, q)

        print '(A)', "MCS-GARCH simulate -> fit sanity check"
        print '(A,I0,A,I0,A,I0)', "Observations: ", n, "  days: ", ndays, "  intraday bins/day: ", nbins
        print '(A)', "Daily variance is supplied exogenously; diurnal variance is estimated by intraday bin."
        print '(A)', "---------------------------------------------------------------"
        print '(A,3X,A,9X,A,8X,A,8X,A,7X,A,5X,A)', "Row", "omega", "alpha", "beta", "persist", "niter", "conv"
        print '(A)', "---------------------------------------------------------------"
        print '(A,1X,ES12.4,2(1X,F10.4),1X,F10.4,5X,A,5X,A)', &
              "TRUE", true_params%omega, true_params%alpha, true_params%beta, &
              true_params%alpha + true_params%beta, "-", "-"
        print '(A,2X,ES12.4,2(1X,F10.4),1X,F10.4,1X,I9,5X,L1)', &
              "FIT", fit%params%omega, fit%params%alpha, fit%params%beta, &
              fit%persist, fit%niter, fit%converged
        print '(A)', "---------------------------------------------------------------"
        print '(A,F12.4)', "logL: ", fit%loglik

        deallocate(y, daily_var, true_diurnal, fit_diurnal, q, bin_id)
    end subroutine run_xmcsgarch

end module xmcsgarch_mod

program xmcsgarch
    use xmcsgarch_mod, only: run_xmcsgarch
    implicit none
    call run_xmcsgarch()
end program xmcsgarch
