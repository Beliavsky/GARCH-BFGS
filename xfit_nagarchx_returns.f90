! Fit NAGARCH-X(1,1) and compare to NAGARCH on daily CC returns.
!
!   NAGARCH-X: h_t = omega + alpha*(r_{t-1} - theta*sqrt(h_{t-1}))^2
!                          + delta*RV_{t-1} + beta*h_{t-1}
!
! RV is the daily sum of squared intraday log-returns (regular session only).
! Usage: xfit_nagarchx_returns [intraday_file] [ntest]

program xfit_nagarchx_returns
    use nagarchx_compare_mod, only: run_nagarchx_compare
    implicit none
    call run_nagarchx_compare()
end program xfit_nagarchx_returns
