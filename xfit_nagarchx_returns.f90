! Fit NAGARCH-X(1,1) and compare to NAGARCH on daily CC returns.
!
!   NAGARCH-X: h_t = omega + alpha*(r_{t-1} - theta*sqrt(h_{t-1}))^2
!                          + delta*RV_{t-1} + beta*h_{t-1}
!
! RV is the daily sum of squared intraday log-returns (regular session only).
! Usage: xfit_nagarchx_returns [file1 [file2 ...]]
!   Files may be *.csv or *.bin intraday price files; glob patterns are expanded.
!   With no arguments the default file_pattern below is used.

program xfit_nagarchx_returns
    use kind_mod,        only: dp
    use date_mod,        only: print_program_header
    use input_files_mod, only: collect_input_filenames, MAX_PATH_LEN
    use nagarchx_compare_mod, only: run_nagarchx_compare
    implicit none
    real(dp) :: t_start, t_end
    character(len=MAX_PATH_LEN), allocatable :: filenames(:)
    character(len=*), parameter :: file_pattern = "c:\python\databento\data_1min\*.bin"

    call print_program_header("xfit_nagarchx_returns.f90")
    call cpu_time(t_start)
    call collect_input_filenames(filenames, &
        file_pattern=file_pattern, &
        default_filenames=[character(len=MAX_PATH_LEN) :: &
            "c:\python\databento\data_1min\spy_1min_databento.csv"])
    call run_nagarchx_compare(filenames)
    deallocate(filenames)
    call cpu_time(t_end)
    write(*,'(a,f10.3)') "Overall elapsed seconds:", t_end - t_start
end program xfit_nagarchx_returns
