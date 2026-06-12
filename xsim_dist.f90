! Simulate iid standardized data from selected distributions for xfit_dist.f90.
!
! The output CSV has an integer row index in the first column and one simulated
! series per distribution in the remaining columns.

program xsim_dist
    use date_mod, only: print_program_header
    use kind_mod, only: dp
    use distributions_mod, only: dist_id_from_name
    use random_mod, only: random_dist_std
    implicit none

    integer, parameter :: default_nobs = 10000
    integer, parameter :: default_seed = 123456
    character(len=*), parameter :: default_output_file = "dist_input.csv"
    character(len=16), parameter :: dist_list(*) = [character(len=16) :: &
        "NORMAL", "T", "GED", "LOGISTIC", "LAPLACE", "SECH", "NIG_SYM", "NIG", "FS_SKEWT"]
    character(len=16), parameter :: col_names(*) = [character(len=16) :: &
        "NORMAL", "T_6", "GED_1P5", "LOGISTIC", "LAPLACE", "SECH", "NIG_SYM_3", "NIG_A3_RHO0P5", &
        "FS_SKEWT_N6_X1P5"]
    real(dp), parameter :: shape(*) = [ &
        0.0_dp, 6.0_dp, 1.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp, 3.0_dp, 6.0_dp]
    real(dp), parameter :: shape2(*) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 1.5_dp]

    character(len=256) :: output_file, arg
    integer :: nobs, seed_val, unit, i, j, dist_id, ios
    integer :: seed_size
    integer, allocatable :: seed(:)
    real(dp) :: t_start, t_end

    call print_program_header("xsim_dist.f90")
    call cpu_time(t_start)

    nobs = default_nobs
    seed_val = default_seed
    output_file = default_output_file

    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *, iostat=ios) nobs
        if (ios /= 0 .or. nobs < 1) then
            print '(A)', "Usage: xsim_dist.exe [nobs>=1] [output_csv] [seed]"
            error stop "xsim_dist: invalid nobs"
        end if
    end if
    if (command_argument_count() >= 2) call get_command_argument(2, output_file)
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read(arg, *, iostat=ios) seed_val
        if (ios /= 0) then
            print '(A)', "Usage: xsim_dist.exe [nobs>=1] [output_csv] [seed]"
            error stop "xsim_dist: invalid seed"
        end if
    end if

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = seed_val
    call random_seed(put=seed)
    deallocate(seed)

    open(newunit=unit, file=output_file, status="replace", action="write")
    write(unit, '(A)', advance="no") "index"
    do j = 1, size(col_names)
        write(unit, '(A,A)', advance="no") ",", trim(col_names(j))
    end do
    write(unit, *)

    do i = 1, nobs
        write(unit, '(I0)', advance="no") i
        do j = 1, size(dist_list)
            dist_id = dist_id_from_name(dist_list(j))
            if (dist_id == 0) error stop "xsim_dist: unsupported distribution"
            write(unit, '(A,ES24.16)', advance="no") ",", random_dist_std(dist_id, shape(j), shape2(j))
        end do
        write(unit, *)
    end do
    close(unit)

    print '(A,A)', "Output file: ", trim(output_file)
    print '(A,I0,A,I0,A,I0)', "Rows: ", nobs, "  distributions: ", size(dist_list), "  seed: ", seed_val
    print '(A)', "Distribution columns: NORMAL T_6 GED_1P5 LOGISTIC LAPLACE SECH NIG_SYM_3 NIG_A3_RHO0P5 FS_SKEWT_N6_X1P5"

    call cpu_time(t_end)
    print '(/,A,F0.3,A)', "elapsed time: ", t_end - t_start, " s"
end program xsim_dist
