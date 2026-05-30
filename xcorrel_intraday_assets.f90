! Compute cross-asset return covariance/correlation matrices at several
! intraday aggregation scales, plus realized volatilities by asset.

module correl_intraday_assets_mod
    use kind_mod, only: dp
    use date_mod, only: print_program_header
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_file, filter_intraday_session, &
                               infer_bar_interval_seconds
    use intraday_correlation_report_mod, only: print_frequency_correlation_report, print_correlation_anomaly_report
    use intraday_returns_mod, only: frequency_label
    use path_utils_mod, only: resolve_filename, basename_without_extension, files_with_extension_in_dir
    use strings_mod, only: uppercase
    implicit none
    private

    logical, parameter :: use_regular_session = .false.
    integer, parameter :: freq_seconds(*) = [integer :: &
        15*60, 30*60, 45*60, 60*60, 90*60, 120*60, 240*60, 480*60, 720*60, 1440*60]

    public :: run_correl_intraday_assets

contains

    ! Read assets and print covariance/correlation matrices by aggregation scale.
    subroutine run_correl_intraday_assets()
        character(len=512), allocatable :: filenames(:)
        character(len=16), allocatable :: asset_names(:)
        character(len=16), allocatable :: freq_labels(:)
        type(ohlcv_series_t), allocatable :: raw(:), bars(:)
        integer, allocatable :: nobs_by_freq(:)
        integer :: i, ifreq, nasset, source_seconds, nfreq
        real(dp), allocatable :: corr_by_freq(:, :, :)
        real(dp) :: t0, t1, read_sec, compute_sec

        call print_program_header("xcorrel_intraday_assets.f90")
        call input_filenames(filenames)
        nasset = size(filenames)
        allocate(asset_names(nasset), raw(nasset), bars(nasset))

        call cpu_time(t0)
        do i = 1, nasset
            call set_asset_name(filenames(i), asset_names(i))
            call read_intraday_prices_file(filenames(i), raw(i))
            if (use_regular_session) then
                call filter_intraday_session(raw(i), bars(i))
            else
                bars(i) = raw(i)
            end if
        end do
        call cpu_time(t1)
        read_sec = t1 - t0

        call cpu_time(t0)
        print '(A)', "Cross-asset intraday return covariance/correlation by time scale"
        print '(A,I0)', "Assets: ", nasset
        if (use_regular_session) then
            print '(A)', "Session: regular"
        else
            print '(A)', "Session: all-hours"
        end if
        print '(A)', "Files:"
        do i = 1, nasset
            print '(2X,A16,1X,A)', asset_names(i), trim(filenames(i))
        end do
        print '(A)', ""
        allocate(corr_by_freq(nasset, nasset, size(freq_seconds)))
        allocate(nobs_by_freq(size(freq_seconds)), freq_labels(size(freq_seconds)))
        corr_by_freq = 0.0_dp
        nobs_by_freq = 0
        freq_labels = ""
        nfreq = 0
        do ifreq = 1, size(freq_seconds)
            source_seconds = infer_bar_interval_seconds(bars(1))
            if (freq_seconds(ifreq) < source_seconds) cycle
            if (mod(freq_seconds(ifreq), source_seconds) /= 0) cycle
            nfreq = nfreq + 1
            freq_labels(nfreq) = frequency_label(freq_seconds(ifreq))
            call print_frequency_correlation_report(freq_seconds(ifreq), source_seconds, asset_names, bars, &
                                                    use_regular_session, corr_by_freq(:, :, nfreq), &
                                                    nobs_by_freq(nfreq))
        end do
        if (nfreq >= 2) then
            call print_correlation_anomaly_report(asset_names, freq_labels(1:nfreq), &
                                                  corr_by_freq(:, :, 1:nfreq), nobs_by_freq(1:nfreq))
        end if
        call cpu_time(t1)
        compute_sec = t1 - t0
        print '(/,A,F10.3)', "Read seconds:    ", read_sec
        print '(A,F10.3)', "Compute seconds: ", compute_sec
    end subroutine run_correl_intraday_assets

    ! Use command-line files if supplied; otherwise use ES/JY/TY continuous futures.
    subroutine input_filenames(filenames)
        character(len=512), allocatable, intent(out) :: filenames(:)
        character(len=512) :: arg, data_dir
        character(len=512), allocatable :: raw_names(:)
        integer :: nargs, i, nfile

        nargs = command_argument_count()
        data_dir = ""
        if (nargs > 0) then
            allocate(raw_names(nargs))
            nfile = 0
            do i = 1, nargs
                call get_command_argument(i, arg)
                if (index(arg, "dir=") == 1) then
                    data_dir = arg(5:)
                else
                    nfile = nfile + 1
                    raw_names(nfile) = arg
                end if
            end do
            if (nfile < 1) then
                if (len_trim(data_dir) < 1) error stop "input_filenames: no input files"
                call files_from_directory(data_dir, filenames)
            else
                allocate(filenames(nfile))
                do i = 1, nfile
                    filenames(i) = resolve_filename(raw_names(i), data_dir)
                end do
            end if
            deallocate(raw_names)
        else
            allocate(filenames(3))
            filenames = [character(len=512) :: &
                "c:\python\intraday_prices\continuous\ES.csv", &
                "c:\python\intraday_prices\continuous\JY.csv", &
                "c:\python\intraday_prices\continuous\TY.csv"]
        end if
    end subroutine input_filenames

    ! Use all BIN files in dir if present; otherwise use price CSV files.
    subroutine files_from_directory(data_dir, filenames)
        character(len=*), intent(in) :: data_dir
        character(len=512), allocatable, intent(out) :: filenames(:)
        character(len=512), allocatable :: files(:)

        call files_with_extension_in_dir(data_dir, ".bin", files)
        if (size(files) > 0) then
            allocate(filenames(size(files)))
            filenames = files
            deallocate(files)
            return
        end if
        call files_with_extension_in_dir(data_dir, ".csv", files)
        if (size(files) < 1) error stop "files_from_directory: no .bin or .csv files found"
        call keep_price_files(files, filenames)
        deallocate(files)
        if (size(filenames) < 1) error stop "files_from_directory: no price files found"
    end subroutine files_from_directory

    ! Skip known summary/report files when expanding a directory.
    subroutine keep_price_files(files_in, files_out)
        character(len=512), intent(in) :: files_in(:)
        character(len=512), allocatable, intent(out) :: files_out(:)
        character(len=512), allocatable :: tmp(:)
        character(len=512) :: base
        integer :: i, n

        allocate(tmp(size(files_in)))
        n = 0
        do i = 1, size(files_in)
            base = uppercase(basename_without_extension(files_in(i)))
            if (index(base, "SUMMARY") > 0) cycle
            n = n + 1
            tmp(n) = files_in(i)
        end do
        allocate(files_out(n))
        if (n > 0) files_out = tmp(1:n)
        deallocate(tmp)
    end subroutine keep_price_files

    ! Store a shortened basename label for a file.
    subroutine set_asset_name(filename, label)
        character(len=*), intent(in) :: filename
        character(len=*), intent(out) :: label
        character(len=512) :: base

        base = basename_without_extension(filename)
        label = ""
        label(1:min(len(label), len_trim(base))) = base(1:min(len(label), len_trim(base)))
    end subroutine set_asset_name

end module correl_intraday_assets_mod

program xcorrel_intraday_assets
    use correl_intraday_assets_mod, only: run_correl_intraday_assets
    implicit none

    call run_correl_intraday_assets()
end program xcorrel_intraday_assets
