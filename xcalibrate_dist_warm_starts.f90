! Calibrate warm-start shape parameters for iid distribution fits.
!
! For each no-shape source distribution, choose target shape parameters for
! t, GED, and NIG by minimizing KL(source || target). Since the source entropy
! is fixed, this minimizes the expected target negative log density.
! Also reports the corresponding two-parameter scale+shape optimum.

program xcalibrate_dist_warm_starts
    use kind_mod, only: dp
    use date_mod, only: print_program_header
    use distributions_mod, only: dist_normal, dist_t, dist_ged, dist_logistic, dist_laplace, dist_sech, dist_nig, &
                                 dist_names, nll_std_val, pdf_normal, pdf_logistic, pdf_laplace, pdf_sech
    implicit none

    integer, parameter :: nsrc = 4
    integer, parameter :: ntgt = 3
    integer, parameter :: source_ids(nsrc) = [dist_normal, dist_logistic, dist_laplace, dist_sech]
    integer, parameter :: target_ids(ntgt) = [dist_t, dist_ged, dist_nig]
    real(dp), parameter :: x_max = 40.0_dp
    integer, parameter :: n_grid = 20000
    real(dp) :: sigma, shape, obj
    integer :: isrc, itgt

    call print_program_header("xcalibrate_dist_warm_starts.f90")

    print '(A)', "KL-calibrated warm-start shapes"
    print '(A)', repeat("-", 58)
    print '(A12,1X,A12,1X,A12,1X,A14)', "source", "target", "shape", "cross_entropy"
    print '(A)', repeat("-", 58)
    do isrc = 1, nsrc
        do itgt = 1, ntgt
            call calibrate_pair(source_ids(isrc), target_ids(itgt), shape, obj)
            print '(A12,1X,A12,1X,F12.6,1X,F14.8)', &
                trim(dist_names(source_ids(isrc))), trim(dist_names(target_ids(itgt))), shape, obj
        end do
    end do
    print '(A)', repeat("-", 58)

    print '(/,A)', "KL-calibrated warm-start sigma and shape"
    print '(A)', repeat("-", 72)
    print '(A12,1X,A12,1X,A12,1X,A12,1X,A14)', "source", "target", "sigma", "shape", "cross_entropy"
    print '(A)', repeat("-", 72)
    do isrc = 1, nsrc
        do itgt = 1, ntgt
            call calibrate_pair_scale_shape(source_ids(isrc), target_ids(itgt), sigma, shape, obj)
            print '(A12,1X,A12,1X,F12.6,1X,F12.6,1X,F14.8)', &
                trim(dist_names(source_ids(isrc))), trim(dist_names(target_ids(itgt))), sigma, shape, obj
        end do
    end do
    print '(A)', repeat("-", 72)

    print '(/,A)', "Fortran rows by source [t, ged, nig]:"
    do isrc = 1, nsrc
        write(*,'(A12,A)', advance='no') trim(dist_names(source_ids(isrc))), " : ["
        do itgt = 1, ntgt
            call calibrate_pair(source_ids(isrc), target_ids(itgt), shape, obj)
            if (itgt > 1) write(*,'(A)', advance='no') ", "
            write(*,'(F10.6)', advance='no') shape
        end do
        write(*,'(A)') "]"
    end do

    print '(/,A)', "Fortran sigma rows by source [t, ged, nig]:"
    do isrc = 1, nsrc
        write(*,'(A12,A)', advance='no') trim(dist_names(source_ids(isrc))), " : ["
        do itgt = 1, ntgt
            call calibrate_pair_scale_shape(source_ids(isrc), target_ids(itgt), sigma, shape, obj)
            if (itgt > 1) write(*,'(A)', advance='no') ", "
            write(*,'(F10.6)', advance='no') sigma
        end do
        write(*,'(A)') "]"
    end do

    print '(/,A)', "Fortran shape rows by source [t, ged, nig]:"
    do isrc = 1, nsrc
        write(*,'(A12,A)', advance='no') trim(dist_names(source_ids(isrc))), " : ["
        do itgt = 1, ntgt
            call calibrate_pair_scale_shape(source_ids(isrc), target_ids(itgt), sigma, shape, obj)
            if (itgt > 1) write(*,'(A)', advance='no') ", "
            write(*,'(F10.6)', advance='no') shape
        end do
        write(*,'(A)') "]"
    end do

contains

    subroutine calibrate_pair(source_id, target_id, shape, obj)
        ! Minimize expected target NLL over the target shape parameter.
        integer, intent(in) :: source_id, target_id
        real(dp), intent(out) :: shape, obj
        real(dp) :: lo, hi

        select case (target_id)
        case (dist_t)
            lo = 2.05_dp
            hi = 99.0_dp
        case (dist_ged)
            lo = 0.5_dp
            hi = 10.0_dp
        case (dist_nig)
            lo = 0.2_dp
            hi = 19.5_dp
        case default
            error stop "calibrate_pair: unsupported target"
        end select
        call golden_min(source_id, target_id, lo, hi, shape, obj)
    end subroutine calibrate_pair

    subroutine calibrate_pair_scale_shape(source_id, target_id, sigma, shape, obj)
        ! Coordinate-minimize expected target NLL over scale and shape.
        integer, intent(in) :: source_id, target_id
        real(dp), intent(out) :: sigma, shape, obj
        real(dp) :: lo, hi, prev_obj, new_obj
        integer :: iter

        call calibrate_pair(source_id, target_id, shape, obj)
        sigma = 1.0_dp
        prev_obj = cross_entropy_scaled(source_id, target_id, sigma, shape)
        do iter = 1, 20
            call golden_min_sigma(source_id, target_id, shape, 0.2_dp, 3.0_dp, sigma, new_obj)
            select case (target_id)
            case (dist_t)
                lo = 2.05_dp
                hi = 99.0_dp
            case (dist_ged)
                lo = 0.5_dp
                hi = 10.0_dp
            case (dist_nig)
                lo = 0.2_dp
                hi = 19.5_dp
            case default
                error stop "calibrate_pair_scale_shape: unsupported target"
            end select
            call golden_min_shape_scaled(source_id, target_id, sigma, lo, hi, shape, new_obj)
            if (abs(prev_obj - new_obj) < 1.0e-10_dp) exit
            prev_obj = new_obj
        end do
        obj = cross_entropy_scaled(source_id, target_id, sigma, shape)
    end subroutine calibrate_pair_scale_shape

    subroutine golden_min(source_id, target_id, lo, hi, xmin, fmin)
        ! Golden-section minimization on a bounded interval.
        integer, intent(in) :: source_id, target_id
        real(dp), intent(in) :: lo, hi
        real(dp), intent(out) :: xmin, fmin
        real(dp), parameter :: gr = 0.6180339887498948482_dp
        real(dp) :: a, b, c, d, fc, fd
        integer :: iter

        a = lo
        b = hi
        c = b - gr*(b - a)
        d = a + gr*(b - a)
        fc = cross_entropy(source_id, target_id, c)
        fd = cross_entropy(source_id, target_id, d)
        do iter = 1, 80
            if (fc < fd) then
                b = d
                d = c
                fd = fc
                c = b - gr*(b - a)
                fc = cross_entropy(source_id, target_id, c)
            else
                a = c
                c = d
                fc = fd
                d = a + gr*(b - a)
                fd = cross_entropy(source_id, target_id, d)
            end if
        end do
        if (fc < fd) then
            xmin = c
            fmin = fc
        else
            xmin = d
            fmin = fd
        end if
    end subroutine golden_min

    subroutine golden_min_sigma(source_id, target_id, shape, lo, hi, xmin, fmin)
        ! Golden-section minimization over target scale.
        integer, intent(in) :: source_id, target_id
        real(dp), intent(in) :: shape, lo, hi
        real(dp), intent(out) :: xmin, fmin
        real(dp), parameter :: gr = 0.6180339887498948482_dp
        real(dp) :: a, b, c, d, fc, fd
        integer :: iter

        a = lo
        b = hi
        c = b - gr*(b - a)
        d = a + gr*(b - a)
        fc = cross_entropy_scaled(source_id, target_id, c, shape)
        fd = cross_entropy_scaled(source_id, target_id, d, shape)
        do iter = 1, 70
            if (fc < fd) then
                b = d
                d = c
                fd = fc
                c = b - gr*(b - a)
                fc = cross_entropy_scaled(source_id, target_id, c, shape)
            else
                a = c
                c = d
                fc = fd
                d = a + gr*(b - a)
                fd = cross_entropy_scaled(source_id, target_id, d, shape)
            end if
        end do
        if (fc < fd) then
            xmin = c
            fmin = fc
        else
            xmin = d
            fmin = fd
        end if
    end subroutine golden_min_sigma

    subroutine golden_min_shape_scaled(source_id, target_id, sigma, lo, hi, xmin, fmin)
        ! Golden-section minimization over target shape with scale fixed.
        integer, intent(in) :: source_id, target_id
        real(dp), intent(in) :: sigma, lo, hi
        real(dp), intent(out) :: xmin, fmin
        real(dp), parameter :: gr = 0.6180339887498948482_dp
        real(dp) :: a, b, c, d, fc, fd
        integer :: iter

        a = lo
        b = hi
        c = b - gr*(b - a)
        d = a + gr*(b - a)
        fc = cross_entropy_scaled(source_id, target_id, sigma, c)
        fd = cross_entropy_scaled(source_id, target_id, sigma, d)
        do iter = 1, 70
            if (fc < fd) then
                b = d
                d = c
                fd = fc
                c = b - gr*(b - a)
                fc = cross_entropy_scaled(source_id, target_id, sigma, c)
            else
                a = c
                c = d
                fc = fd
                d = a + gr*(b - a)
                fd = cross_entropy_scaled(source_id, target_id, sigma, d)
            end if
        end do
        if (fc < fd) then
            xmin = c
            fmin = fc
        else
            xmin = d
            fmin = fd
        end if
    end subroutine golden_min_shape_scaled

    function cross_entropy(source_id, target_id, shape) result(val)
        ! Simpson integration of E_source[-log target_density(X)].
        integer, intent(in) :: source_id, target_id
        real(dp), intent(in) :: shape
        real(dp) :: val, dx, x, w, dens
        integer :: i

        dx = 2.0_dp*x_max / real(n_grid, dp)
        val = 0.0_dp
        do i = 0, n_grid
            x = -x_max + dx*real(i, dp)
            if (i == 0 .or. i == n_grid) then
                w = 1.0_dp
            else if (mod(i, 2) == 0) then
                w = 2.0_dp
            else
                w = 4.0_dp
            end if
            dens = source_pdf(source_id, x)
            if (dens > 0.0_dp) val = val + w*dens*nll_std_val(x, target_id, shape)
        end do
        val = val*dx/3.0_dp
    end function cross_entropy

    function cross_entropy_scaled(source_id, target_id, sigma, shape) result(val)
        ! Simpson integration of E_source[-log target_density_sigma(X)].
        integer, intent(in) :: source_id, target_id
        real(dp), intent(in) :: sigma, shape
        real(dp) :: val, dx, x, z, w, dens, sig
        integer :: i

        sig = max(sigma, 1.0e-12_dp)
        dx = 2.0_dp*x_max / real(n_grid, dp)
        val = 0.0_dp
        do i = 0, n_grid
            x = -x_max + dx*real(i, dp)
            if (i == 0 .or. i == n_grid) then
                w = 1.0_dp
            else if (mod(i, 2) == 0) then
                w = 2.0_dp
            else
                w = 4.0_dp
            end if
            dens = source_pdf(source_id, x)
            z = x / sig
            if (dens > 0.0_dp) val = val + w*dens*(log(sig) + nll_std_val(z, target_id, shape))
        end do
        val = val*dx/3.0_dp
    end function cross_entropy_scaled

    function source_pdf(source_id, x) result(dens)
        ! PDF for a no-shape source distribution.
        integer, intent(in) :: source_id
        real(dp), intent(in) :: x
        real(dp) :: dens

        select case (source_id)
        case (dist_normal)
            dens = pdf_normal(x)
        case (dist_logistic)
            dens = pdf_logistic(x)
        case (dist_laplace)
            dens = pdf_laplace(x)
        case (dist_sech)
            dens = pdf_sech(x)
        case default
            dens = 0.0_dp
        end select
    end function source_pdf

end program xcalibrate_dist_warm_starts
