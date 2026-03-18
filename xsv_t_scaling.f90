! Simulate from SV(1)-lev-t and fit all four SV variants.
! Since data are SV-lev-t-generated, SV-lev-t should achieve the highest logL.
! Branched from xsv_lev_scaling.f90.

program xsv_t_scaling
    use kind_mod,    only: dp
    use sv_mod,      only: sv_simulate, sv_set_data, sv_set_types, sv_np, sv_obj, &
                           sv_sym_inv_transform, sv_lev_inv_transform, &
                           sv_t_inv_transform,   sv_lev_t_inv_transform, &
                           sv_transform, &
                           proc_sv, proc_sv_lev, dist_normal, dist_t, &
                           n_proc, n_dist, model_names
    use bfgs_module, only: bfgs_minimize
    implicit none

    ! ---- tuneable parameters ----
    integer,  parameter :: n_sizes  = 2
    integer,  parameter :: n_start  = 250
    integer,  parameter :: n_mult   = 6
    integer,  parameter :: seed_val = 42
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol     = 1.0e-7_dp

    ! ---- true SV(1)-lev-t parameters ----
    real(dp), parameter :: true_mu        =  0.0_dp
    real(dp), parameter :: true_phi       =  0.97_dp
    real(dp), parameter :: true_sigma_eta =  0.15_dp
    real(dp), parameter :: true_rho       = -0.6_dp
    real(dp), parameter :: true_nu        =  8.0_dp   ! t degrees of freedom

    integer, parameter :: nmod = n_proc * n_dist   ! 4 models

    real(dp), allocatable :: y(:), p(:), p0(:)
    real(dp) :: mu, phi, sigma_eta, rho, nu, f_opt
    logical  :: converged
    integer  :: np, n_iter, i, n, iproc, idist, imod

    integer  :: sizes(n_sizes)
    real(dp) :: logls(nmod, n_sizes)
    character(len=10) :: mnames(nmod)   ! flat list of model names for printing

    character(len=*), parameter :: hdr_fmt  = "(I9, 4(A12))"
    character(len=*), parameter :: row_fmt  = "(I9, 4(F12.1))"
    character(len=*), parameter :: delt_fmt = "(9X, 4(F12.1))"

    ! build flat name list in imod order
    do iproc = proc_sv, proc_sv_lev
        do idist = dist_normal, dist_t
            mnames((iproc-1)*n_dist + idist) = model_names(iproc, idist)
        end do
    end do

    n = n_start
    do i = 1, n_sizes
        sizes(i) = n
        n = n * n_mult
    end do

    do i = 1, n_sizes
        n = sizes(i)
        allocate(y(n))
        call sv_simulate(true_mu, true_phi, true_sigma_eta, true_rho, dist_t, true_nu, n, seed_val, y)
        call sv_set_data(y, n)

        do iproc = proc_sv, proc_sv_lev
            do idist = dist_normal, dist_t

                imod = (iproc-1)*n_dist + idist
                call sv_set_types(iproc, idist)
                np = sv_np()
                allocate(p(np), p0(np))

                select case (iproc*10 + idist)
                case (proc_sv*10     + dist_normal)
                    call sv_sym_inv_transform(0.0_dp, 0.90_dp, 0.30_dp, p0)
                case (proc_sv*10     + dist_t)
                    call sv_t_inv_transform(0.0_dp, 0.90_dp, 0.30_dp, 8.0_dp, p0)
                case (proc_sv_lev*10 + dist_normal)
                    call sv_lev_inv_transform(0.0_dp, 0.90_dp, 0.30_dp, -0.3_dp, p0)
                case (proc_sv_lev*10 + dist_t)
                    call sv_lev_t_inv_transform(0.0_dp, 0.90_dp, 0.30_dp, -0.3_dp, 8.0_dp, p0)
                end select

                p = p0
                call bfgs_minimize(sv_obj, p, np, max_iter, gtol, f_opt, n_iter, converged)
                call sv_transform(p, mu, phi, sigma_eta, rho, nu)
                logls(imod, i) = -n * f_opt

                deallocate(p, p0)

            end do
        end do

        deallocate(y)
    end do

    print '(A)', ""
    print '(A)', " SV(1)-lev-t data: all four SV models"
    print '(A,I0,A,I0,A)', " Sizes: ", n_start, " x ", n_mult, "^(0,1,2,...)"
    print '(A)', " True: mu=0.00  phi=0.97  sigma_eta=0.15  rho=-0.60  nu=8.0"
    print '(A)', ""

    ! header: n, then one column per model
    write(*, hdr_fmt) 0, (trim(mnames(imod)), imod=1, nmod)
    print '(A)', repeat("-", 9 + 4*12)

    do i = 1, n_sizes
        write(*, row_fmt) sizes(i), (logls(imod, i), imod=1, nmod)
    end do

    print '(A)', ""
    print '(A)', " Delta logL relative to SV-N:"
    write(*, hdr_fmt) 0, (trim(mnames(imod)), imod=1, nmod)
    print '(A)', repeat("-", 9 + 4*12)
    do i = 1, n_sizes
        write(*, delt_fmt) (logls(imod, i) - logls(1, i), imod=1, nmod)
    end do
    print '(A)', ""

end program xsv_t_scaling
