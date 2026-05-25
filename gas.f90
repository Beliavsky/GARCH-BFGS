! GAS(1,1) volatility models with Normal and Student-t innovations.
!
! The log-variance f_t = log(h_t) evolves as:
!   f_{t+1} = omega + (alpha + gamma*I(u_t < 0))*s_t + beta*f_t
!
! where s_t is the score of log p(y_t | f_t) w.r.t. f_t, and
! gamma = 0 for the symmetric model (proc_gas).
! beta is constrained to (-1,1) via a tanh transform.
!
! Score functions (unscaled; alpha absorbs the scale):
!   Normal: s_t = (u_t^2 - 1) / 2
!   t:      s_t = (nu*u_t^2 - (nu-2)) / (2*(nu-2 + u_t^2))   [bounded]
!
! Unconditional log-variance: E[f] = omega/(1-beta)  =>  h_unc = exp(omega/(1-beta))
!
! Parameter layout (unconstrained vector p):
!   proc_gas  + dist_normal: [omega, alpha, atanh(beta)]           np=3
!   proc_gas  + dist_t:      [omega, alpha, atanh(beta), p_nu]     np=4
!   proc_agas + dist_normal: [omega, alpha, gamma, atanh(beta)]    np=4
!   proc_agas + dist_t:      [omega, alpha, gamma, atanh(beta), p_nu] np=5

module gas_mod
    use kind_mod,       only: dp
    use math_const_mod, only: pi, log_sqrt_2pi
    use random_mod,     only: random_normal, random_t_std
    implicit none
    private
    public :: gas_set_data, gas_set_types, gas_np, gas_obj, gas_simulate
    public :: gas_sym_inv_transform, gas_asym_inv_transform, gas_transform
    public :: proc_gas, proc_agas, dist_normal, dist_t, n_proc, n_dist
    public :: proc_names, dist_names, has_shape

    integer, parameter :: proc_gas    = 1  ! symmetric GAS(1,1)
    integer, parameter :: proc_agas   = 2  ! asymmetric: larger score weight for u_t < 0
    integer, parameter :: dist_normal = 1
    integer, parameter :: dist_t      = 2
    integer, parameter :: n_proc      = proc_agas
    integer, parameter :: n_dist      = dist_t

    character(len=8), parameter :: proc_names(n_proc) = &
        ["GAS     ", "AGAS    "]
    character(len=8), parameter :: dist_names(n_dist) = &
        ["Normal  ", "t       "]
    logical, parameter :: has_shape(n_dist) = [.false., .true.]

    real(dp), allocatable, save :: y_mod(:)  ! observation data
    integer,               save :: n_mod     ! number of observations
    integer,               save :: iproc_mod ! current process type
    integer,               save :: idist_mod ! current distribution type

contains

    subroutine gas_set_data(y, n)
        ! Store observation data for use by gas_obj.
        real(dp), intent(in) :: y(:)  ! demeaned returns
        integer,  intent(in) :: n     ! number of observations
        if (allocated(y_mod)) deallocate(y_mod)
        allocate(y_mod(n))
        y_mod = y(1:n)
        n_mod = n
    end subroutine gas_set_data

    subroutine gas_set_types(iproc, idist)
        ! Set process and distribution type for the next gas_np / gas_obj call.
        integer, intent(in) :: iproc  ! proc_gas or proc_agas
        integer, intent(in) :: idist  ! dist_normal or dist_t
        iproc_mod = iproc
        idist_mod = idist
    end subroutine gas_set_types

    function gas_np() result(np)
        ! Number of unconstrained parameters for the current process/distribution.
        integer :: np
        np = 3                                       ! omega, alpha, atanh(beta)
        if (iproc_mod == proc_agas) np = np + 1     ! gamma
        if (idist_mod == dist_t)    np = np + 1     ! p_nu
    end function gas_np

    subroutine gas_sym_inv_transform(omega, alpha, beta, p)
        ! Map constrained symmetric GAS parameters to unconstrained vector p(1:3).
        real(dp), intent(in)  :: omega  ! log-variance intercept (unconstrained)
        real(dp), intent(in)  :: alpha  ! score coefficient (unconstrained)
        real(dp), intent(in)  :: beta   ! persistence, |beta| < 1
        real(dp), intent(out) :: p(3)   ! [omega, alpha, atanh(beta)]
        p(1) = omega
        p(2) = alpha
        p(3) = atanh(beta)
    end subroutine gas_sym_inv_transform

    subroutine gas_asym_inv_transform(omega, alpha, gamma, beta, p)
        ! Map constrained asymmetric GAS parameters to unconstrained vector p(1:4).
        real(dp), intent(in)  :: omega   ! log-variance intercept (unconstrained)
        real(dp), intent(in)  :: alpha   ! score coefficient (unconstrained)
        real(dp), intent(in)  :: gamma   ! asymmetry coefficient (unconstrained)
        real(dp), intent(in)  :: beta    ! persistence, |beta| < 1
        real(dp), intent(out) :: p(4)    ! [omega, alpha, gamma, atanh(beta)]
        p(1) = omega
        p(2) = alpha
        p(3) = gamma
        p(4) = atanh(beta)
    end subroutine gas_asym_inv_transform

    subroutine gas_transform(p, omega, alpha, gamma, beta)
        ! Map unconstrained vector p to constrained GAS parameters.
        real(dp), intent(in)  :: p(:)    ! unconstrained parameter vector
        real(dp), intent(out) :: omega   ! log-variance intercept
        real(dp), intent(out) :: alpha   ! score coefficient
        real(dp), intent(out) :: gamma   ! asymmetry coefficient (0 for proc_gas)
        real(dp), intent(out) :: beta    ! persistence in (-1,1)
        omega = p(1)
        alpha = p(2)
        if (iproc_mod == proc_agas) then
            gamma = p(3)
            beta  = tanh(p(4))
        else
            gamma = 0.0_dp
            beta  = tanh(p(3))
        end if
    end subroutine gas_transform

    subroutine gas_obj(p, np, f, g)
        ! Average NLL and analytical gradient w.r.t. unconstrained parameters.
        !
        ! Gradients are computed by propagating filter sensitivities delta_t^k =
        ! d(f_t)/d(p_k) alongside the filter, accumulating g_k -= s_t * delta_t^k.
        ! The lnorm_t normalising constant derivative uses a single scalar forward
        ! difference to avoid implementing the digamma function.
        integer,  intent(in)  :: np      ! number of parameters
        real(dp), intent(in)  :: p(np)   ! unconstrained parameters
        real(dp), intent(out) :: f       ! average NLL
        real(dp), intent(out) :: g(np)   ! gradient of average NLL

        real(dp), parameter :: h_fd = 1.0e-7_dp  ! step for dlnorm_t_dnu

        real(dp) :: omega, alpha, gamma, beta, beta2
        real(dp) :: nu, lnorm_t, dnu_dpnu, dlnorm_t_dnu
        real(dp) :: ft, ft_old, u2, ct, It, st, At
        real(dp) :: ds_dft, ds_dnu, dl_dnu_t, q
        real(dp) :: logl
        real(dp) :: delta(np)
        integer  :: t, i_beta

        call gas_transform(p, omega, alpha, gamma, beta)
        beta2  = beta**2
        i_beta = merge(4, 3, iproc_mod == proc_agas)

        nu           = 0.0_dp
        lnorm_t      = 0.0_dp
        dlnorm_t_dnu = 0.0_dp
        dnu_dpnu     = 0.0_dp
        if (idist_mod == dist_t) then
            nu       = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(np)))
            dnu_dpnu = (nu-2.0_dp) * (100.0_dp-nu) / 98.0_dp
            lnorm_t  = log_gamma(0.5_dp*(nu+1.0_dp)) - log_gamma(0.5_dp*nu) &
                     - 0.5_dp*log(pi*(nu-2.0_dp))
            dlnorm_t_dnu = (log_gamma(0.5_dp*(nu+h_fd+1.0_dp)) - log_gamma(0.5_dp*(nu+h_fd)) &
                          - 0.5_dp*log(pi*(nu+h_fd-2.0_dp)) - lnorm_t) / h_fd
        end if

        ! Initialise filter and sensitivities at unconditional log-variance
        ft            = omega / (1.0_dp - beta)
        delta         = 0.0_dp
        delta(1)      = 1.0_dp / (1.0_dp - beta)                     ! d(f_1)/d(omega)
        delta(i_beta) = omega * (1.0_dp + beta) / (1.0_dp - beta)    ! d(f_1)/d(p_beta)

        logl = 0.0_dp
        g    = 0.0_dp

        do t = 1, n_mod
            ft = min(max(ft, -500.0_dp), 500.0_dp)
            u2 = y_mod(t)**2 / exp(ft)
            It = merge(1.0_dp, 0.0_dp, y_mod(t) < 0.0_dp)
            ct = alpha + gamma * It

            select case (idist_mod)
            case (dist_normal)
                logl     = logl - log_sqrt_2pi - 0.5_dp*(ft + u2)
                st       = 0.5_dp*(u2 - 1.0_dp)
                At       = ct*(-0.5_dp*u2) + beta

                g        = g - st*delta
                ft_old   = ft
                ft       = omega + ct*st + beta*ft_old
                delta(1)      = 1.0_dp                + At*delta(1)
                delta(2)      = st                    + At*delta(2)
                if (iproc_mod == proc_agas) &
                    delta(3)  = st*It                 + At*delta(3)
                delta(i_beta) = (1.0_dp-beta2)*ft_old + At*delta(i_beta)

            case (dist_t)
                q        = nu - 2.0_dp + u2
                logl     = logl + lnorm_t - 0.5_dp*ft &
                               - 0.5_dp*(nu+1.0_dp)*log(1.0_dp + u2/(nu-2.0_dp))
                st       = (nu*u2 - (nu-2.0_dp)) / (2.0_dp*q)
                ds_dft   = -(nu+1.0_dp)*(nu-2.0_dp)*u2 / (2.0_dp*q**2)
                ds_dnu   = u2*(u2 - 3.0_dp) / (2.0_dp*q**2)
                dl_dnu_t = dlnorm_t_dnu - 0.5_dp*log(1.0_dp + u2/(nu-2.0_dp)) &
                         + 0.5_dp*(nu+1.0_dp)*u2 / ((nu-2.0_dp)*q)
                At       = ct*ds_dft + beta

                g        = g - st*delta
                g(np)    = g(np) - dl_dnu_t*dnu_dpnu
                ft_old   = ft
                ft       = omega + ct*st + beta*ft_old
                delta(1)      = 1.0_dp                    + At*delta(1)
                delta(2)      = st                        + At*delta(2)
                if (iproc_mod == proc_agas) &
                    delta(3)  = st*It                     + At*delta(3)
                delta(i_beta) = (1.0_dp-beta2)*ft_old     + At*delta(i_beta)
                delta(np)     = ct*ds_dnu*dnu_dpnu        + At*delta(np)
            end select
        end do

        f = -logl / n_mod
        g =     g / n_mod
    end subroutine gas_obj

    subroutine gas_simulate(omega, alpha, gamma_coef, beta, idist, nu, n, seed_val, y)
        ! Simulate n observations from a GAS(1,1) process.
        ! f_1 is initialised at the unconditional log-variance omega/(1-beta).
        ! Innovations are drawn from the standardised distribution (unit variance).
        ! Set gamma_coef=0 for the symmetric (proc_gas) variant.
        real(dp), intent(in)  :: omega      ! log-variance intercept
        real(dp), intent(in)  :: alpha      ! score coefficient
        real(dp), intent(in)  :: gamma_coef ! asymmetry weight (0 for symmetric)
        real(dp), intent(in)  :: beta       ! persistence, |beta| < 1
        integer,  intent(in)  :: idist      ! dist_normal or dist_t
        real(dp), intent(in)  :: nu         ! t degrees of freedom (ignored for dist_normal)
        integer,  intent(in)  :: n          ! number of observations
        integer,  intent(in)  :: seed_val   ! RNG seed
        real(dp), intent(out) :: y(n)       ! simulated returns

        real(dp) :: ft, h, eps, u2, It, ct, st, q
        integer  :: t, sz
        integer, allocatable :: seed_arr(:)

        call random_seed(size=sz)
        allocate(seed_arr(sz))
        seed_arr = seed_val
        call random_seed(put=seed_arr)
        deallocate(seed_arr)

        ft = omega / (1.0_dp - beta)
        do t = 1, n
            ft  = min(max(ft, -500.0_dp), 500.0_dp)
            h   = exp(ft)
            if (idist == dist_t) then
                eps = random_t_std(nu)
            else
                eps = random_normal()
            end if
            y(t) = sqrt(h) * eps
            u2   = eps**2
            It   = merge(1.0_dp, 0.0_dp, eps < 0.0_dp)
            ct   = alpha + gamma_coef * It

            if (idist == dist_normal) then
                st = 0.5_dp * (u2 - 1.0_dp)
            else
                q  = nu - 2.0_dp + u2
                st = (nu*u2 - (nu - 2.0_dp)) / (2.0_dp * q)
            end if

            ft = omega + ct*st + beta*ft
        end do
    end subroutine gas_simulate

end module gas_mod
