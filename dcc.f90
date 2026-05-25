! DCC/ADCC correlation modules for two-stage multivariate GARCH estimation.
!
! All three stage-2 models share the same state (standardized residuals z_t,
! Q_bar = E[z_t z_t'], N_bar = E[n_t n_t'] where n_t = min(z_t, 0)).
! Stage-2 NLL routines use numerical gradients (2-4 parameters; cheap).
!
! DCC(1,1) Normal -- np=2: a, b
!   Q_t  = (1-a-b)*Q_bar + a*z_{t-1}z' + b*Q_{t-1}
!   NLL  = (0.5/T)*sum_t [log|R_t| + z_t'R_t^{-1}z_t - z_t'z_t]
!
! ADCC(1,1) Normal -- np=3: a, b, g
!   Q_t  = (1-a-b)*Q_bar - g*N_bar + a*z_{t-1}z' + b*Q_{t-1} + g*n_{t-1}n'
!   NLL  = same formula
!
! ADCC(1,1) Student-t -- np=4: a, b, g, nu>2
!   Q_t  = ADCC recursion
!   NLL  = (1/T)*sum_t [C_nu - 0.5*log|R_t| + 0.5*(nu+N)*log(1+z'R^{-1}z/(nu-2))]
!   where C_nu = log_gamma(nu/2) - log_gamma((nu+N)/2) + N/2*log(pi*(nu-2))
!   (sign: NLL to be minimised)
!
! Parameter transforms (multinomial logit ensures positivity and sum < 1):
!   DCC   : a = e^p1/D, b = e^p2/D, D = 1+e^p1+e^p2
!   ADCC  : a,b,g via D = 1+e^p1+e^p2+e^p3
!   ADCC-t: same a,b,g; nu = 2 + exp(p4)

module dcc_mod
    use kind_mod,   only: dp
    use linalg_mod, only: chol_factor, chol_logdet, chol_solve_vec
    implicit none
    private

    integer,               save :: dc_na = 0    ! n_assets
    integer,               save :: dc_nt = 0    ! n_obs
    real(dp), allocatable, save :: dc_z(:,:)    ! (na, nt) standardised residuals
    real(dp), allocatable, save :: dc_qbar(:,:) ! (na, na) sample covariance of z_t
    real(dp), allocatable, save :: dc_nbar(:,:) ! (na, na) sample covariance of n_t

    public :: dcc_set_resid
    public :: dcc_obj,    dcc_transform,    dcc_inv_transform
    public :: adcc_obj,   adcc_transform,   adcc_inv_transform
    public :: adcc_t_obj, adcc_t_transform, adcc_t_inv_transform

contains

    ! Store z_t = r_t/sqrt(h_t) and pre-compute Q_bar, N_bar.
    subroutine dcc_set_resid(z, na, nt)
        integer,  intent(in) :: na, nt     ! n_assets, n_obs
        real(dp), intent(in) :: z(na, nt)  ! z(:,t) = standardised return vector at t
        real(dp) :: zt(na), nv(na)
        integer  :: t, i, j
        if (allocated(dc_z))    deallocate(dc_z)
        if (allocated(dc_qbar)) deallocate(dc_qbar)
        if (allocated(dc_nbar)) deallocate(dc_nbar)
        allocate(dc_z(na,nt), dc_qbar(na,na), dc_nbar(na,na))
        dc_z    = z
        dc_na   = na
        dc_nt   = nt
        dc_qbar = 0.0_dp
        dc_nbar = 0.0_dp
        do t = 1, nt
            zt = z(:,t)
            nv = min(zt, 0.0_dp)
            do j = 1, na
                do i = 1, na
                    dc_qbar(i,j) = dc_qbar(i,j) + zt(i)*zt(j)
                    dc_nbar(i,j) = dc_nbar(i,j) + nv(i)*nv(j)
                end do
            end do
        end do
        dc_qbar = dc_qbar / real(nt, dp)
        dc_nbar = dc_nbar / real(nt, dp)
    end subroutine dcc_set_resid

    ! ── Parameter transforms ──────────────────────────────────────────────────

    ! DCC: p(1:2) -> (a, b) via multinomial logit.
    subroutine dcc_transform(p, np, a, b)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: a, b
        real(dp) :: denom
        denom = 1.0_dp + exp(p(1)) + exp(p(2))
        a = exp(p(1)) / denom
        b = exp(p(2)) / denom
    end subroutine dcc_transform

    subroutine dcc_inv_transform(a, b, p)
        real(dp), intent(in)  :: a, b
        real(dp), intent(out) :: p(2)
        real(dp), parameter :: eps = 1.0e-10_dp
        real(dp) :: lo
        lo   = max(1.0_dp - a - b, eps)
        p(1) = log(max(a, eps) / lo)
        p(2) = log(max(b, eps) / lo)
    end subroutine dcc_inv_transform

    ! ADCC: p(1:3) -> (a, b, g) via multinomial logit.
    subroutine adcc_transform(p, np, a, b, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: a, b, g
        real(dp) :: denom
        denom = 1.0_dp + exp(p(1)) + exp(p(2)) + exp(p(3))
        a = exp(p(1)) / denom
        b = exp(p(2)) / denom
        g = exp(p(3)) / denom
    end subroutine adcc_transform

    subroutine adcc_inv_transform(a, b, g, p)
        real(dp), intent(in)  :: a, b, g
        real(dp), intent(out) :: p(3)
        real(dp), parameter :: eps = 1.0e-10_dp
        real(dp) :: lo
        lo   = max(1.0_dp - a - b - g, eps)
        p(1) = log(max(a, eps) / lo)
        p(2) = log(max(b, eps) / lo)
        p(3) = log(max(g, eps) / lo)
    end subroutine adcc_inv_transform

    ! ADCC-t: p(1:3) -> (a,b,g) as above; p(4) -> nu = 2 + exp(p4).
    subroutine adcc_t_transform(p, np, a, b, g, nu)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: a, b, g, nu
        call adcc_transform(p(1:3), 3, a, b, g)
        nu = 2.0_dp + exp(p(4))
    end subroutine adcc_t_transform

    subroutine adcc_t_inv_transform(a, b, g, nu, p)
        real(dp), intent(in)  :: a, b, g, nu
        real(dp), intent(out) :: p(4)
        call adcc_inv_transform(a, b, g, p(1:3))
        p(4) = log(max(nu - 2.0_dp, 1.0e-10_dp))
    end subroutine adcc_t_inv_transform

    ! ── Shared inner kernel ───────────────────────────────────────────────────

    ! Given Q_t, compute R_t, its Cholesky L_t, log|R_t|, z'R^{-1}z, z'z.
    subroutine rt_stats(qt_in, na, zt, logdetR, zRz, zz, ok)
        integer,  intent(in)  :: na
        real(dp), intent(in)  :: qt_in(na,na), zt(na)
        real(dp), intent(out) :: logdetR, zRz, zz
        logical,  intent(out) :: ok
        real(dp) :: Rt(na,na), Lt(na,na), Rinvzt(na), qd(na)
        real(dp), parameter :: reg = 1.0e-8_dp
        integer  :: i, j
        do i = 1, na
            qd(i) = sqrt(max(qt_in(i,i), reg))
        end do
        do j = 1, na
            do i = 1, na
                Rt(i,j) = qt_in(i,j) / (qd(i)*qd(j))
            end do
        end do
        do i = 1, na
            Rt(i,i) = 1.0_dp   ! enforce exact unit diagonal
        end do
        call chol_factor(Rt, na, Lt, ok)
        if (.not. ok) return
        logdetR = chol_logdet(Lt, na)
        call chol_solve_vec(Lt, na, zt, Rinvzt)
        zRz = dot_product(zt, Rinvzt)
        zz  = dot_product(zt, zt)
    end subroutine rt_stats

    ! ── NLL routines ─────────────────────────────────────────────────────────

    ! DCC(1,1) Normal stage-2 correlation NLL/T.
    subroutine dcc_nll(p, np, f)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp) :: a, b, Qt(dc_na,dc_na), zt(dc_na)
        real(dp) :: logdetR, zRz, zz
        integer  :: t, i, j
        logical  :: ok
        call dcc_transform(p, np, a, b)
        Qt = dc_qbar
        f  = 0.0_dp
        do t = 1, dc_nt
            zt = dc_z(:,t)
            call rt_stats(Qt, dc_na, zt, logdetR, zRz, zz, ok)
            if (.not. ok) then; f = 1.0e30_dp; return; end if
            f = f + logdetR + zRz - zz
            do j = 1, dc_na
                do i = 1, dc_na
                    Qt(i,j) = (1.0_dp-a-b)*dc_qbar(i,j) + a*zt(i)*zt(j) + b*Qt(i,j)
                end do
            end do
        end do
        f = 0.5_dp * f / real(dc_nt, dp)
    end subroutine dcc_nll

    ! ADCC(1,1) Normal stage-2 correlation NLL/T.
    subroutine adcc_nll(p, np, f)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp) :: a, b, g, Qt(dc_na,dc_na), zt(dc_na), nv(dc_na)
        real(dp) :: logdetR, zRz, zz
        integer  :: t, i, j
        logical  :: ok
        call adcc_transform(p, np, a, b, g)
        Qt = dc_qbar
        f  = 0.0_dp
        do t = 1, dc_nt
            zt = dc_z(:,t)
            nv = min(zt, 0.0_dp)
            call rt_stats(Qt, dc_na, zt, logdetR, zRz, zz, ok)
            if (.not. ok) then; f = 1.0e30_dp; return; end if
            f = f + logdetR + zRz - zz
            do j = 1, dc_na
                do i = 1, dc_na
                    Qt(i,j) = (1.0_dp-a-b)*dc_qbar(i,j) - g*dc_nbar(i,j) &
                              + a*zt(i)*zt(j) + b*Qt(i,j) + g*nv(i)*nv(j)
                end do
            end do
        end do
        f = 0.5_dp * f / real(dc_nt, dp)
    end subroutine adcc_nll

    ! ADCC(1,1) Student-t stage-2 NLL/T.
    ! Includes log-Gamma terms and the heavy-tail factor; excludes -0.5*sum_i log(h_it)
    ! (that constant is handled by stage-1 accounting in the main program).
    subroutine adcc_t_nll(p, np, f)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp) :: a, b, g, nu, Qt(dc_na,dc_na), zt(dc_na), nv(dc_na)
        real(dp) :: logdetR, zRz, zz, rn, cn
        integer  :: t, i, j
        logical  :: ok
        real(dp), parameter :: pi = acos(-1.0_dp)
        call adcc_t_transform(p, np, a, b, g, nu)
        rn = real(dc_na, dp)
        ! Per-obs constant (sign: contributes positively to NLL to be minimised)
        cn = log_gamma(0.5_dp*nu) - log_gamma(0.5_dp*(nu+rn)) &
             + 0.5_dp*rn*log(pi*(nu-2.0_dp))
        Qt = dc_qbar
        f  = 0.0_dp
        do t = 1, dc_nt
            zt = dc_z(:,t)
            nv = min(zt, 0.0_dp)
            call rt_stats(Qt, dc_na, zt, logdetR, zRz, zz, ok)
            if (.not. ok) then; f = 1.0e30_dp; return; end if
            f = f + cn + 0.5_dp*logdetR &
                  + 0.5_dp*(nu+rn)*log(1.0_dp + zRz/(nu-2.0_dp))
            do j = 1, dc_na
                do i = 1, dc_na
                    Qt(i,j) = (1.0_dp-a-b)*dc_qbar(i,j) - g*dc_nbar(i,j) &
                              + a*zt(i)*zt(j) + b*Qt(i,j) + g*nv(i)*nv(j)
                end do
            end do
        end do
        f = f / real(dc_nt, dp)
    end subroutine adcc_t_nll

    ! ── BFGS wrappers with central-difference numerical gradient ─────────────

    subroutine dcc_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: pf(np), pb(np), ff, fb
        integer  :: j
        call dcc_nll(p, np, f)
        do j = 1, np
            pf = p; pf(j) = pf(j) + h
            pb = p; pb(j) = pb(j) - h
            call dcc_nll(pf, np, ff)
            call dcc_nll(pb, np, fb)
            g(j) = (ff - fb) / (2.0_dp * h)
        end do
    end subroutine dcc_obj

    subroutine adcc_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: pf(np), pb(np), ff, fb
        integer  :: j
        call adcc_nll(p, np, f)
        do j = 1, np
            pf = p; pf(j) = pf(j) + h
            pb = p; pb(j) = pb(j) - h
            call adcc_nll(pf, np, ff)
            call adcc_nll(pb, np, fb)
            g(j) = (ff - fb) / (2.0_dp * h)
        end do
    end subroutine adcc_obj

    subroutine adcc_t_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: pf(np), pb(np), ff, fb
        integer  :: j
        call adcc_t_nll(p, np, f)
        do j = 1, np
            pf = p; pf(j) = pf(j) + h
            pb = p; pb(j) = pb(j) - h
            call adcc_t_nll(pf, np, ff)
            call adcc_t_nll(pb, np, fb)
            g(j) = (ff - fb) / (2.0_dp * h)
        end do
    end subroutine adcc_t_obj

end module dcc_mod
