! Pseudo-random samplers for standardised distributions used as GARCH innovations.
!
! All scalar samplers are zero-argument or shape-parameter-only functions.
! Array overloads append integer dimension arguments (numpy-style):
!
!   z = random_normal()           ! scalar
!   z = random_normal(n)          ! rank-1, size n
!   z = random_normal(n, m)       ! rank-2, shape (n, m)
!
!   x = random_t_std(nu)          ! scalar; nu > 2
!   x = random_t_std(nu, n)       ! rank-1
!   x = random_t_std(nu, n, m)    ! rank-2
!
!   (same pattern for random_ged_std, random_laplace_std, random_logistic_std,
!    random_sech, random_nig_sym, random_nig, random_fs_skewt)
!
!   g = random_gamma(a)           ! scalar Gamma(a,1), a > 0 — building block
!   g = random_gamma(a, n)        ! rank-1
!   g = random_gamma(a, n, m)     ! rank-2
!
! Notes:
!   random_laplace_std, random_logistic_std, random_sech use inverse-CDF
!   (vectorised for array draws).  random_normal uses vectorised Box-Muller.
!   random_gamma uses Marsaglia-Tsang (rejection — looped for arrays).
!   random_nig_sym uses the normal-variance-mean mixture representation
!   with an InverseGaussian(1, alp^2) mixing variable (beta=0 special case).
!   random_nig(alp, bet) is the general skewed NIG (alp > |bet|):
!     gam2 = alp^2-bet^2,  V = (gam2/alp^2)*W,  W ~ InvGauss(1, gam2^2/alp^2),
!     x = bet*(V - gam2/alp^2) + sqrt(V)*Z.  Reduces to random_nig_sym when bet=0.
!   random_vg_sym(nu): symmetric VG, nu > 0.
!     G ~ Gamma(1/nu, nu),  x = sqrt(G)*Z.
!   random_vg(nu, rho): general skewed VG, nu > 0, rho in (-1,1).
!     theta = rho/sqrt(nu),  sigma = sqrt(1-rho^2),
!     G ~ Gamma(1/nu, nu),  x = theta*(G-1) + sigma*sqrt(G)*Z.

module random_mod
    use kind_mod,          only: dp
    use math_const_mod,    only: pi, two_pi, sqrt2, sqrt3
    use distributions_mod, only: ged_lambda, dist_normal, dist_t, dist_ged, dist_logistic, &
                                 dist_laplace, dist_sech, dist_nig_sym, dist_nig_gen, dist_fs_skewt
    implicit none
    private

    public :: random_normal
    public :: random_gamma
    public :: random_t_std
    public :: random_ged_std
    public :: random_laplace_std
    public :: random_logistic_std
    public :: random_sech
    public :: random_nig_sym
    public :: random_nig
    public :: random_fs_skewt
    public :: random_vg_sym
    public :: random_vg
    public :: random_dist_std

    interface random_normal
        module procedure random_normal_0
        module procedure random_normal_1
        module procedure random_normal_2
    end interface

    interface random_gamma
        module procedure random_gamma_0
        module procedure random_gamma_1
        module procedure random_gamma_2
    end interface

    interface random_t_std
        module procedure random_t_std_0
        module procedure random_t_std_1
        module procedure random_t_std_2
    end interface

    interface random_ged_std
        module procedure random_ged_std_0
        module procedure random_ged_std_1
        module procedure random_ged_std_2
    end interface

    interface random_laplace_std
        module procedure random_laplace_std_0
        module procedure random_laplace_std_1
        module procedure random_laplace_std_2
    end interface

    interface random_logistic_std
        module procedure random_logistic_std_0
        module procedure random_logistic_std_1
        module procedure random_logistic_std_2
    end interface

    interface random_sech
        module procedure random_sech_0
        module procedure random_sech_1
        module procedure random_sech_2
    end interface

    interface random_nig_sym
        module procedure random_nig_sym_0
        module procedure random_nig_sym_1
        module procedure random_nig_sym_2
    end interface

    interface random_nig
        module procedure random_nig_0
        module procedure random_nig_1
        module procedure random_nig_2
    end interface

    interface random_fs_skewt
        module procedure random_fs_skewt_0
        module procedure random_fs_skewt_1
        module procedure random_fs_skewt_2
    end interface

    interface random_vg_sym
        module procedure random_vg_sym_0
        module procedure random_vg_sym_1
        module procedure random_vg_sym_2
    end interface

    interface random_vg
        module procedure random_vg_0
        module procedure random_vg_1
        module procedure random_vg_2
    end interface

contains

    function random_dist_std(dist_id, shape, shape2) result(x)
        ! Draw one standardized variate from a supported distribution id.
        integer, intent(in) :: dist_id
        real(dp), intent(in), optional :: shape, shape2
        real(dp) :: x, s1, s2

        s1 = 0.0_dp
        s2 = 0.0_dp
        if (present(shape)) s1 = shape
        if (present(shape2)) s2 = shape2

        select case (dist_id)
        case (dist_normal)
            x = random_normal_0()
        case (dist_t)
            x = random_t_std_0(s1)
        case (dist_ged)
            x = random_ged_std_0(s1)
        case (dist_logistic)
            x = random_logistic_std_0()
        case (dist_laplace)
            x = random_laplace_std_0()
        case (dist_sech)
            x = random_sech_0()
        case (dist_nig_sym)
            x = random_nig_sym_0(s1)
        case (dist_nig_gen)
            x = random_nig_0(s1, s1*s2)
        case (dist_fs_skewt)
            x = random_fs_skewt_0(s1, s2)
        case default
            error stop "random_dist_std: unsupported distribution"
        end select
    end function random_dist_std

    ! ── Normal ────────────────────────────────────────────────────────────────

    function random_normal_0() result(z)
        ! One N(0,1) draw via Box-Muller.
        real(dp) :: z
        real(dp) :: u1, u2
        do
            call random_number(u1)
            if (u1 > 0.0_dp) exit
        end do
        call random_number(u2)
        z = sqrt(-2.0_dp * log(u1)) * cos(two_pi * u2)
    end function random_normal_0

    function random_normal_1(n) result(z)
        ! n N(0,1) draws via vectorised Box-Muller.
        integer, intent(in) :: n   ! number of draws
        real(dp) :: z(n), u1(n), u2(n)
        call random_number(u1)
        where (u1 == 0.0_dp) u1 = tiny(1.0_dp)
        call random_number(u2)
        z = sqrt(-2.0_dp * log(u1)) * cos(two_pi * u2)
    end function random_normal_1

    function random_normal_2(n, m) result(z)
        ! n×m array of N(0,1) draws via vectorised Box-Muller.
        integer, intent(in) :: n   ! number of rows
        integer, intent(in) :: m   ! number of columns
        real(dp) :: z(n, m), u1(n, m), u2(n, m)
        call random_number(u1)
        where (u1 == 0.0_dp) u1 = tiny(1.0_dp)
        call random_number(u2)
        z = sqrt(-2.0_dp * log(u1)) * cos(two_pi * u2)
    end function random_normal_2

    ! ── Gamma ─────────────────────────────────────────────────────────────────

    recursive function random_gamma_0(a) result(x)
        ! One Gamma(a,1) draw via Marsaglia-Tsang. Works for all a > 0.
        real(dp), intent(in) :: a   ! shape parameter
        real(dp) :: x
        real(dp) :: d, c, z, v, u
        if (a < 1.0_dp) then
            ! Gamma(a) = Gamma(a+1) * U^(1/a)
            x = random_gamma_0(a + 1.0_dp)
            call random_number(u)
            x = x * u**(1.0_dp / a)
            return
        end if
        d = a - 1.0_dp / 3.0_dp
        c = 1.0_dp / sqrt(9.0_dp * d)
        do
            z = random_normal_0()
            v = (1.0_dp + c*z)**3
            if (v > 0.0_dp) then
                call random_number(u)
                if (log(u) < 0.5_dp*z**2 + d*(1.0_dp - v + log(v))) then
                    x = d * v
                    return
                end if
            end if
        end do
    end function random_gamma_0

    function random_gamma_1(a, n) result(x)
        ! n Gamma(a,1) draws.
        real(dp), intent(in) :: a   ! shape parameter
        integer,  intent(in) :: n   ! number of draws
        real(dp) :: x(n)
        integer  :: i
        do i = 1, n
            x(i) = random_gamma_0(a)
        end do
    end function random_gamma_1

    function random_gamma_2(a, n, m) result(x)
        ! n×m array of Gamma(a,1) draws.
        real(dp), intent(in) :: a   ! shape parameter
        integer,  intent(in) :: n   ! number of rows
        integer,  intent(in) :: m   ! number of columns
        real(dp) :: x(n, m)
        integer  :: i, j
        do j = 1, m
            do i = 1, n
                x(i, j) = random_gamma_0(a)
            end do
        end do
    end function random_gamma_2

    ! ── Student-t ─────────────────────────────────────────────────────────────

    function random_t_std_0(nu) result(z)
        ! One draw from standardised t(nu): mean 0, variance 1. Requires nu > 2.
        real(dp), intent(in) :: nu   ! degrees of freedom (> 2)
        real(dp) :: z
        real(dp) :: g
        g = random_gamma_0(0.5_dp * nu)
        z = random_normal_0() / sqrt(2.0_dp * g / nu) * sqrt((nu - 2.0_dp) / nu)
    end function random_t_std_0

    function random_t_std_1(nu, n) result(z)
        ! n draws from standardised t(nu).
        real(dp), intent(in) :: nu   ! degrees of freedom (> 2)
        integer,  intent(in) :: n    ! number of draws
        real(dp) :: z(n)
        integer  :: i
        do i = 1, n
            z(i) = random_t_std_0(nu)
        end do
    end function random_t_std_1

    function random_t_std_2(nu, n, m) result(z)
        ! n×m array of draws from standardised t(nu).
        real(dp), intent(in) :: nu   ! degrees of freedom (> 2)
        integer,  intent(in) :: n    ! number of rows
        integer,  intent(in) :: m    ! number of columns
        real(dp) :: z(n, m)
        integer  :: i, j
        do j = 1, m
            do i = 1, n
                z(i, j) = random_t_std_0(nu)
            end do
        end do
    end function random_t_std_2

    ! ── GED ───────────────────────────────────────────────────────────────────

    function random_ged_std_0(nu) result(x)
        ! One draw from standardised GED(nu): mean 0, variance 1.
        ! Uses: x = lam*(2G)^(1/nu)*sign(u-0.5),  G ~ Gamma(1/nu, 1),  u ~ U(0,1).
        real(dp), intent(in) :: nu   ! shape (> 0; nu=2: Normal, nu=1: Laplace)
        real(dp) :: x
        real(dp) :: g, u
        g = random_gamma_0(1.0_dp / nu)
        call random_number(u)
        x = ged_lambda(nu) * (2.0_dp * g)**(1.0_dp / nu) * sign(1.0_dp, u - 0.5_dp)
    end function random_ged_std_0

    function random_ged_std_1(nu, n) result(x)
        ! n draws from standardised GED(nu).
        real(dp), intent(in) :: nu   ! shape parameter
        integer,  intent(in) :: n    ! number of draws
        real(dp) :: x(n)
        integer  :: i
        do i = 1, n
            x(i) = random_ged_std_0(nu)
        end do
    end function random_ged_std_1

    function random_ged_std_2(nu, n, m) result(x)
        ! n×m array of draws from standardised GED(nu).
        real(dp), intent(in) :: nu   ! shape parameter
        integer,  intent(in) :: n    ! number of rows
        integer,  intent(in) :: m    ! number of columns
        real(dp) :: x(n, m)
        integer  :: i, j
        do j = 1, m
            do i = 1, n
                x(i, j) = random_ged_std_0(nu)
            end do
        end do
    end function random_ged_std_2

    ! ── Laplace ───────────────────────────────────────────────────────────────

    function random_laplace_std_0() result(x)
        ! One draw from standardised Laplace: mean 0, variance 1.
        ! Inverse CDF: x = -sign(u-0.5)*log(1-2|u-0.5|)/sqrt(2).
        real(dp) :: x
        real(dp) :: u
        call random_number(u)
        u = max(min(u, 1.0_dp - tiny(1.0_dp)), tiny(1.0_dp))
        x = -sign(1.0_dp, u - 0.5_dp) * log(1.0_dp - 2.0_dp*abs(u - 0.5_dp)) / sqrt2
    end function random_laplace_std_0

    function random_laplace_std_1(n) result(x)
        ! n draws from standardised Laplace via vectorised inverse CDF.
        integer, intent(in) :: n   ! number of draws
        real(dp) :: x(n), u(n)
        call random_number(u)
        u = max(min(u, 1.0_dp - tiny(1.0_dp)), tiny(1.0_dp))
        x = -sign(1.0_dp, u - 0.5_dp) * log(1.0_dp - 2.0_dp*abs(u - 0.5_dp)) / sqrt2
    end function random_laplace_std_1

    function random_laplace_std_2(n, m) result(x)
        ! n×m array of draws from standardised Laplace via vectorised inverse CDF.
        integer, intent(in) :: n   ! number of rows
        integer, intent(in) :: m   ! number of columns
        real(dp) :: x(n, m), u(n, m)
        call random_number(u)
        u = max(min(u, 1.0_dp - tiny(1.0_dp)), tiny(1.0_dp))
        x = -sign(1.0_dp, u - 0.5_dp) * log(1.0_dp - 2.0_dp*abs(u - 0.5_dp)) / sqrt2
    end function random_laplace_std_2

    ! ── Logistic ──────────────────────────────────────────────────────────────

    function random_logistic_std_0() result(x)
        ! One draw from standardised Logistic: mean 0, variance 1.
        ! Inverse CDF: x = log(u/(1-u)) * sqrt(3)/pi.
        real(dp) :: x
        real(dp) :: u
        call random_number(u)
        u = max(min(u, 1.0_dp - tiny(1.0_dp)), tiny(1.0_dp))
        x = log(u / (1.0_dp - u)) * sqrt3 / pi
    end function random_logistic_std_0

    function random_logistic_std_1(n) result(x)
        ! n draws from standardised Logistic via vectorised inverse CDF.
        integer, intent(in) :: n   ! number of draws
        real(dp) :: x(n), u(n)
        call random_number(u)
        u = max(min(u, 1.0_dp - tiny(1.0_dp)), tiny(1.0_dp))
        x = log(u / (1.0_dp - u)) * sqrt3 / pi
    end function random_logistic_std_1

    function random_logistic_std_2(n, m) result(x)
        ! n×m array of draws from standardised Logistic via vectorised inverse CDF.
        integer, intent(in) :: n   ! number of rows
        integer, intent(in) :: m   ! number of columns
        real(dp) :: x(n, m), u(n, m)
        call random_number(u)
        u = max(min(u, 1.0_dp - tiny(1.0_dp)), tiny(1.0_dp))
        x = log(u / (1.0_dp - u)) * sqrt3 / pi
    end function random_logistic_std_2

    ! ── Hyperbolic secant ─────────────────────────────────────────────────────

    function random_sech_0() result(x)
        ! One draw from the hyperbolic secant distribution: mean 0, variance 1.
        ! Inverse CDF: x = (2/pi)*log(tan(pi*u/2)).
        real(dp) :: x
        real(dp) :: u
        call random_number(u)
        u = max(min(u, 1.0_dp - tiny(1.0_dp)), tiny(1.0_dp))
        x = (2.0_dp / pi) * log(tan(0.5_dp * pi * u))
    end function random_sech_0

    function random_sech_1(n) result(x)
        ! n draws from the hyperbolic secant distribution via vectorised inverse CDF.
        integer, intent(in) :: n   ! number of draws
        real(dp) :: x(n), u(n)
        call random_number(u)
        u = max(min(u, 1.0_dp - tiny(1.0_dp)), tiny(1.0_dp))
        x = (2.0_dp / pi) * log(tan(0.5_dp * pi * u))
    end function random_sech_1

    function random_sech_2(n, m) result(x)
        ! n×m array of draws from the hyperbolic secant distribution.
        integer, intent(in) :: n   ! number of rows
        integer, intent(in) :: m   ! number of columns
        real(dp) :: x(n, m), u(n, m)
        call random_number(u)
        u = max(min(u, 1.0_dp - tiny(1.0_dp)), tiny(1.0_dp))
        x = (2.0_dp / pi) * log(tan(0.5_dp * pi * u))
    end function random_sech_2

    ! ── Symmetric NIG ─────────────────────────────────────────────────────────

    function random_nig_sym_0(alp) result(x)
        ! One draw from symmetric unit-variance NIG(alp). Requires alp > 0.
        ! Uses normal variance-mean mixture: x = sqrt(V)*N(0,1),
        ! V ~ InverseGaussian(1, alp^2) via Michael-Schucany-Haas (1976).
        real(dp), intent(in) :: alp   ! tail/shape parameter (> 0)
        real(dp) :: x
        real(dp) :: v, y, u, lam2
        lam2 = alp**2
        y    = random_normal_0()**2
        v    = 1.0_dp + y / (2.0_dp*lam2) &
               - sqrt(4.0_dp*lam2*y + y**2) / (2.0_dp*lam2)
        call random_number(u)
        if (u > 1.0_dp / (1.0_dp + v)) v = 1.0_dp / v
        x = random_normal_0() * sqrt(v)
    end function random_nig_sym_0

    function random_nig_sym_1(alp, n) result(x)
        ! n draws from symmetric unit-variance NIG(alp).
        real(dp), intent(in) :: alp   ! tail/shape parameter (> 0)
        integer,  intent(in) :: n     ! number of draws
        real(dp) :: x(n)
        integer  :: i
        do i = 1, n
            x(i) = random_nig_sym_0(alp)
        end do
    end function random_nig_sym_1

    function random_nig_sym_2(alp, n, m) result(x)
        ! n×m array of draws from symmetric unit-variance NIG(alp).
        real(dp), intent(in) :: alp   ! tail/shape parameter (> 0)
        integer,  intent(in) :: n     ! number of rows
        integer,  intent(in) :: m     ! number of columns
        real(dp) :: x(n, m)
        integer  :: i, j
        do j = 1, m
            do i = 1, n
                x(i, j) = random_nig_sym_0(alp)
            end do
        end do
    end function random_nig_sym_2

    ! ── General NIG ───────────────────────────────────────────────────────────

    function random_nig_0(alp, bet) result(x)
        ! One draw from general unit-variance NIG(alp, bet).  Requires alp > |bet|.
        ! gam2 = alp^2-bet^2;  V = (gam2/alp^2)*W,  W ~ InvGauss(1, gam2^2/alp^2);
        ! x = bet*(V - gam2/alp^2) + sqrt(V)*Z.
        real(dp), intent(in) :: alp   ! tail parameter (> |bet|)
        real(dp), intent(in) :: bet   ! skewness parameter
        real(dp) :: x
        real(dp) :: alp2, gam2, mu_v, lam2, w, y, u
        alp2 = alp**2
        gam2 = alp2 - bet**2           ! gamma^2
        mu_v = gam2 / alp2             ! E[V] = gamma^2/alpha^2
        lam2 = gam2 * gam2 / alp2     ! InvGauss shape = gamma^4/alpha^2
        y = random_normal_0()**2
        w = 1.0_dp + y/(2.0_dp*lam2) - sqrt(4.0_dp*lam2*y + y**2)/(2.0_dp*lam2)
        call random_number(u)
        if (u > 1.0_dp/(1.0_dp + w)) w = 1.0_dp / w   ! W ~ InvGauss(1, lam2)
        w = mu_v * w                   ! V = mu_v*W ~ InvGauss(mu_v, mu_v^2*lam2)
        x = bet*(w - mu_v) + random_normal_0() * sqrt(w)
    end function random_nig_0

    function random_nig_1(alp, bet, n) result(x)
        ! n draws from general unit-variance NIG(alp, bet).
        real(dp), intent(in) :: alp   ! tail parameter (> |bet|)
        real(dp), intent(in) :: bet   ! skewness parameter
        integer,  intent(in) :: n     ! number of draws
        real(dp) :: x(n)
        integer  :: i
        do i = 1, n
            x(i) = random_nig_0(alp, bet)
        end do
    end function random_nig_1

    function random_nig_2(alp, bet, n, m) result(x)
        ! n×m array of draws from general unit-variance NIG(alp, bet).
        real(dp), intent(in) :: alp   ! tail parameter (> |bet|)
        real(dp), intent(in) :: bet   ! skewness parameter
        integer,  intent(in) :: n     ! number of rows
        integer,  intent(in) :: m     ! number of columns
        real(dp) :: x(n, m)
        integer  :: i, j
        do j = 1, m
            do i = 1, n
                x(i, j) = random_nig_0(alp, bet)
            end do
        end do
    end function random_nig_2

    function random_fs_skewt_0(nu, xi) result(x)
        real(dp), intent(in) :: nu, xi
        real(dp) :: x
        real(dp) :: u, zabs, xip, xim, m1, raw_mean, raw_second, raw_sd, y

        xip = max(xi, 1.0e-8_dp)
        xim = 1.0_dp / xip
        call random_number(u)
        zabs = abs(random_t_std_0(nu))
        if (u < xip / (xip + xim)) then
            y = xip * zabs
        else
            y = -xim * zabs
        end if
        m1 = sqrt(nu - 2.0_dp) * exp(log_gamma(0.5_dp*(nu - 1.0_dp)) - &
             0.5_dp*log(pi) - log_gamma(0.5_dp*nu))
        raw_mean = m1 * (xip - xim)
        raw_second = (xip**3 + xim**3) / (xip + xim)
        raw_sd = sqrt(max(raw_second - raw_mean**2, 1.0e-12_dp))
        x = (y - raw_mean) / raw_sd
    end function random_fs_skewt_0

    function random_fs_skewt_1(nu, xi, n) result(x)
        real(dp), intent(in) :: nu, xi
        integer, intent(in) :: n
        real(dp) :: x(n)
        integer :: i
        do i = 1, n
            x(i) = random_fs_skewt_0(nu, xi)
        end do
    end function random_fs_skewt_1

    function random_fs_skewt_2(nu, xi, n, m) result(x)
        real(dp), intent(in) :: nu, xi
        integer, intent(in) :: n, m
        real(dp) :: x(n, m)
        integer :: i, j
        do j = 1, m
            do i = 1, n
                x(i, j) = random_fs_skewt_0(nu, xi)
            end do
        end do
    end function random_fs_skewt_2

    ! ── Symmetric VG ──────────────────────────────────────────────────────────

    function random_vg_sym_0(nu) result(x)
        ! One draw from symmetric unit-variance VG(nu). Requires nu > 0.
        ! G ~ Gamma(1/nu, nu);  x = sqrt(G)*Z.
        real(dp), intent(in) :: nu   ! shape parameter (> 0)
        real(dp) :: x
        real(dp) :: g
        g = random_gamma_0(1.0_dp / nu) * nu   ! Gamma(1/nu, nu): scale by nu
        x = random_normal_0() * sqrt(g)
    end function random_vg_sym_0

    function random_vg_sym_1(nu, n) result(x)
        ! n draws from symmetric unit-variance VG(nu).
        real(dp), intent(in) :: nu   ! shape parameter (> 0)
        integer,  intent(in) :: n    ! number of draws
        real(dp) :: x(n)
        integer  :: i
        do i = 1, n
            x(i) = random_vg_sym_0(nu)
        end do
    end function random_vg_sym_1

    function random_vg_sym_2(nu, n, m) result(x)
        ! n×m array of draws from symmetric unit-variance VG(nu).
        real(dp), intent(in) :: nu   ! shape parameter (> 0)
        integer,  intent(in) :: n    ! number of rows
        integer,  intent(in) :: m    ! number of columns
        real(dp) :: x(n, m)
        integer  :: i, j
        do j = 1, m
            do i = 1, n
                x(i, j) = random_vg_sym_0(nu)
            end do
        end do
    end function random_vg_sym_2

    ! ── General VG ────────────────────────────────────────────────────────────

    function random_vg_0(nu, rho) result(x)
        ! One draw from general unit-variance VG(nu, rho).
        ! rho = theta*sqrt(nu) in (-1,1); theta = rho/sqrt(nu);
        ! sigma = sqrt(1-rho^2);  G ~ Gamma(1/nu, nu).
        ! x = theta*(G-1) + sigma*sqrt(G)*Z.
        real(dp), intent(in) :: nu    ! shape parameter (> 0)
        real(dp), intent(in) :: rho   ! skewness in (-1,1)
        real(dp) :: x
        real(dp) :: g, theta, sigma
        theta = rho / sqrt(nu)
        sigma = sqrt(1.0_dp - rho**2)
        g     = random_gamma_0(1.0_dp / nu) * nu
        x     = theta*(g - 1.0_dp) + sigma*sqrt(g)*random_normal_0()
    end function random_vg_0

    function random_vg_1(nu, rho, n) result(x)
        ! n draws from general unit-variance VG(nu, rho).
        real(dp), intent(in) :: nu    ! shape parameter (> 0)
        real(dp), intent(in) :: rho   ! skewness in (-1,1)
        integer,  intent(in) :: n     ! number of draws
        real(dp) :: x(n)
        integer  :: i
        do i = 1, n
            x(i) = random_vg_0(nu, rho)
        end do
    end function random_vg_1

    function random_vg_2(nu, rho, n, m) result(x)
        ! n×m array of draws from general unit-variance VG(nu, rho).
        real(dp), intent(in) :: nu    ! shape parameter (> 0)
        real(dp), intent(in) :: rho   ! skewness in (-1,1)
        integer,  intent(in) :: n     ! number of rows
        integer,  intent(in) :: m     ! number of columns
        real(dp) :: x(n, m)
        integer  :: i, j
        do j = 1, m
            do i = 1, n
                x(i, j) = random_vg_0(nu, rho)
            end do
        end do
    end function random_vg_2

end module random_mod
