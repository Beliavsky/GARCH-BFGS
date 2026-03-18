! Special mathematical functions at double precision.
!
! digamma(x)           — d/dx log(Gamma(x)),  x > 0
! trigamma(x)          — d²/dx² log(Gamma(x)) = d/dx digamma(x),  x > 0
! bessel_k01(x,lk1,r) — log(K1(x)) and K0(x)/K1(x) for x > 0,
!                        where K0, K1 are modified Bessel functions of
!                        the second kind (not available as Fortran intrinsics).
! bessel_k_nu(nu, x)  — log(K_{nu}(x)) for x > 0, arbitrary real nu >= 0.
!                        Uses forward recurrence from K0, K1 for nu >= 0.

module special_mod
    use kind_mod,       only: dp
    use math_const_mod, only: pi
    implicit none
    private

    public :: digamma, trigamma, bessel_k01, bessel_k_nu

contains

    pure elemental function digamma(x) result(psi)
        ! Digamma function psi(x) = d/dx log(Gamma(x)) for x > 0.
        ! Uses recurrence psi(x) = psi(x+1) - 1/x to shift to x >= 8,
        ! then the asymptotic expansion (Abramowitz & Stegun 6.3.18).
        real(dp), intent(in) :: x   ! argument (> 0)
        real(dp) :: psi, xr
        xr  = x
        psi = 0.0_dp
        do while (xr < 8.0_dp)
            psi = psi - 1.0_dp / xr
            xr  = xr + 1.0_dp
        end do
        psi = psi + log(xr)               &
              - 1.0_dp/(  2.0_dp * xr   ) &
              - 1.0_dp/( 12.0_dp * xr**2) &
              + 1.0_dp/(120.0_dp * xr**4) &
              - 1.0_dp/(252.0_dp * xr**6)
    end function digamma

    pure elemental function trigamma(x) result(psi1)
        ! Trigamma function psi'(x) = d/dx digamma(x) for x > 0.
        ! Uses recurrence psi'(x) = psi'(x+1) + 1/x² to shift to x >= 8,
        ! then the asymptotic expansion (Abramowitz & Stegun 6.4.12).
        real(dp), intent(in) :: x   ! argument (> 0)
        real(dp) :: psi1, xr
        xr   = x
        psi1 = 0.0_dp
        do while (xr < 8.0_dp)
            psi1 = psi1 + 1.0_dp / xr**2
            xr   = xr + 1.0_dp
        end do
        psi1 = psi1 + 1.0_dp /        xr    &
               + 1.0_dp / ( 2.0_dp * xr**2) &
               + 1.0_dp / ( 6.0_dp * xr**3) &
               - 1.0_dp / (30.0_dp * xr**5) &
               + 1.0_dp / (42.0_dp * xr**7) &
               - 1.0_dp / (30.0_dp * xr**9)
    end function trigamma

    pure elemental subroutine bessel_k01(x, lk1, ratio)
        ! log(K1(x)) and K0(x)/K1(x) for x > 0.
        ! K0, K1 are modified Bessel functions of the second kind.
        ! Uses A&S 9.8 polynomial approximations for x <= 2 and 2 < x <= 600,
        ! and asymptotic expansion for x > 600 to avoid underflow.
        real(dp), intent(in)  :: x      ! argument (> 0)
        real(dp), intent(out) :: lk1    ! log(K1(x))
        real(dp), intent(out) :: ratio  ! K0(x)/K1(x)
        real(dp) :: t, p, bk0, bk1, bi0, bi1, lhx

        if (x > 600.0_dp) then
            t     = 1.0_dp / x
            lk1   = 0.5_dp*(log(pi) - log(2.0_dp) - log(x)) - x &
                    + log(1.0_dp + t*(0.375_dp + t*(-0.1171875_dp + t*0.1025390625_dp)))
            ratio = 1.0_dp - 0.5_dp*t
        else if (x <= 2.0_dp) then
            ! A&S 9.8.1, 9.8.3: I0, I1 (series with t=(x/3.75)^2)
            t   = (x / 3.75_dp)**2
            bi0 = 1.0_dp + t*(3.5156229_dp + t*(3.0899424_dp + &
                  t*(1.2067492_dp + t*(0.2659732_dp + t*(0.0360768_dp + &
                  t*0.0045813_dp)))))
            bi1 = x*(0.5_dp + t*(0.87890594_dp + t*(0.51498869_dp + &
                  t*(0.15084934_dp + t*(0.02658733_dp + t*(0.00301532_dp + &
                  t*0.00032411_dp))))))
            lhx = log(0.5_dp * x)
            ! A&S 9.8.5: K0 (t=(x/2)^2)
            t   = (0.5_dp * x)**2
            bk0 = -lhx*bi0 + (-0.57721566_dp + t*(0.42278420_dp + &
                  t*(0.23069756_dp + t*(0.03488590_dp + t*(0.00262698_dp + &
                  t*(0.00010750_dp + t*0.00000740_dp))))))
            ! A&S 9.8.7: K1
            bk1 = lhx*bi1 + (1.0_dp/x)*(1.0_dp + t*(0.15443144_dp + &
                  t*(-0.67278579_dp + t*(-0.18156897_dp + t*(-0.01919402_dp + &
                  t*(-0.00110404_dp + t*(-0.00004686_dp)))))))
            lk1   = log(bk1)
            ratio = bk0 / bk1
        else
            ! A&S 9.8.6, 9.8.8: K0, K1 (t=2/x; prefactor exp(-x)/sqrt(x))
            t   = 2.0_dp / x
            p   = exp(-x) / sqrt(x)
            bk0 = p * (1.25331414_dp + t*(-0.07832358_dp + t*(0.02189568_dp + &
                  t*(-0.01062446_dp + t*(0.00587872_dp + t*(-0.00251540_dp + &
                  t*0.00053208_dp))))))
            bk1 = p * (1.25331414_dp + t*(0.23498619_dp + t*(-0.03655620_dp + &
                  t*(0.01504268_dp + t*(-0.00780353_dp + t*(0.00325614_dp + &
                  t*(-0.00068245_dp)))))))
            lk1   = log(bk1)
            ratio = bk0 / bk1
        end if
    end subroutine bessel_k01

    pure function bessel_k_nu(nu, x) result(lknu)
        ! log(K_{nu}(x)) for x > 0 and nu >= 0 (real order).
        ! For integer or half-integer orders: forward recurrence from K0, K1
        ! via K_{n+1}(x) = (2n/x)*K_n(x) + K_{n-1}(x) up to floor(nu),
        ! then linear interpolation in order for the fractional remainder.
        ! For large x the forward recurrence is numerically stable (K grows).
        real(dp), intent(in) :: nu   ! order (>= 0)
        real(dp), intent(in) :: x    ! argument (> 0)
        real(dp) :: lknu
        real(dp) :: lk0, lk1_val, ratio, km1, k0r, k1r, kn
        integer  :: n, j
        real(dp) :: frac

        ! Get log(K0), K1/K0 ratio (= 1/ratio from bessel_k01 which returns K0/K1)
        call bessel_k01(x, lk1_val, ratio)
        ! ratio = K0/K1 from bessel_k01; lk1_val = log(K1)
        lk0 = lk1_val + log(ratio)   ! log(K0) = log(K1) + log(K0/K1)

        if (nu < 0.5_dp) then
            ! Interpolate between K0 (order 0) and K1 (order 1)
            lknu = (1.0_dp - nu) * lk0 + nu * lk1_val
            return
        end if

        n    = int(nu)               ! integer part
        frac = nu - real(n, dp)      ! fractional part in [0, 1)

        ! Forward recurrence: K_{n+1} = (2n/x)*K_n + K_{n-1}
        ! Work in log-scaled values to avoid overflow; carry the log of K_n
        ! and the ratio K_n / K_{n-1}.
        ! Use direct values for small n to avoid log precision issues.
        km1 = exp(lk0)     ! K_{n-1} = K_0
        k0r = exp(lk1_val) ! K_n     = K_1
        do j = 1, n - 1
            kn  = (2.0_dp * real(j, dp) / x) * k0r + km1
            km1 = k0r
            k0r = kn
        end do
        ! k0r = K_n (= K_{floor(nu)} when n >= 1)
        ! km1 = K_{n-1}

        if (n == 0) then
            ! nu < 1: interpolate between K0, K1
            lknu = (1.0_dp - frac) * lk0 + frac * lk1_val
        else if (frac < 1.0e-14_dp) then
            lknu = log(k0r)
        else
            ! Interpolate in order: log K_{n+frac} ≈ (1-frac)*log K_n + frac*log K_{n+1}
            k1r  = (2.0_dp * real(n, dp) / x) * k0r + km1   ! K_{n+1}
            lknu = (1.0_dp - frac) * log(k0r) + frac * log(k1r)
        end if
    end function bessel_k_nu

end module special_mod
