! Mathematical constants at double precision.
! All modules that need these constants should use this module
! rather than defining their own copies.

module math_const_mod
    use kind_mod, only: dp
    implicit none
    private

    real(dp), parameter, public :: pi           = 3.14159265358979323846_dp
    real(dp), parameter, public :: two_pi       = 2.0_dp * pi
    real(dp), parameter, public :: sqrt2        = 1.41421356237309504880_dp
    real(dp), parameter, public :: sqrt3        = 1.73205080756887729353_dp
    real(dp), parameter, public :: log2         = 0.69314718055994530942_dp
    real(dp), parameter, public :: half_log2    = 0.5_dp * log2            ! 0.5*log(2)
    real(dp), parameter, public :: log_sqrt_2pi = 0.91893853320467274178_dp ! 0.5*log(2*pi)

end module math_const_mod
