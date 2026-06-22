program test_var_univariate
use random_mod, only: random_normal
use var_univariate_mod, only: dp, var_empirical, var_empirical_interp, &
   var_extrapolated_empirical, var_harrell_davis, var_trimmed_harrell_davis, &
   var_rectangular, var_triangular, var_kernel_gaussian, var_kernel_epanechnikov, &
   var_gaussian, var_student_t, var_cornish_fisher, var_gaussian_mixture, &
   var_evt_gpd_mom, var_evt_hill_weissman, var_mc_gaussian_hd, &
   var_mc_student_t_hd, var_log_to_arith
implicit none
integer, parameter :: n = 2500
real(dp)           :: r(n), alpha, mu, sig
real(dp)           :: wmix(2), mmix(2), smix(2)
integer            :: i, k
real(dp)           :: z

alpha = 0.99_dp
k = 125

call random_seed()
do i = 1, n
   z = random_normal()
   r(i) = 0.0003_dp + 0.012_dp*z
   if (mod(i, 97) == 0) r(i) = r(i) - 0.05_dp
end do

mu = sum(r)/real(n, dp)
sig = sqrt(sum((r - mu)**2)/real(n - 1, dp))
wmix = [0.92_dp, 0.08_dp]
mmix = [mu, mu - 0.025_dp]
smix = [0.010_dp, 0.035_dp]

write (*,'(a, f8.4)') "alpha = ", alpha
write (*,'(a, f12.6)') "empirical                  ", var_empirical(r, alpha)
write (*,'(a, f12.6)') "empirical interp           ", var_empirical_interp(r, alpha)
write (*,'(a, f12.6)') "extrapolated empirical     ", var_extrapolated_empirical(r, alpha)
write (*,'(a, f12.6)') "harrell davis              ", var_harrell_davis(r, alpha)
write (*,'(a, f12.6)') "trimmed harrell davis      ", var_trimmed_harrell_davis(r, alpha, 0.95_dp)
write (*,'(a, f12.6)') "rectangular                ", var_rectangular(r, alpha, 3)
write (*,'(a, f12.6)') "triangular                 ", var_triangular(r, alpha, 3)
write (*,'(a, f12.6)') "kernel gaussian            ", var_kernel_gaussian(r, alpha)
write (*,'(a, f12.6)') "kernel epanechnikov        ", var_kernel_epanechnikov(r, alpha)
write (*,'(a, f12.6)') "gaussian                   ", var_gaussian(r, alpha)
write (*,'(a, f12.6)') "student t nu=5             ", var_student_t(r, alpha, 5.0_dp)
write (*,'(a, f12.6)') "cornish fisher             ", var_cornish_fisher(r, alpha)
write (*,'(a, f12.6)') "gaussian mixture           ", var_gaussian_mixture(wmix, mmix, smix, alpha)
write (*,'(a, i0, a, f12.6)') "evt gpd mom k=", k, "           ", var_evt_gpd_mom(r, alpha, k)
write (*,'(a, i0, a, f12.6)') "evt hill k=", k, "              ", var_evt_hill_weissman(r, alpha, k)
write (*,'(a, f12.6)') "mc gaussian hd             ", var_mc_gaussian_hd(mu, sig, alpha, 5000)
write (*,'(a, f12.6)') "mc student t hd nu=5       ", var_mc_student_t_hd(mu, sig, 5.0_dp, alpha, 5000)
write (*,'(a, f12.6)') "arith from gaussian logvar ", var_log_to_arith(var_gaussian(r, alpha))
end program test_var_univariate
