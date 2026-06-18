! by ChatGPT 5.5 based on articles by Roman Rubsamen at
! https://portfoliooptimizer.io/blog/value-at-risk-estimation-improved-historical-and-simulation-based-estimates-with-the-harrell-davis-quantile-estimator/
! https://portfoliooptimizer.io/blog/value-at-risk-univariate-estimation-methods/
module var_univariate_mod
use random_mod, only: random_normal, random_t_std
implicit none
private
public :: dp
public :: var_empirical, var_empirical_interp, var_extrapolated_empirical
public :: var_harrell_davis, var_trimmed_harrell_davis
public :: var_rectangular, var_triangular
public :: var_kernel_gaussian, var_kernel_epanechnikov
public :: var_gaussian, var_student_t, var_cornish_fisher, var_gaussian_mixture
public :: var_evt_gpd_mom, var_evt_hill_weissman
public :: var_mc_gaussian_hd, var_mc_student_t_hd
public :: var_log_to_arith

integer, parameter :: dp = kind(1.0d0)
contains

function var_empirical(r, alpha) result(var)
   ! VaR = -r_(floor(n*(1-alpha)) + 1), with r sorted ascending.
   real(dp), intent(in) :: r(:), alpha
   real(dp)             :: var
   real(dp), allocatable :: x(:)
   integer              :: n, j

   call check_alpha(alpha)
   n = size(r)
   call require(n >= 1, "var_empirical: empty input")

   x = r
   call sort_real(x)
   j = int(floor(real(n, dp)*(1.0_dp - alpha))) + 1
   j = max(1, min(n, j))
   var = -x(j)
end function var_empirical

function var_empirical_interp(r, alpha) result(var)
   ! Linear interpolation using p = 1-alpha and h = (n+1)*p.
   ! For h <= 1 or h >= n, clamp to the end observations.
   real(dp), intent(in) :: r(:), alpha
   real(dp)             :: var
   real(dp), allocatable :: x(:)
   integer              :: n, j
   real(dp)             :: h, gamma

   call check_alpha(alpha)
   n = size(r)
   call require(n >= 1, "var_empirical_interp: empty input")

   x = r
   call sort_real(x)
   h = real(n + 1, dp)*(1.0_dp - alpha)

   if (h <= 1.0_dp) then
      var = -x(1)
   else if (h >= real(n, dp)) then
      var = -x(n)
   else
      j = int(floor(h))
      gamma = h - real(j, dp)
      var = -((1.0_dp - gamma)*x(j) + gamma*x(j + 1))
   end if
end function var_empirical_interp

function var_extrapolated_empirical(r, alpha) result(var)
   ! Hutson-style extrapolated empirical quantile, lower-tail form.
   real(dp), intent(in) :: r(:), alpha
   real(dp)             :: var
   real(dp), allocatable :: x(:)
   integer              :: n
   real(dp)             :: p

   call check_alpha(alpha)
   n = size(r)
   call require(n >= 2, "var_extrapolated_empirical: need at least two observations")

   x = r
   call sort_real(x)
   p = 1.0_dp - alpha

   if (p <= 1.0_dp/real(n + 1, dp)) then
      var = -x(1) - (x(2) - x(1))*log(real(n + 1, dp)*p)
   else if (p < real(n, dp)/real(n + 1, dp)) then
      var = var_empirical_interp(x, alpha)
   else
      var = -x(n) + (x(n) - x(n - 1))*log(real(n + 1, dp)*alpha)
   end if
end function var_extrapolated_empirical

function var_harrell_davis(r, alpha) result(var)
   ! Harrell-Davis VaR = -sum_j w_j r_(j), w_j = I_{j/n}(a,b)-I_{(j-1)/n}(a,b),
   ! a = (n+1)*(1-alpha), b = (n+1)*alpha.
   real(dp), intent(in) :: r(:), alpha
   real(dp)             :: var
   real(dp), allocatable :: x(:)
   integer              :: n

   call check_alpha(alpha)
   n = size(r)
   call require(n >= 1, "var_harrell_davis: empty input")

   x = r
   call sort_real(x)
   var = -harrell_davis_quantile_sorted(x, 1.0_dp - alpha)
end function var_harrell_davis

function var_trimmed_harrell_davis(r, alpha, mass_width) result(var)
   ! Trimmed Harrell-Davis VaR. mass_width is the retained central beta-probability
   ! mass for the HD weights, e.g. 0.95. Weights outside the beta CDF interval
   ! [(1-mass_width)/2, (1+mass_width)/2] are set to zero and renormalized.
   real(dp), intent(in) :: r(:), alpha, mass_width
   real(dp)             :: var
   real(dp), allocatable :: x(:), w(:)
   integer              :: n, j
   real(dp)             :: p, a, b, c0, c1, lo, hi, s

   call check_alpha(alpha)
   call require(mass_width > 0.0_dp .and. mass_width <= 1.0_dp, &
                "var_trimmed_harrell_davis: mass_width must be in (0,1]")
   n = size(r)
   call require(n >= 1, "var_trimmed_harrell_davis: empty input")

   x = r
   call sort_real(x)
   p = 1.0_dp - alpha
   a = real(n + 1, dp)*p
   b = real(n + 1, dp)*(1.0_dp - p)
   lo = beta_inv(0.5_dp*(1.0_dp - mass_width), a, b)
   hi = beta_inv(0.5_dp*(1.0_dp + mass_width), a, b)

   allocate(w(n))
   s = 0.0_dp
   do j = 1, n
      c0 = max(real(j - 1, dp)/real(n, dp), lo)
      c1 = min(real(j, dp)/real(n, dp), hi)
      if (c1 > c0) then
         w(j) = beta_reg(c1, a, b) - beta_reg(c0, a, b)
      else
         w(j) = 0.0_dp
      end if
      s = s + w(j)
   end do
   call require(s > 0.0_dp, "var_trimmed_harrell_davis: zero retained weight")
   var = -sum(w*x)/s
end function var_trimmed_harrell_davis

function var_rectangular(r, alpha, half_width) result(var)
   ! Rectangular L-estimator: average order statistics in a window around the
   ! empirical quantile location.
   real(dp), intent(in) :: r(:), alpha
   integer, intent(in)  :: half_width
   real(dp)             :: var
   real(dp), allocatable :: x(:)
   integer              :: n, j, j1, j2

   call check_alpha(alpha)
   n = size(r)
   call require(n >= 1, "var_rectangular: empty input")
   call require(half_width >= 0, "var_rectangular: negative half_width")

   x = r
   call sort_real(x)
   j = int(floor(real(n, dp)*(1.0_dp - alpha))) + 1
   j = max(1, min(n, j))
   j1 = max(1, j - half_width)
   j2 = min(n, j + half_width)
   var = -sum(x(j1:j2))/real(j2 - j1 + 1, dp)
end function var_rectangular

function var_triangular(r, alpha, half_width) result(var)
   ! Triangular L-estimator: triangular weights around the empirical quantile
   ! location, with radius half_width.
   real(dp), intent(in) :: r(:), alpha
   integer, intent(in)  :: half_width
   real(dp)             :: var
   real(dp), allocatable :: x(:)
   integer              :: n, j, j1, j2, i
   real(dp)             :: w, s, q

   call check_alpha(alpha)
   n = size(r)
   call require(n >= 1, "var_triangular: empty input")
   call require(half_width >= 0, "var_triangular: negative half_width")

   x = r
   call sort_real(x)
   j = int(floor(real(n, dp)*(1.0_dp - alpha))) + 1
   j = max(1, min(n, j))
   j1 = max(1, j - half_width)
   j2 = min(n, j + half_width)
   q = 0.0_dp
   s = 0.0_dp
   do i = j1, j2
      w = real(half_width + 1 - abs(i - j), dp)
      q = q + w*x(i)
      s = s + w
   end do
   var = -q/s
end function var_triangular

function var_kernel_gaussian(r, alpha, h) result(var)
   ! Gaussian-kernel smoothed VaR. Solves mean Phi((x-r_i)/h) = 1-alpha.
   ! If h <= 0, use Silverman's rule: 1.06*sd*n^(-1/5).
   real(dp), intent(in) :: r(:), alpha
   real(dp), intent(in), optional :: h
   real(dp)             :: var
   real(dp)             :: bw, lo, hi, mid, p, sx, xmin, xmax
   integer              :: n, iter

   call check_alpha(alpha)
   n = size(r)
   call require(n >= 2, "var_kernel_gaussian: need at least two observations")

   bw = optional_bandwidth(r, h)
   p = 1.0_dp - alpha
   xmin = minval(r)
   xmax = maxval(r)
   sx = sample_sd(r)
   lo = xmin - 12.0_dp*max(bw, sx, 1.0e-12_dp)
   hi = xmax + 12.0_dp*max(bw, sx, 1.0e-12_dp)

   do iter = 1, 120
      mid = 0.5_dp*(lo + hi)
      if (kernel_cdf_gaussian(r, mid, bw) < p) then
         lo = mid
      else
         hi = mid
      end if
   end do

   var = -0.5_dp*(lo + hi)
end function var_kernel_gaussian

function var_kernel_epanechnikov(r, alpha, h) result(var)
   ! Epanechnikov-kernel smoothed VaR. Solves the corresponding smoothed CDF.
   ! If h <= 0, use a Silverman-style bandwidth.
   real(dp), intent(in) :: r(:), alpha
   real(dp), intent(in), optional :: h
   real(dp)             :: var
   real(dp)             :: bw, lo, hi, mid, p, xmin, xmax
   integer              :: n, iter

   call check_alpha(alpha)
   n = size(r)
   call require(n >= 2, "var_kernel_epanechnikov: need at least two observations")

   bw = optional_bandwidth(r, h)
   p = 1.0_dp - alpha
   xmin = minval(r)
   xmax = maxval(r)
   lo = xmin - bw
   hi = xmax + bw

   do iter = 1, 120
      mid = 0.5_dp*(lo + hi)
      if (kernel_cdf_epanechnikov(r, mid, bw) < p) then
         lo = mid
      else
         hi = mid
      end if
   end do

   var = -0.5_dp*(lo + hi)
end function var_kernel_epanechnikov

function var_gaussian(r, alpha) result(var)
   ! Gaussian VaR = -mu - sigma*z_(1-alpha).
   real(dp), intent(in) :: r(:), alpha
   real(dp)             :: var, mu, sig, z

   call check_alpha(alpha)
   call require(size(r) >= 2, "var_gaussian: need at least two observations")

   mu = sum(r)/real(size(r), dp)
   sig = sample_sd(r)
   z = norm_inv(1.0_dp - alpha)
   var = -mu - sig*z
end function var_gaussian

function var_student_t(r, alpha, nu) result(var)
   ! Student-t parametric VaR. Uses sample mean and sd and a standardized t_nu
   ! innovation with variance 1, so nu must be greater than 2.
   real(dp), intent(in) :: r(:), alpha, nu
   real(dp)             :: var, mu, sig, z

   call check_alpha(alpha)
   call require(size(r) >= 2, "var_student_t: need at least two observations")
   call require(nu > 2.0_dp, "var_student_t: nu must be greater than 2")

   mu = sum(r)/real(size(r), dp)
   sig = sample_sd(r)
   z = student_t_inv(1.0_dp - alpha, nu)*sqrt((nu - 2.0_dp)/nu)
   var = -mu - sig*z
end function var_student_t

function var_cornish_fisher(r, alpha) result(var)
   ! Modified VaR using a Cornish-Fisher adjusted normal quantile.
   real(dp), intent(in) :: r(:), alpha
   real(dp)             :: var, mu, sig, z, zcf, skew, exk

   call check_alpha(alpha)
   call require(size(r) >= 4, "var_cornish_fisher: need at least four observations")

   mu = sum(r)/real(size(r), dp)
   sig = sample_sd(r)
   skew = sample_skew(r)
   exk = sample_excess_kurt(r)
   z = norm_inv(1.0_dp - alpha)
   zcf = z + (z*z - 1.0_dp)*skew/6.0_dp + &
         (z**3 - 3.0_dp*z)*exk/24.0_dp - &
         (2.0_dp*z**3 - 5.0_dp*z)*skew*skew/36.0_dp
   var = -mu - sig*zcf
end function var_cornish_fisher

function var_gaussian_mixture(w, mu, sig, alpha) result(var)
   ! Parametric Gaussian-mixture VaR for supplied mixture parameters.
   real(dp), intent(in) :: w(:), mu(:), sig(:), alpha
   real(dp)             :: var
   integer              :: m, iter
   real(dp)             :: p, sw, lo, hi, mid, scale

   call check_alpha(alpha)
   m = size(w)
   call require(m >= 1, "var_gaussian_mixture: empty mixture")
   call require(size(mu) == m .and. size(sig) == m, "var_gaussian_mixture: shape mismatch")
   call require(all(w >= 0.0_dp), "var_gaussian_mixture: negative weight")
   call require(all(sig > 0.0_dp), "var_gaussian_mixture: non-positive sigma")
   sw = sum(w)
   call require(sw > 0.0_dp, "var_gaussian_mixture: zero total weight")

   p = 1.0_dp - alpha
   scale = max(maxval(sig), 1.0e-12_dp)
   lo = minval(mu - 12.0_dp*sig) - 12.0_dp*scale
   hi = maxval(mu + 12.0_dp*sig) + 12.0_dp*scale

   do iter = 1, 140
      mid = 0.5_dp*(lo + hi)
      if (gaussian_mixture_cdf(mid, w, mu, sig)/sw < p) then
         lo = mid
      else
         hi = mid
      end if
   end do
   var = -0.5_dp*(lo + hi)
end function var_gaussian_mixture

function var_evt_gpd_mom(r, alpha, k) result(var)
   ! EVT/GPD VaR on losses = -returns. Uses top k losses and MOM GPD fit.
   ! Practical code should study sensitivity to k and prefer MLE/PWM if available.
   real(dp), intent(in) :: r(:), alpha
   integer, intent(in)  :: k
   real(dp)             :: var
   real(dp), allocatable :: loss(:), y(:)
   integer              :: n, i
   real(dp)             :: p, u, m1, m2, xi, sig, ratio

   call check_alpha(alpha)
   n = size(r)
   call require(n >= 3, "var_evt_gpd_mom: need at least three observations")
   call require(k >= 2 .and. k < n, "var_evt_gpd_mom: require 2 <= k < n")

   loss = -r
   call sort_real(loss)
   u = loss(n - k)
   allocate(y(k))
   do i = 1, k
      y(i) = loss(n - k + i) - u
   end do

   m1 = sum(y)/real(k, dp)
   m2 = sum((y - m1)**2)/real(k, dp)
   call require(m1 > 0.0_dp .and. m2 > 0.0_dp, "var_evt_gpd_mom: degenerate tail excesses")

   xi = 0.5_dp*(1.0_dp - (m1*m1)/m2)
   sig = 0.5_dp*m1*(1.0_dp + (m1*m1)/m2)
   call require(sig > 0.0_dp, "var_evt_gpd_mom: non-positive scale estimate")

   p = 1.0_dp - alpha
   ratio = real(k, dp)/(real(n, dp)*p)
   call require(ratio > 0.0_dp, "var_evt_gpd_mom: invalid tail probability")

   if (abs(xi) < 1.0e-10_dp) then
      var = u + sig*log(ratio)
   else
      var = u + (sig/xi)*(ratio**xi - 1.0_dp)
   end if
end function var_evt_gpd_mom

function var_evt_hill_weissman(r, alpha, k) result(var)
   ! EVT/Hill-Weissman VaR on positive losses = -returns.
   ! Implements the common Hill form: VaR = u*(k/(n*p))^hill_xi.
   real(dp), intent(in) :: r(:), alpha
   integer, intent(in)  :: k
   real(dp)             :: var
   real(dp), allocatable :: loss(:)
   integer              :: n, i
   real(dp)             :: p, u, hill_xi, ratio

   call check_alpha(alpha)
   n = size(r)
   call require(n >= 3, "var_evt_hill_weissman: need at least three observations")
   call require(k >= 1 .and. k < n, "var_evt_hill_weissman: require 1 <= k < n")

   loss = -r
   call sort_real(loss)
   u = loss(n - k)
   call require(u > 0.0_dp, "var_evt_hill_weissman: threshold loss must be positive")

   hill_xi = 0.0_dp
   do i = 1, k
      call require(loss(n - i + 1) > 0.0_dp, "var_evt_hill_weissman: tail losses must be positive")
      hill_xi = hill_xi + log(loss(n - i + 1)/u)
   end do
   hill_xi = hill_xi/real(k, dp)
   call require(hill_xi > 0.0_dp, "var_evt_hill_weissman: non-positive Hill estimate")

   p = 1.0_dp - alpha
   ratio = real(k, dp)/(real(n, dp)*p)
   var = u*ratio**hill_xi
end function var_evt_hill_weissman

function var_mc_gaussian_hd(mu, sig, alpha, nsim) result(var)
   ! Simulation-based Gaussian VaR using the Harrell-Davis quantile estimator.
   real(dp), intent(in) :: mu, sig, alpha
   integer, intent(in)  :: nsim
   real(dp)             :: var
   real(dp), allocatable :: r(:)
   integer              :: i

   call check_alpha(alpha)
   call require(sig > 0.0_dp, "var_mc_gaussian_hd: sig must be positive")
   call require(nsim >= 2, "var_mc_gaussian_hd: nsim must be at least 2")

   allocate(r(nsim))
   do i = 1, nsim
      r(i) = mu + sig*random_normal()
   end do
   var = var_harrell_davis(r, alpha)
end function var_mc_gaussian_hd

function var_mc_student_t_hd(mu, sig, nu, alpha, nsim) result(var)
   ! Simulation-based standardized Student-t VaR using the Harrell-Davis quantile estimator.
   real(dp), intent(in) :: mu, sig, nu, alpha
   integer, intent(in)  :: nsim
   real(dp)             :: var
   real(dp), allocatable :: r(:)
   integer              :: i

   call check_alpha(alpha)
   call require(sig > 0.0_dp, "var_mc_student_t_hd: sig must be positive")
   call require(nu > 2.0_dp, "var_mc_student_t_hd: nu must be greater than 2")
   call require(nsim >= 2, "var_mc_student_t_hd: nsim must be at least 2")

   allocate(r(nsim))
   do i = 1, nsim
      r(i) = mu + sig*random_t_std(nu)
   end do
   var = var_harrell_davis(r, alpha)
end function var_mc_student_t_hd

pure elemental function var_log_to_arith(var_log) result(var_arith)
   ! Convert log-return loss VaR to arithmetic-return loss VaR.
   real(dp), intent(in) :: var_log
   real(dp)             :: var_arith
   var_arith = 1.0_dp - exp(-var_log)
end function var_log_to_arith

function harrell_davis_quantile_sorted(x, p) result(q)
   real(dp), intent(in) :: x(:), p
   real(dp)             :: q, a, b, w
   integer              :: n, j

   call require(p > 0.0_dp .and. p < 1.0_dp, "harrell_davis_quantile_sorted: p must be in (0,1)")
   n = size(x)
   a = real(n + 1, dp)*p
   b = real(n + 1, dp)*(1.0_dp - p)
   q = 0.0_dp
   do j = 1, n
      w = beta_reg(real(j, dp)/real(n, dp), a, b) - &
          beta_reg(real(j - 1, dp)/real(n, dp), a, b)
      q = q + w*x(j)
   end do
end function harrell_davis_quantile_sorted

function optional_bandwidth(r, h) result(bw)
   real(dp), intent(in) :: r(:)
   real(dp), intent(in), optional :: h
   real(dp)             :: bw
   integer              :: n

   n = size(r)
   if (present(h)) then
      bw = h
   else
      bw = 0.0_dp
   end if
   if (bw <= 0.0_dp) then
      bw = 1.06_dp*sample_sd(r)*real(n, dp)**(-0.2_dp)
   end if
   call require(bw > 0.0_dp, "optional_bandwidth: bandwidth is not positive")
end function optional_bandwidth

function kernel_cdf_gaussian(x, q, h) result(cdf)
   real(dp), intent(in) :: x(:), q, h
   real(dp)             :: cdf
   integer              :: i

   cdf = 0.0_dp
   do i = 1, size(x)
      cdf = cdf + norm_cdf((q - x(i))/h)
   end do
   cdf = cdf/real(size(x), dp)
end function kernel_cdf_gaussian

function kernel_cdf_epanechnikov(x, q, h) result(cdf)
   real(dp), intent(in) :: x(:), q, h
   real(dp)             :: cdf, u
   integer              :: i

   cdf = 0.0_dp
   do i = 1, size(x)
      u = (q - x(i))/h
      if (u <= -1.0_dp) then
         cdf = cdf
      else if (u >= 1.0_dp) then
         cdf = cdf + 1.0_dp
      else
         cdf = cdf + 0.5_dp + 0.75_dp*(u - u**3/3.0_dp)
      end if
   end do
   cdf = cdf/real(size(x), dp)
end function kernel_cdf_epanechnikov

function gaussian_mixture_cdf(x, w, mu, sig) result(cdf)
   real(dp), intent(in) :: x, w(:), mu(:), sig(:)
   real(dp)             :: cdf
   integer              :: i

   cdf = 0.0_dp
   do i = 1, size(w)
      cdf = cdf + w(i)*norm_cdf((x - mu(i))/sig(i))
   end do
end function gaussian_mixture_cdf

pure function norm_cdf(x) result(p)
   real(dp), intent(in) :: x
   real(dp)             :: p
   p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
end function norm_cdf

pure function norm_inv(p) result(x)
   ! Peter J. Acklam's rational approximation for the standard normal quantile.
   real(dp), intent(in) :: p
   real(dp)             :: x, q, r
   real(dp), parameter  :: a1 = -3.969683028665376d+01, a2 =  2.209460984245205d+02
   real(dp), parameter  :: a3 = -2.759285104469687d+02, a4 =  1.383577518672690d+02
   real(dp), parameter  :: a5 = -3.066479806614716d+01, a6 =  2.506628277459239d+00
   real(dp), parameter  :: b1 = -5.447609879822406d+01, b2 =  1.615858368580409d+02
   real(dp), parameter  :: b3 = -1.556989798598866d+02, b4 =  6.680131188771972d+01
   real(dp), parameter  :: b5 = -1.328068155288572d+01
   real(dp), parameter  :: c1 = -7.784894002430293d-03, c2 = -3.223964580411365d-01
   real(dp), parameter  :: c3 = -2.400758277161838d+00, c4 = -2.549732539343734d+00
   real(dp), parameter  :: c5 =  4.374664141464968d+00, c6 =  2.938163982698783d+00
   real(dp), parameter  :: d1 =  7.784695709041462d-03, d2 =  3.224671290700398d-01
   real(dp), parameter  :: d3 =  2.445134137142996d+00, d4 =  3.754408661907416d+00
   real(dp), parameter  :: plow = 0.02425_dp, phigh = 1.0_dp - plow

   if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
   else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
   else if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
          ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_dp)
   else if (p <= phigh) then
      q = p - 0.5_dp
      r = q*q
      x = (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6)*q / &
          (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1.0_dp)
   else
      q = sqrt(-2.0_dp*log(1.0_dp - p))
      x = -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
           ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_dp)
   end if
end function norm_inv

function student_t_cdf(x, nu) result(p)
   real(dp), intent(in) :: x, nu
   real(dp)             :: p, y

   call require(nu > 0.0_dp, "student_t_cdf: nu must be positive")
   if (abs(x) <= 0.0_dp) then
      p = 0.5_dp
   else
      y = beta_reg(nu/(nu + x*x), 0.5_dp*nu, 0.5_dp)
      if (x > 0.0_dp) then
         p = 1.0_dp - 0.5_dp*y
      else
         p = 0.5_dp*y
      end if
   end if
end function student_t_cdf

function student_t_inv(p, nu) result(x)
   real(dp), intent(in) :: p, nu
   real(dp)             :: x, lo, hi, mid
   integer              :: iter

   call require(p > 0.0_dp .and. p < 1.0_dp, "student_t_inv: p must be in (0,1)")
   call require(nu > 0.0_dp, "student_t_inv: nu must be positive")

   lo = -1.0_dp
   hi = 1.0_dp
   do while (student_t_cdf(lo, nu) > p)
      lo = 2.0_dp*lo
      call require(abs(lo) < 1.0e8_dp, "student_t_inv: lower bracket failed")
   end do
   do while (student_t_cdf(hi, nu) < p)
      hi = 2.0_dp*hi
      call require(abs(hi) < 1.0e8_dp, "student_t_inv: upper bracket failed")
   end do

   do iter = 1, 120
      mid = 0.5_dp*(lo + hi)
      if (student_t_cdf(mid, nu) < p) then
         lo = mid
      else
         hi = mid
      end if
   end do
   x = 0.5_dp*(lo + hi)
end function student_t_inv

function beta_reg(x, a, b) result(bt)
   ! Regularized incomplete beta I_x(a,b), based on a continued fraction.
   real(dp), intent(in) :: x, a, b
   real(dp)             :: bt, front

   call require(a > 0.0_dp .and. b > 0.0_dp, "beta_reg: a and b must be positive")
   call require(x >= 0.0_dp .and. x <= 1.0_dp, "beta_reg: x must be in [0,1]")

   if (x <= 0.0_dp) then
      bt = 0.0_dp
   else if (x >= 1.0_dp) then
      bt = 1.0_dp
   else
      front = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
                  a*log(x) + b*log(1.0_dp - x))
      if (x < (a + 1.0_dp)/(a + b + 2.0_dp)) then
         bt = front*betacf(a, b, x)/a
      else
         bt = 1.0_dp - front*betacf(b, a, 1.0_dp - x)/b
      end if
      bt = max(0.0_dp, min(1.0_dp, bt))
   end if
end function beta_reg

function beta_inv(p, a, b) result(x)
   real(dp), intent(in) :: p, a, b
   real(dp)             :: x, lo, hi, mid
   integer              :: iter

   call require(p >= 0.0_dp .and. p <= 1.0_dp, "beta_inv: p must be in [0,1]")
   if (p <= 0.0_dp) then
      x = 0.0_dp
      return
   else if (p >= 1.0_dp) then
      x = 1.0_dp
      return
   end if

   lo = 0.0_dp
   hi = 1.0_dp
   do iter = 1, 120
      mid = 0.5_dp*(lo + hi)
      if (beta_reg(mid, a, b) < p) then
         lo = mid
      else
         hi = mid
      end if
   end do
   x = 0.5_dp*(lo + hi)
end function beta_inv

function betacf(a, b, x) result(cf)
   real(dp), intent(in) :: a, b, x
   real(dp)             :: cf
   integer, parameter   :: maxit = 300
   real(dp), parameter  :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
   integer              :: m, m2
   real(dp)             :: aa, c, d, del, h, qab, qam, qap

   qab = a + b
   qap = a + 1.0_dp
   qam = a - 1.0_dp
   c = 1.0_dp
   d = 1.0_dp - qab*x/qap
   if (abs(d) < fpmin) d = fpmin
   d = 1.0_dp/d
   h = d
   do m = 1, maxit
      m2 = 2*m
      aa = real(m, dp)*(b - real(m, dp))*x/((qam + real(m2, dp))*(a + real(m2, dp)))
      d = 1.0_dp + aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      h = h*d*c
      aa = -(a + real(m, dp))*(qab + real(m, dp))*x/((a + real(m2, dp))*(qap + real(m2, dp)))
      d = 1.0_dp + aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      del = d*c
      h = h*del
      if (abs(del - 1.0_dp) <= eps) exit
   end do
   cf = h
end function betacf

subroutine sort_real(x)
   real(dp), intent(inout) :: x(:)
   call quicksort(x, 1, size(x))
end subroutine sort_real

recursive subroutine quicksort(x, left, right)
   real(dp), intent(inout) :: x(:)
   integer, intent(in)     :: left, right
   integer                 :: i, j
   real(dp)                :: pivot, tmp

   if (left >= right) return
   i = left
   j = right
   pivot = x((left + right)/2)
   do
      do while (x(i) < pivot)
         i = i + 1
      end do
      do while (pivot < x(j))
         j = j - 1
      end do
      if (i <= j) then
         tmp = x(i); x(i) = x(j); x(j) = tmp
         i = i + 1
         j = j - 1
      end if
      if (i > j) exit
   end do
   if (left < j) call quicksort(x, left, j)
   if (i < right) call quicksort(x, i, right)
end subroutine quicksort

function sample_sd(x) result(sd)
   real(dp), intent(in) :: x(:)
   real(dp)             :: sd, mu
   integer              :: n

   n = size(x)
   call require(n >= 2, "sample_sd: need at least two observations")
   mu = sum(x)/real(n, dp)
   sd = sqrt(sum((x - mu)**2)/real(n - 1, dp))
end function sample_sd

function sample_skew(x) result(skew)
   real(dp), intent(in) :: x(:)
   real(dp)             :: skew, mu, s
   integer              :: n

   n = size(x)
   mu = sum(x)/real(n, dp)
   s = sqrt(sum((x - mu)**2)/real(n, dp))
   call require(s > 0.0_dp, "sample_skew: zero variance")
   skew = sum(((x - mu)/s)**3)/real(n, dp)
end function sample_skew

function sample_excess_kurt(x) result(exk)
   real(dp), intent(in) :: x(:)
   real(dp)             :: exk, mu, s
   integer              :: n

   n = size(x)
   mu = sum(x)/real(n, dp)
   s = sqrt(sum((x - mu)**2)/real(n, dp))
   call require(s > 0.0_dp, "sample_excess_kurt: zero variance")
   exk = sum(((x - mu)/s)**4)/real(n, dp) - 3.0_dp
end function sample_excess_kurt

subroutine check_alpha(alpha)
   real(dp), intent(in) :: alpha
   call require(alpha > 0.0_dp .and. alpha < 1.0_dp, "alpha must be in (0,1)")
end subroutine check_alpha

subroutine require(ok, msg)
   logical, intent(in)          :: ok
   character(len=*), intent(in) :: msg
   if (.not. ok) then
      write (*,*) trim(msg)
      error stop
   end if
end subroutine require

end module var_univariate_mod
