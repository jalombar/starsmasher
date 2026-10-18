module riemann_data
!     Arrays for the Riemann-solver dissipation, nav=4 (see riemann.f).
!     They are allocated only when nav=4, on first use.
  implicit none

!     gam1p(i): particle i's adiabatic exponent Gamma_1, set in pressure.f
  real*8, allocatable :: gam1p(:)
!     rgrad(:,i): SPH gradients at particle i of rho (1:3), p (4:6) and
!     the velocity, d v_b/d x_a in 6+3*(b-1)+a
  real*8, allocatable :: rgrad(:,:)
end module riemann_data
