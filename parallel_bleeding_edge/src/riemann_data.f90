module riemann_data
!     Arrays for the Riemann-solver dissipation, nav=4 (see riemann.f).
!     They are allocated only when nav=4, on first use.
  implicit none

!     gam1p(i): particle i's adiabatic exponent Gamma_1, set in pressure.f
  real*8, allocatable :: gam1p(:)
end module riemann_data
