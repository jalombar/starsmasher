!***********************************************************************
!     The softened pair interaction of the gravity, for the particles whose
!     separation is inside the softening of one of them (r^2 below the
!     larger of (2h_i)^2 and (2h_j)^2).  Outside it the interaction is
!     Newtonian, G m_j/r^2, and the callers do that themselves.
!
!     This is the near-pair formula of get_gravity_using_cpus (cpu_grav.f)
!     and of the kernels of the GPU library (grav_force_direct.cu) in the
!     gather form: the acceleration of the particle i from j is
!     gacc*(x_j - x_i) and the potential of i gets gpot.  The same
!     expressions in the same order as cpu_grav.f for nkernel = 0 (cubic
!     spline), 1 (Wendland 3,3) and 2 (Wendland C4), and the same rules for
!     a pair at zero distance: no force, and the potential of the particle
!     itself, which a point particle (flag 0) does not have.
!***********************************************************************
module gravpair
  implicit none
  private
  public :: grav_soft
  integer, parameter :: dp = kind(1.d0)

contains

!     nk: nkernel.  r2: squared distance of i and j.  twohpi, twohpj: 2h of the two.  amj: mass of j.
!     flagj: 1 for a gas particle, 0 for a point particle.  selffac: 1 if j is i itself, 1/2 for another particle at zero distance
  pure subroutine grav_soft(nk, r2, twohpi, twohpj, amj, flagj, selffac, gacc, gpot)
    integer, intent(in) :: nk
    real(dp), intent(in) :: r2, twohpi, twohpj, amj, flagj, selffac
    real(dp), intent(out) :: gacc, gpot
    real(dp) :: rinv, rinv2, mrinv1, mrinv3, invq1, invq2, qq1, qq2, q21, q22, q31, q32, invq21, invq22, invq31, invq32
    real(dp) :: acc1, acc2, pot1, pot2, g1, g2, mj11, mj12, mj21, mj22, fourhpi2, fourhpj2, cself
    fourhpi2 = twohpi*twohpi
    fourhpj2 = twohpj*twohpj
    mj11 = amj/twohpi
    mj12 = amj/twohpj
    mj21 = mj11/fourhpi2
    mj22 = mj12/fourhpj2
    if (r2 <= 0.d0) then
       gacc = 0.d0
       select case (nk)
       case (0)
          cself = 1.4d0
       case (1)
          cself = 1.9140625d0
       case default
          cself = 1.71875d0
       end select
       gpot = -cself*selffac*flagj*(mj11 + mj12)
       return
    end if
    rinv = 1/sqrt(r2)
    rinv2 = rinv*rinv
    mrinv1 = rinv*amj
    mrinv3 = rinv2*mrinv1
    invq1 = rinv*twohpi
    invq2 = rinv*twohpj
    qq1 = 1/invq1
    qq2 = 1/invq2
    q21 = qq1*qq1
    q22 = qq2*qq2
    q31 = qq1*q21
    q32 = qq2*q22
    invq21 = invq1*invq1
    invq22 = invq2*invq2
    invq31 = invq1*invq21
    invq32 = invq2*invq22
    if (r2 <= fourhpi2) then
       g1 = 1.d0
    else
       g1 = 0.d0
    end if
    if (r2 <= fourhpj2) then
       g2 = 1.d0
    else
       g2 = 0.d0
    end if
    select case (nk)
    case (0)
       ! the cubic spline
       if (qq1 < 0.5d0) then
          acc1 = 10.666666666667d0 + q21*(32.0d0*qq1 - 38.4d0)
          pot1 = -2.8d0 + q21*(5.333333333333d0 + q21*(6.4d0*qq1 - 9.6d0))
       else
          acc1 = 21.333333333333d0 - 48.0d0*qq1 + 38.4d0*q21 - 10.666666666667d0*q31 - 0.066666666667d0*invq31
          pot1 = -3.2d0 + 0.066666666667d0*invq1 + q21*(10.666666666667d0 + qq1*(-16.0d0 + qq1*(9.6d0 - 2.133333333333d0*qq1)))
       end if
       if (qq2 < 0.5d0) then
          acc2 = 10.666666666667d0 + q22*(32.0d0*qq2 - 38.4d0)
          pot2 = -2.8d0 + q22*(5.333333333333d0 + q22*(6.4d0*qq2 - 9.6d0))
       else
          acc2 = 21.333333333333d0 - 48.0d0*qq2 + 38.4d0*q22 - 10.666666666667d0*q32 - 0.066666666667d0*invq32
          pot2 = -3.2d0 + 0.066666666667d0*invq2 + q22*(10.666666666667d0 + qq2*(-16.0d0 + qq2*(9.6d0 - 2.133333333333d0*qq2)))
       end if
       gacc = 0.5d0*(g1*mj21*acc1 + (1.0d0 - g1)*mrinv3 + g2*mj22*acc2 + (1.0d0 - g2)*mrinv3)
       gpot = 0.5d0*(g1*mj11*pot1 + (g1 - 1.0d0)*mrinv1 + g2*mj12*pot2 + (g2 - 1.0d0)*mrinv1)
    case (1)
       ! the Wendland W_{3,3} kernel
       acc1 = 28.4375d0 + q21*((-187.6875d0) + q21*(804.375d0 + q21*((-4379.375d0) + qq1*(9009.d0 + qq1*((-8957.8125d0) &
            + qq1*(5005.d0 + qq1*((-1515.9375d0) + 195.d0*qq1)))))))
       acc2 = 28.4375d0 + q22*((-187.6875d0) + q22*(804.375d0 + q22*((-4379.375d0) + qq2*(9009.d0 + qq2*((-8957.8125d0) &
            + qq2*(5005.d0 + qq2*((-1515.9375d0) + 195.d0*qq2)))))))
       if (r2 <= fourhpi2) then
          pot1 = (-3.828125d0) + q21*(14.21875d0 + q21*((-46.921875d0) + q21*(134.0625d0 + q21*((-547.421875d0) &
               + qq1*(1001.d0 + qq1*((-895.78125d0) + qq1*(455.d0 + qq1*((-126.328125d0) + 15.d0*qq1))))))))
       else
          pot1 = 0.d0
       end if
       if (r2 <= fourhpj2) then
          pot2 = (-3.828125d0) + q22*(14.21875d0 + q22*((-46.921875d0) + q22*(134.0625d0 + q22*((-547.421875d0) &
               + qq2*(1001.d0 + qq2*((-895.78125d0) + qq2*(455.d0 + qq2*((-126.328125d0) + 15.d0*qq2))))))))
       else
          pot2 = 0.d0
       end if
       gacc = 0.5d0*(g1*mj21*acc1 + (1.0d0 - g1)*mrinv3 + g2*mj22*acc2 + (1.0d0 - g2)*mrinv3)
       gpot = 0.5d0*(mj11*pot1 + (g1 - 1.0d0)*mrinv1 + mj12*pot2 + (g2 - 1.0d0)*mrinv1)
    case default
       ! the Wendland C4 kernel
       acc1 = 20.625d0 + q21*((-115.5d0) + q21*(618.75d0 + qq1*((-1155.d0) + qq1*(962.5d0 + qq1*((-396.d0) + 65.625d0*qq1)))))
       acc2 = 20.625d0 + q22*((-115.5d0) + q22*(618.75d0 + qq2*((-1155.d0) + qq2*(962.5d0 + qq2*((-396.d0) + 65.625d0*qq2)))))
       if (r2 <= fourhpi2) then
          pot1 = (-3.4375d0) + q21*(10.3125d0 + q21*((-28.875d0) + q21*(103.125d0 + qq1*((-165.d0) &
               + qq1*(120.3125d0 + qq1*((-44.d0) + 6.5625d0*qq1))))))
       else
          pot1 = 0.d0
       end if
       if (r2 <= fourhpj2) then
          pot2 = (-3.4375d0) + q22*(10.3125d0 + q22*((-28.875d0) + q22*(103.125d0 + qq2*((-165.d0) &
               + qq2*(120.3125d0 + qq2*((-44.d0) + 6.5625d0*qq2))))))
       else
          pot2 = 0.d0
       end if
       gacc = 0.5d0*(g1*mj21*acc1 + (1.0d0 - g1)*mrinv3 + g2*mj22*acc2 + (1.0d0 - g2)*mrinv3)
       gpot = 0.5d0*(mj11*pot1 + (g1 - 1.0d0)*mrinv1 + mj12*pot2 + (g2 - 1.0d0)*mrinv1)
    end select
  end subroutine grav_soft

end module gravpair
