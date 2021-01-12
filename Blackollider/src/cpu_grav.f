      subroutine get_gravity_using_cpus
      include 'starsmasher.h'
      integer i,j,jlower,jupper
      real*8 dr1,dr2,dr3
      real*8 r2,rinv,rinv2,mrinv1,mrinv3,gacc,gpot
      real*8 invq1,qq1,q21,q31,invq21,invq31,acc1,
     $     pot1,g1,mj11,mj21
      real*8 invq2,qq2,q22,q32,invq22,invq32,acc2,
     $     pot2,g2,mj12,mj22
      real*8 dgx, dgy, dgz
      real*8 ami, amj, massratio
      real*8 twohpi, twohpj, fourhpi2, fourhpj2

c      write(69,*) 'particle',ntot,'has u,x,y,z,m,h=',
c     $     u(ntot),x(ntot),y(ntot),z(ntot),am(ntot),hp(ntot)

      do i=1,ntot
         gx(i)=0.d0
         gy(i)=0.d0
         gz(i)=0.d0
         grpot(i)=0.d0
      enddo

      if(nselfgravity.eq.0)then
         jlower=ntot
         if(u(jlower).ne.0.d0)then
            write(69,*)'last particle not a point particle?'
            stop
         endif
      endif
      jupper=ntot

      do i=ngrav_lower,ngrav_upper

         if(i.eq.ntot) jlower=1 ! so any last point particle interacts with everything

         if(nselfgravity.eq.1) jlower=i

         ami=am(i)
         twohpi=2*hp(i)
         fourhpi2=twohpi*twohpi

         if(nkernel.eq.0) then
c     nkernel=0 is the cubic spline kernel

            do j=jlower,jupper
               dr1=x(j)-x(i)
               dr2=y(j)-y(i)
               dr3=z(j)-z(i)
               r2=dr1**2+dr2**2+dr3**2
               rinv  = 1/sqrt(r2)
               rinv2 = rinv   * rinv
               amj=am(j)
               twohpj=2*hp(j)
               fourhpj2=twohpj*twohpj
               massratio=ami/amj
               mrinv1 = rinv   * amj
               mrinv3 = rinv2  * mrinv1
               if(r2.ge.max(fourhpi2,fourhpj2))then
                  dgx=mrinv3 * dr1
                  dgy=mrinv3 * dr2
                  dgz=mrinv3 * dr3
                  gx(i) = gx(i) + dgx
                  gx(j) = gx(j) - dgx*massratio
                  gy(i) = gy(i) + dgy
                  gy(j) = gy(j) - dgy*massratio
                  gz(i) = gz(i) + dgz
                  gz(j) = gz(j) - dgz*massratio
                  grpot(i) = grpot(i) -mrinv1
                  grpot(j) = grpot(j) -mrinv1*massratio
               else
                  if(r2.gt.0)then
                     invq1= rinv * twohpi
                     invq2= rinv * twohpj
                     qq1= 1/invq1
                     qq2= 1/invq2
                  else
                     invq1= 0.d0
                     invq2= 0.d0
                     qq1= 0.d0
                     qq2= 0.d0
                  endif
                  q21=qq1*qq1
                  q22=qq2*qq2
                  q31=qq1*q21
                  q32=qq2*q22
                  invq21=invq1*invq1
                  invq22=invq2*invq2
                  invq31=invq1*invq21
                  invq32=invq2*invq22
                  if(qq1.lt.0.5d0) then
                     acc1=
     $                    10.666666666667d0+q21*(32.0d0*qq1-38.4d0)
                     pot1=     -2.8d0+q21*(5.333333333333d0+q21*
     $                    (6.4d0*qq1-9.6d0))
                  else
                     acc1=
     $                    21.333333333333d0-48.0d0*qq1+38.4d0*q21
     $                    -10.666666666667d0*q31-0.066666666667d0*invq31
                     pot1=    
     $                    -3.2d0+0.066666666667d0*invq1+q21*
     $                    (10.666666666667d0+qq1*(-16.0d0+qq1*
     $                    (9.6d0-2.133333333333d0*qq1)))
                  endif
                  if(qq2.lt.0.5d0) then
                     acc2=
     $                    10.666666666667d0+q22*(32.0d0*qq2-38.4d0)
                     pot2=     -2.8d0+q22*(5.333333333333d0+q22*
     $                    (6.4d0*qq2-9.6d0))
                  else
                     acc2=
     $                    21.333333333333d0-48.0d0*qq2+38.4d0*q22
     $                    -10.666666666667d0*q32-0.066666666667d0*invq32
                     pot2=     
     $                    -3.2d0+0.066666666667d0*invq2+q22*
     $                    (10.666666666667d0+qq2*(-16.0d0+qq2*
     $                    (9.6d0-2.133333333333d0*qq2)))
                  endif
                  mj11=amj/twohpi
                  mj12=amj/twohpj
                  mj21=mj11/fourhpi2
                  mj22=mj12/fourhpj2
                  if(r2.le.fourhpi2)then
                     g1=1
                  else
                     g1=0
                  endif
                  if(r2.le.fourhpj2)then
                     g2=1
                  else
                     g2=0
                  endif
                  if(r2.gt.0)then
                     gacc=0.5d0*(g1*mj21*acc1+(1.0d0-g1)*mrinv3+
     $                    g2*mj22*acc2+(1.0d0-g2)*mrinv3)
                     gpot=0.5d0*(g1*mj11*pot1+(g1-1.0d0)*mrinv1+
     $                    g2*mj12*pot2+(g2-1.0d0)*mrinv1)
                  else
                     gacc=0.d0
                     if(u(i).ne.0)then
                        gpot=-0.7d0*(mj11 + mj12) ! 0.7 instead of 1.4 so that we don't double count below when i=j
                     else
                        gpot=0.d0
                     endif
                  endif
                  dgx=gacc * dr1
                  dgy=gacc * dr2
                  dgz=gacc * dr3
                  gx(i) = gx(i) + dgx
                  gx(j) = gx(j) - dgx*massratio
                  gy(i) = gy(i) + dgy
                  gy(j) = gy(j) - dgy*massratio
                  gz(i) = gz(i) + dgz
                  gz(j) = gz(j) - dgz*massratio
                  grpot(i) = grpot(i) + gpot
                  grpot(j) = grpot(j) + gpot*massratio
               endif
            enddo
         else if (nkernel.eq.1) then
c     nkernel=1 is the Wendland W_{3,3} kernel

            do j=jlower,jupper
               dr1=x(j)-x(i)
               dr2=y(j)-y(i)
               dr3=z(j)-z(i)
               r2=dr1**2+dr2**2+dr3**2
               rinv  = 1/sqrt(r2)
               rinv2 = rinv   * rinv
               amj=am(j)
               twohpj=2*hp(j)
               fourhpj2=twohpj*twohpj
               massratio=ami/amj
               mrinv1 = rinv   * amj
               mrinv3 = rinv2  * mrinv1
               if(r2.ge.max(fourhpi2,fourhpj2))then
                  dgx=mrinv3 * dr1
                  dgy=mrinv3 * dr2
                  dgz=mrinv3 * dr3
                  gx(i) = gx(i) + dgx
                  gx(j) = gx(j) - dgx*massratio
                  gy(i) = gy(i) + dgy
                  gy(j) = gy(j) - dgy*massratio
                  gz(i) = gz(i) + dgz
                  gz(j) = gz(j) - dgz*massratio
                  grpot(i) = grpot(i) -mrinv1
                  grpot(j) = grpot(j) -mrinv1*massratio
               else
                  if(r2.gt.0)then
                     invq1= rinv * twohpi
                     invq2= rinv * twohpj
                     qq1= 1/invq1
                     qq2= 1/invq2
                  else
                     invq1= 0.d0
                     invq2= 0.d0
                     qq1= 0.d0
                     qq2= 0.d0
                  endif
!     Relate to the CUDA code grav_force_direct.cu...
!     qq1 = q.x,  qq2 = q.y
!     q21 = q2.x, q22 = q2.y
                  q21=qq1*qq1 
                  q22=qq2*qq2
                  q31=qq1*q21
                  q32=qq2*q22
                  invq21=invq1*invq1
                  invq22=invq2*invq2
                  invq31=invq1*invq21
                  invq32=invq2*invq22
                  acc1=28.4375d0 + q21 * ( (-187.6875d0) + q21 * (
     $                 804.375d0 + q21 * ( (-4379.375d0) + qq1 * (
     $                 9009.d0 + qq1 * ( (-8957.8125d0) + qq1 * (5005.d0
     $                 + qq1 * ( (-1515.9375d0) + 195.d0*qq1)))))))
                  acc2=28.4375d0 + q22 * ( (-187.6875d0) + q22 * (
     $                 804.375d0 + q22 * ( (-4379.375d0) + qq2 * (
     $                 9009.d0 + qq2 * ( (-8957.8125d0) + qq2 * (5005.d0
     $                 + qq2 * ( (-1515.9375d0) + 195.d0*qq2)))))))
                  mj11=amj/twohpi
                  mj12=amj/twohpj
                  mj21=mj11/fourhpi2
                  mj22=mj12/fourhpj2
                  if(r2.le.fourhpi2)then
                     g1=1
                     pot1=(-3.828125d0) + q21 * (14.21875d0 + q21 * (
     $                    (-46.921875d0) + q21 * (134.0625d0 + q21 * (
     $                    (-547.421875d0) + qq1 * (1001.d0 + qq1 * (
     $                    (-895.78125d0) + qq1 * (455.d0 + qq1 * (
     $                    (-126.328125d0) + 15.d0*qq1))))))))
                  else
                     g1=0
                     pot1=0
                  endif
                  if(r2.le.fourhpj2)then
                     g2=1
                     pot2=(-3.828125d0) + q22 * (14.21875d0 + q22 * (
     $                    (-46.921875d0) + q22 * (134.0625d0 + q22 * (
     $                    (-547.421875d0) + qq2 * (1001.d0 + qq2 * (
     $                    (-895.78125d0) + qq2 * (455.d0 + qq2 * (
     $                    (-126.328125d0) + 15.d0*qq2))))))))
                  else
                     g2=0
                     pot2=0
                  endif
                  if(r2.gt.0)then
                     gacc=0.5d0*(g1*mj21*acc1+(1.0d0-g1)*mrinv3+
     $                    g2*mj22*acc2+(1.0d0-g2)*mrinv3)
                     gpot=0.5d0*(mj11*pot1+(g1-1.0d0)*mrinv1+
     $                    mj12*pot2+(g2-1.0d0)*mrinv1)
                  else
                     gacc=0.d0
                     if(u(i).ne.0)then
                        gpot=-0.5d0*1.9140625d0*(mj11 + mj12) ! 0.5 factor is so that we don't double count below when i=j
                     else
                        gpot=0.d0
                     endif
                  endif
                  dgx=gacc * dr1
                  dgy=gacc * dr2
                  dgz=gacc * dr3
                  gx(i) = gx(i) + dgx
                  gx(j) = gx(j) - dgx*massratio
                  gy(i) = gy(i) + dgy
                  gy(j) = gy(j) - dgy*massratio
                  gz(i) = gz(i) + dgz
                  gz(j) = gz(j) - dgz*massratio
                  grpot(i) = grpot(i) + gpot
                  grpot(j) = grpot(j) + gpot*massratio
               endif
            enddo
         else
c     nkernel=2 is the Wendland C4 kernel

            do j=jlower,jupper
               dr1=x(j)-x(i)
               dr2=y(j)-y(i)
               dr3=z(j)-z(i)
               r2=dr1**2+dr2**2+dr3**2
               rinv  = 1/sqrt(r2)
               rinv2 = rinv   * rinv
               amj=am(j)
               twohpj=2*hp(j)
               fourhpj2=twohpj*twohpj
               massratio=ami/amj
               mrinv1 = rinv   * amj
               mrinv3 = rinv2  * mrinv1
               if(r2.ge.max(fourhpi2,fourhpj2))then
                  dgx=mrinv3 * dr1
                  dgy=mrinv3 * dr2
                  dgz=mrinv3 * dr3
                  gx(i) = gx(i) + dgx
                  gx(j) = gx(j) - dgx*massratio
                  gy(i) = gy(i) + dgy
                  gy(j) = gy(j) - dgy*massratio
                  gz(i) = gz(i) + dgz
                  gz(j) = gz(j) - dgz*massratio
                  grpot(i) = grpot(i) -mrinv1
                  grpot(j) = grpot(j) -mrinv1*massratio
               else
                  if(r2.gt.0)then
                     invq1= rinv * twohpi
                     invq2= rinv * twohpj
                     qq1= 1/invq1
                     qq2= 1/invq2
                  else
                     invq1= 0.d0
                     invq2= 0.d0
                     qq1= 0.d0
                     qq2= 0.d0
                  endif
!     Relate to the CUDA code grav_force_direct.cu...
!     qq1 = q.x,  qq2 = q.y
!     q21 = q2.x, q22 = q2.y
                  q21=qq1*qq1 
                  q22=qq2*qq2
                  q31=qq1*q21
                  q32=qq2*q22
                  invq21=invq1*invq1
                  invq22=invq2*invq2
                  invq31=invq1*invq21
                  invq32=invq2*invq22
                  acc1=20.625d0 + q21 * ( (-115.5d0) + q21 * (618.75d0 +
     $                 qq1 * ( (-1155.d0) + qq1 * (962.5d0 + qq1 * (
     $                 (-396.d0) + 65.625d0 * qq1)))))
                  acc2=20.625d0 + q22 * ( (-115.5d0) + q22 * (618.75d0 +
     $                 qq2 * ( (-1155.d0) + qq2 * (962.5d0 + qq2 * (
     $                 (-396.d0) + 65.625d0 * qq2)))))
                  mj11=amj/twohpi
                  mj12=amj/twohpj
                  mj21=mj11/fourhpi2
                  mj22=mj12/fourhpj2
                  if(r2.le.fourhpi2)then
                     g1=1
                     pot1=(-3.4375d0) + q21 * (10.3125d0 + q21 * (
     $                    (-28.875d0)  + q21 * (103.125d0 + qq1 *(
     $                    (-165.d0) + qq1 * (120.3125d0 + qq1 * ( (-44.d0)
     $                    + 6.5625d0*qq1))))))
                  else
                     g1=0
                     pot1=0
                  endif
                  if(r2.le.fourhpj2)then
                     g2=1
                     pot2=(-3.4375d0) + q22 * (10.3125d0 + q22 * (
     $                    (-28.875d0)  + q22 * (103.125d0 + qq2 *(
     $                    (-165.d0) + qq2 * (120.3125d0 + qq2 * ( (-44.d0)
     $                    + 6.5625d0*qq2))))))
                  else
                     g2=0
                     pot2=0
                  endif
                  if(r2.gt.0)then
                     gacc=0.5d0*(g1*mj21*acc1+(1.0d0-g1)*mrinv3+
     $                    g2*mj22*acc2+(1.0d0-g2)*mrinv3)
                     gpot=0.5d0*(mj11*pot1+(g1-1.0d0)*mrinv1+
     $                    mj12*pot2+(g2-1.0d0)*mrinv1)
                  else
                     gacc=0.d0
                     if(u(i).ne.0)then
                        gpot=-0.5d0*1.71875d0*(mj11 + mj12) ! 0.5 factor is so that we don't double count below when i=j
                     else
                        gpot=0.d0
                     endif
                  endif
                  dgx=gacc * dr1
                  dgy=gacc * dr2
                  dgz=gacc * dr3
                  gx(i) = gx(i) + dgx
                  gx(j) = gx(j) - dgx*massratio
                  gy(i) = gy(i) + dgy
                  gy(j) = gy(j) - dgy*massratio
                  gz(i) = gz(i) + dgz
                  gz(j) = gz(j) - dgz*massratio
                  grpot(i) = grpot(i) + gpot
                  grpot(j) = grpot(j) + gpot*massratio
               endif
            enddo
         endif
      enddo

      return
      end

      subroutine set_nusegpus
      implicit none
      integer nintvar,neos,nusegpus,nselfgravity,ncooling
      common/integration/nintvar,neos,nusegpus,nselfgravity,ncooling
      nusegpus=0
      return
      end

      subroutine firsthalf_grav_forces(ntot, ngrav_lower, mygravlength,
     $     x, y, z, am, range, q, nkernel)
      integer ntot,ngrav_lower,mygravlength,q,nkernel
      real*8 x(ntot),y(ntot),z(ntot),am(ntot),range(ntot)
      return
      end

      subroutine lasthalf_grav_forces(ntot, gx, gy, gz, grpot)
      integer ntot
      real*8 gx(ntot),gy(ntot),gz(ntot),grpot(ntot)
      return
      end

      subroutine gpu_init_dev(i,theta_angle)
      integer i
      real*8 theta_angle
      return
      end
