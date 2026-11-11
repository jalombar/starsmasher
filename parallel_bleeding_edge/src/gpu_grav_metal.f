      subroutine set_nusegpus
      implicit none
      integer nintvar,neos,nusegpus,nselfgravity,ncooling,nkernel
      common/integration/nintvar,neos,nusegpus,nselfgravity,ncooling,
     $     nkernel
      nusegpus=1
      nselfgravity=1
      return
      end

      subroutine get_gravity_using_cpus
      return
      end

      subroutine get_gpu_count(actual_gpu_count,myrank_in)
      implicit none
      integer actual_gpu_count,myrank_in
      actual_gpu_count=1
      return
      end

      subroutine gpu_set_ndev(i,n,np)
      integer i,n,np
      return
      end

      subroutine firsthalf_grav_forces_far(ntot, ngrav_lower, mygravlength,
     $     x, y, z, am, range, q, nkernel, ok)
      integer ntot,ngrav_lower,mygravlength,q,nkernel,ok
      real*8 x(ntot),y(ntot),z(ntot),am(ntot),range(ntot)
      ok=0
      return
      end

      subroutine firsthalf_grav_forces_near(ntot, x, y, z, am, range,
     $     nkernel)
      integer ntot,nkernel
      real*8 x(ntot),y(ntot),z(ntot),am(ntot),range(ntot)
      return
      end

      subroutine gravity_far_abort
      return
      end
