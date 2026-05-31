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

      subroutine adjust_ngravprocs_for_backend(ngravprocs_in,nprocs_in,
     $     ppn_in)
      integer ngravprocs_in,nprocs_in,ppn_in
      return
      end
