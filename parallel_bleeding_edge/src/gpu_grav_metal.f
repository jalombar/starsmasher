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
