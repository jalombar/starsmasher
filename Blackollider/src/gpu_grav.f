      subroutine set_nusegpus
      implicit none
      integer nintvar,neos,nusegpus,nselfgravity
      common/integration/nintvar,neos,nusegpus,nselfgravity
      nusegpus=1
      nselfgravity=1
      return
      end

      subroutine get_gravity_using_cpus
      return
      end
