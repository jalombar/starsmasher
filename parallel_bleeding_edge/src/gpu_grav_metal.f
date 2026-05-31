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
      implicit none
      integer ngravprocs_in,nprocs_in,ppn_in
      integer n_lower,n_upper,myrank,nprocs,qthreads,q
      integer rankspernode,nnodes
      common/nlimits/n_lower,n_upper,nprocs,myrank,qthreads,q

      rankspernode=max(ppn_in,1)
      nnodes=max((nprocs_in+rankspernode-1)/rankspernode,1)
      if(ngravprocs_in.gt.nnodes) then
         if(myrank.eq.0) then
            write(69,*)'apple silicon metal backend uses one gpu per node'
            write(69,*)'resetting ngravprocs from',ngravprocs_in,
     $           'to',nnodes
         endif
         ngravprocs_in=nnodes
      endif
      return
      end
