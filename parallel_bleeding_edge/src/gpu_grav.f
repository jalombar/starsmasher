      subroutine set_nusegpus
      implicit none
      integer nintvar,neos,nusegpus,nselfgravity,ncooling,nkernel
      common/integration/nintvar,neos,nusegpus,nselfgravity,ncooling,
     $     nkernel
      nusegpus=1
      nselfgravity=1
      return
      end

      subroutine get_gpu_count(actual_gpu_count,myrank_in)
      implicit none
      integer actual_gpu_count,myrank_in,ios
      character*20 gpucfilename
      character*80 command

      actual_gpu_count=0
      write(gpucfilename,'(A,I4.4)')'temp_gpu_count_',myrank_in
      command='nvidia-smi -L 2>/dev/null | wc -l > '
     $     //trim(gpucfilename)
      call system(trim(command))
      open(unit=400+myrank_in,file=gpucfilename,status='old',
     $     action='read',iostat=ios)
      if(ios.eq.0) then
         read(400+myrank_in,*,iostat=ios)actual_gpu_count
         close(400+myrank_in)
         call system('rm -f '//trim(gpucfilename))
      endif
      if(ios.ne.0) actual_gpu_count=0
      return
      end

c     the number of multiprocessors of the GPU of this process (the smallest of those it sees), which the automatic choice of ngpusph
c     takes as a coarse hint of how fast the GPU is (init.f, auto_gpuclass); 0 if it cannot be read.  The C side is in
c     SPHgrav_lib2/grav_force_direct.cu; the CPU build and the Metal one say 0
      subroutine get_gpu_sm_count(nsm,myrank_in)
      implicit none
      integer nsm,myrank_in
      nsm=0
      call gpu_sm_count(nsm)
      return
      end

      subroutine get_gravity_using_cpus
      return
      end
