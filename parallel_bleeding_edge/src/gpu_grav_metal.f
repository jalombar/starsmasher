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

c     the number of multiprocessors of the GPU: 0, not known (the Metal library does not say it)
      subroutine get_gpu_sm_count(nsm,myrank_in)
      implicit none
      integer nsm,myrank_in
      nsm=0
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

c     ngravtree=2 (the tree gravity on the GPUs) needs the CUDA library: gravtree.f90 calls these, init.f asks
c     gpu_tree_available first
      subroutine gpu_tree_available(ok)
      integer ok
      ok=0
      return
      end

      subroutine gpu_tree_reserve(n, nnode, pa, pb, pneed, pnf,
     $     pni, pwork, pout)
      integer n,nnode
      stop 'ngravtree=2 needs the CUDA (NVIDIA GPU) build'
      end

      subroutine gpu_tree_launch(n, nnode, nwork, nkernel, theta2)
      integer n,nnode,nwork,nkernel
      stop 'ngravtree=2 needs the CUDA (NVIDIA GPU) build'
      end

      subroutine gpu_tree_wait(ierr)
      integer ierr
      stop 'ngravtree=2 needs the CUDA (NVIDIA GPU) build'
      end

c     SPH on the GPU (ngpusph): the CUDA library has it, this build does not (gpusph.f90 asks gpu_sph_available first)
      subroutine gpu_sph_available(ok)
      integer ok
      ok=0
      return
      end

      subroutine gpu_sph_setup(localrank,ngpuperrank,ngravprocs,ierr)
      integer localrank,ngpuperrank,ngravprocs,ierr
      ierr=1
      return
      end

      subroutine gpu_sph_tables(ntab,gtab,dgdhtab,wtab,dwdhtab,
     $     dphidhtab,ctab,ierr)
      integer ntab,ierr
      real*8 gtab(*),dgdhtab(*),wtab(*),dwdhtab(*),dphidhtab(*),ctab
      ierr=1
      return
      end

      subroutine gpu_sph_rhoh(n,nlow,nup,nmaxarr,x,y,z,am,u,hp,divv,
     $     rho,omega,om0,psi,wn,zeta,nn,first,sflag,list,listcap,ltot,
     $     dt,hfloor,nnopt,nselfg,skin,skiplist,work,ierr)
      integer n,nlow,nup,nmaxarr,nn(*),first(*),sflag(*),list(*)
      integer listcap,ltot,nnopt,nselfg,skiplist,ierr
      integer*8 work
      real*8 x(*),y(*),z(*),am(*),u(*),hp(*),divv(*),rho(*),omega(*)
      real*8 om0(*),psi(*),wn(*),zeta(*),dt,hfloor,skin
      ltot=0
      ierr=1
      return
      end

      subroutine gpu_sph_getlists(list,ierr)
      integer list(*),ierr
      ierr=1
      return
      end

      subroutine gpu_sph_fearly(n,vx,vy,vz,u,ierr)
      integer n,ierr
      real*8 vx(*),vy(*),vz(*),u(*)
      ierr=1
      return
      end

      subroutine gpu_sph_ftables(ntab,dwtab,dgtab,ierr)
      integer ntab,ierr
      real*8 dwtab(*),dgtab(*)
      ierr=1
      return
      end

      subroutine gpu_sph_forces(n,nmaxarr,nt,nloop,contig,tp,tiw,x,y,z,
     $     vx,vy,vz,am,u,hp,rho,por2,om,om0,psi,nn,first,list,listcap,
     $     act,blocks,alpha1,beta1,hfloor,gam,nav,nintvar,relaxnoheat,
     $     udotother,ngr,useb,prec,mvx,mvy,mvz,mgx,mgy,mgz,udot,uijmax,
     $     divv,nlarge,work,ierr)
      integer n,nmaxarr,nt,nloop,contig,tp(*),tiw(*),nn(*),first(*)
      integer list(*),listcap,act(*),blocks,nav,nintvar,relaxnoheat
      integer udotother,ngr,useb,prec,nlarge,ierr
      integer*8 work
      real*8 x(*),y(*),z(*),vx(*),vy(*),vz(*),am(*),u(*),hp(*),rho(*)
      real*8 por2(*),om(*),om0(*),psi(*),alpha1,beta1,hfloor,gam
      real*8 mvx(*),mvy(*),mvz(*),mgx(*),mgy(*),mgz(*),udot(*)
      real*8 uijmax(*),divv(*)
      nlarge=0
      ierr=1
      return
      end

c     the integration on the device (gpusph.f90, gi_*): the CUDA library has it, this build does not (gpusph.f90 never asks for it here: gpu_sph_available says 0)
      subroutine gpu_sph_int_flags(fl)
      integer fl
      fl=0
      return
      end

      subroutine gpu_sph_int_hash14(i,w,out)
      integer i
      integer*8 out
      real*8 w(*)
      out=0
      return
      end

      subroutine gpu_sph_int_stats(s)
      integer*8 s(*)
      integer k
      do k=1,16
         s(k)=0
      enddo
      return
      end

      subroutine gpu_sph_int_gind(g,n)
      integer g(*),n
      return
      end

      subroutine gpu_sph_int_go(n,nlow,nup,dt,dth,dthnew,skin,skiplist,
     $     hfloor,nnopt,nselfg,ngrav,vxo,vyo,vzo,uo,vxdot,vydot,vzdot,
     $     udot,uijmax,divv,status)
      integer n,nlow,nup,skiplist,nnopt,nselfg,ngrav,status
      real*8 dt,dth,dthnew,skin,hfloor
      real*8 vxo(*),vyo(*),vzo(*),uo(*),vxdot(*),vydot(*),vzdot(*)
      real*8 udot(*),uijmax(*),divv(*)
      status=1
      return
      end

      subroutine gpu_sph_int_check(n,nlow,nup,dt,skin,skiplist,hfloor,
     $     nnopt,nselfg,hosthash,x,y,z,vx,vy,vz,u,vxo,vyo,vzo,uo,hp,
     $     divv,am,adopt)
      integer n,nlow,nup,skiplist,nnopt,nselfg,adopt
      integer*8 hosthash
      real*8 dt,skin,hfloor
      real*8 x(*),y(*),z(*),vx(*),vy(*),vz(*),u(*),vxo(*),vyo(*)
      real*8 vzo(*),uo(*),hp(*),divv(*),am(*)
      adopt=0
      return
      end

      subroutine gpu_sph_int_drop
      return
      end

      subroutine gpu_sph_int_arm(armed,n,nlow,nup,dth,skin,skiplist,
     $     hfloor,nnopt,nselfg,ngrav,cn,dtforce,vxo,vyo,vzo,uo,vxdot,
     $     vydot,vzdot,udot,uijmax,divv,cb)
      integer armed,n,nlow,nup,skiplist,nnopt,nselfg,ngrav
      integer*8 cb
      real*8 dth,skin,hfloor,cn(*),dtforce
      real*8 vxo(*),vyo(*),vzo(*),uo(*),vxdot(*),vydot(*),vzdot(*)
      real*8 udot(*),uijmax(*),divv(*)
      return
      end

      subroutine gpu_sph_ndev(nd)
      integer nd
      nd=0
      return
      end
