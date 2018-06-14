      implicit none

      integer nmax,nnmax
      parameter(nmax=1600000,nnmax=170)
      integer n_lower,n_upper,myrank,nprocs
      integer qthreads,q,gflag
      real*8 mbh,reat,starmass,starradius
      common/massblackhole/mbh,reat
      common/starsize/starmass,starradius
      common/nlimits/n_lower,n_upper,nprocs,myrank,qthreads,q
      common/flag/gflag
      integer ngrav_lower,ngrav_upper,ngravprocs
      common/ngravlimits/ngrav_lower,ngrav_upper,ngravprocs
      real*8 vxdotsm(nmax),vydotsm(nmax),vzdotsm(nmax)
      integer n,nnopt,ntot,nout,nit,nitpot,ngr,nrelax,nav,nintvar
      common/mine/vxdotsm,vydotsm,vzdotsm
      integer ntab
      parameter(ntab=100000)                                           
      real*8 pi,alpha,beta
      real*8 trelax,sep0,t,tf,dtout,dt,dth,hco,hfloor,tscanon,sepfinal
      real*8 dgdhtab(ntab), dgtab(ntab), gtab(ntab)
      real*8 wtab(ntab),dwtab(ntab),ctab,dwdhtab(ntab),dphidhtab(ntab)
      real*8 x(nmax),y(nmax),z(nmax),vx(nmax),vy(nmax),vz(nmax)  
      real*8 am(nmax),hp(nmax),u(nmax),rho(nmax),por2(nmax)
      real*8 grpot(nmax),meanmolecular(nmax)
      real*8 zeta(nmax),bonet_omega(nmax), bonet_0mega(nmax) 
      real*8 bonet_psi(nmax), bonet_wn(nmax), max_vsig(nmax)
      real*8 vxdot(nmax),vydot(nmax),vzdot(nmax),udot(nmax)
      real*8 gx(nmax),gy(nmax),gz(nmax),tthermal(nmax),ueq(nmax)
      integer nn(nmax),list(nnmax*nmax),first(nmax), ngb_j(nmax)
      real*8 rp,vinf2,equalmass,treloff,tresplintmuoff,bimpact 
      real time0,time1,time2,time3,seconds
      parameter (pi=3.1415926535897932384626d0)
      real*8 arad,boltz,erg,cm,gram,sec,kelvin,qconst
      real*8 crad2,sigma,crad,planck
      parameter(gram=1.d0,sec=1.d0,cm=1.d0,kelvin=1.d0)
      parameter(erg=gram*cm**2/sec**2)
      parameter(boltz=1.380658d-16*erg/kelvin)
      parameter(crad=2.997924580d+10*cm/sec)
      parameter(planck=6.6260755d-27*gram*cm**2/sec)
      parameter(crad2=crad**2)
      parameter(sigma=pi**2*boltz*(boltz*2d0*pi/planck)**3/60d0/crad2)
      parameter(arad=4.0d0*sigma/crad,qconst=1.5d0*boltz/arad)
      real*8 munit,runit,punit,gravconst,redge
      parameter(gravconst = 6.67390d-08)
      common/units/munit,runit,punit
      common/artvis/ alpha,beta,nav
      integer computeexclusivemode,ppn
      common/grav/ ngr,computeexclusivemode
      real*8 omega2,omega_spin
      common/relaxp/ trelax,bimpact,vinf2,equalmass,treloff,tresplintmuoff,omega2,omega_spin,nrelax,nitpot,stellarevolutioncodetype
      common/intpar/ dt,n,nnopt,hco,hfloor,dth,ntot
      common/out/ t,tf,dtout,nout,nit     
      common/wtabul/dgdhtab,dgtab,gtab,wtab,dwtab,ctab,dwdhtab,dphidhtab
      common/part/ x,y,z,vx,vy,vz,am,hp,u,rho,por2,grpot,zeta,bonet_omega,bonet_0mega, bonet_psi, bonet_wn, max_vsig
      common/dyn/ vxdot,vydot,vzdot,udot,gx,gy,gz,tthermal,ueq
      integer kdm
      common/neigh/ nn,list,first, ngb_j
      common/radpress/ meanmolecular,redge
      real*8 cn1,cn2,cn3,cn4,cn5,cn6,cn7
      real*8 e0, semimajoraxis
      integer neos,nusegpus,nselfgravity,ncooling,nkernel
      real*8 gam,teq,tjumpahead
      character*255 startfile1,startfile2,eosfile,opacityfile,profilefile
      logical throwaway
      integer stellarevolutioncodetype	
      common/inputfilenames/startfile1,startfile2,eosfile,opacityfile,profilefile
      common/courantnumbers/ cn1,cn2,cn3,cn4,cn5,cn6,cn7
      common/integration/nintvar,neos,nusegpus,nselfgravity,ncooling,nkernel
      parameter(kdm=5000)
      integer ntypes,cc(nmax)
      parameter(ntypes=32)
      common/softening/cc
      integer n1
      common/binaryp/ sep0,rp, tscanon,sepfinal,n1,throwaway
      logical gonedynamic
      real*8 omeg
      common/rotating/omeg,gonedynamic
      integer nprocsmax
      parameter(nprocsmax=128)
      integer displs(nprocsmax),recvcounts(nprocsmax)
      common/displacements/displs,recvcounts,ppn
      integer ngravprocsmax
      parameter(ngravprocsmax=32)
      integer gravdispls(ngravprocsmax),gravrecvcounts(ngravprocsmax)
      common/gravdisplacements/gravdispls,gravrecvcounts
      common/adiabaticindex/ gam
      integer maxnumx
      parameter(maxnumx=10)
