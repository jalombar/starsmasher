      implicit none
      integer nmax,nnmax,ntab,ngr,nav
      real*8 pi,eta2,alpha,beta,nax
      PARAMETER (NMAX=420000,NNMAX=128)
c      PARAMETER (NMAX=301000,NNMAX=128)
c      PARAMETER (NMAX=221000,NNMAX=1000)
c      PARAMETER (NMAX=180000,NNMAX=100)
      PARAMETER (NTAB=100000)                                           
      PARAMETER (PI=3.1415926535897932384626d0)
      COMMON/ARTVIS/ ETA2,ALPHA,BETA,NAV
      COMMON/GRAV/ NGR
      real*8 omega2,trelax,dt,gam,hmin,hmax
      integer nrelax,n,nnopt,nout,nit
      real*8 t,tf,dtout
      COMMON/RELAXP/ TRELAX,omega2,NRELAX
      COMMON/INTPAR/ DT,GAM,HMIN,HMAX,N,NNOPT
      COMMON/OUT/ T,TF,DTOUT,NOUT,NIT     
      real*8 X(NMAX),Y(NMAX),Z(NMAX),VX(NMAX),VY(NMAX),VZ(NMAX),  
     $      AM(NMAX),HP(NMAX),A(NMAX),RHO(NMAX),POR2(NMAX),
     $      GRPOT(nmax),wmeanmolecular(nmax),BB(NMAX),AA(NMAX)
      COMMON/PART/ X,Y,Z,VX,VY,VZ,  
     $      AM,HP,A,RHO,POR2,
     $      GRPOT,wmeanmolecular,BB,AA
      real*8 VXDOT(NMAX),VYDOT(NMAX),VZDOT(NMAX),ADOT(NMAX),       
     $      GX(NMAX),GY(NMAX),GZ(NMAX)                                  
      COMMON/DYN/ VXDOT,VYDOT,VZDOT,ADOT,       
     $      GX,GY,GZ                                  
      real*8 kelvin,gram,sec,cm,erg,boltz,crad,planck,crad2,sigma,arad,
     $      qconst
      PARAMETER(gram=1.d0,sec=1.d0,cm=1.d0,kelvin=1.d0)
      PARAMETER(erg=gram*cm**2/sec**2)
      PARAMETER(boltz=1.380658d-16*erg/kelvin)
      parameter(crad = 2.997924580d+10*cm/sec,  
     1      planck = 6.6260755d-27*gram*cm**2/sec)
      parameter( crad2 = crad**2,
     3     sigma = pi**2*boltz*(boltz*2.0d0*pi/planck)**3/60.0d0/crad2,
     4     arad  = 4.0d0*sigma/crad, qconst = 1.5d0*boltz/arad)
      real*8 MUNIT,RUNIT,PUNIT,gravconst
      parameter(gravconst = 6.67390d-08)
      common/units/munit,runit,punit
c      parameter(RUNIT=6.9599d10,MUNIT=1.9891d33)
c      parameter(PUNIT=gravconst*(MUNIT/RUNIT**2)**2)
      integer CC(NMAX)
      COMMON/SOFTENING/ CC
      real*8 wtab(ntab),dwtab(ntab),ctab,dwdhtab(ntab)
      common/wtabul/ wtab,dwtab,ctab,dwdhtab
