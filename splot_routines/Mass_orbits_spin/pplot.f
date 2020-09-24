      INCLUDE 'spha.h'                                          
      CHARACTER*12 outfilename
      integer outfilenum,iform,iout,ni1,ni2,ni3,innit,nnit
      logical foundlastoutfile
      real*8 anglexdeg,angleydeg,anglezdeg
      common/azimuth/anglexdeg,angleydeg,anglezdeg
      integer ipp
      common/ppflag/ipp
      integer numfilters,numfiltersmax
      parameter (numfiltersmax=15)
      character*2 filtername(numfiltersmax)
      real*8 wavelength(numfiltersmax),absoluteflux(numfiltersmax)
      common /splotfilters/ wavelength,absoluteflux,numfilters
      common /filternames/ filtername
      real*8 tdesired, tnow, rl1last, tlast, rl1, rl2last, rl2
      character*255 cbeg_ferg
      integer i_mol,itype,kalex,kalex_avail
      real*8 starmass,starradius
      integer i
      real*8 divv(nmax)
      common/commdivv/divv
      real*8 sep0,
     $      bimpact,e0,semimajoraxis,vinf2,equalmass,treloff,
     $      tresplintmuoff,
     $      tscanon,sepfinal,
     $      mbh,
     $      cn1,cn2,cn3,cn4,cn5,cn6,cn7,
     $      omega_spin,reat
      integer nitpot,nintvar,ngravprocs,
     $      qthreads,Gflag,
     $      computeexclusivemode,
     $      ppn,neos,nselfgravity,nusegpus,ncooling
      real*8 teq, tjumpahead
      character*255 startfile1,startfile2,eosfile,opacityfile
      namelist/input/ tf,dtout,n,nnopt,nav,alpha,beta,
     $      ngr,hmin,hmax,nrelax,trelax,sep0,
     $      bimpact,e0,semimajoraxis,vinf2,equalmass,treloff,
     $      tresplintmuoff,
     $      nitpot,tscanon,sepfinal,nintvar,ngravprocs,
     $      qthreads,Gflag,mbh,runit,munit,
     $      cn1,cn2,cn3,cn4,cn5,cn6,cn7,computeexclusivemode,
     $      ppn,omega_spin,neos,nselfgravity,gam,reat,nusegpus,
     $      nselfgravity,ncooling,teq,tjumpahead,
     $      startfile1,startfile2,eosfile,opacityfile,throwaway
      common/integration/nintvar,neos,nusegpus,nselfgravity,ncooling
      common/inputfilenames/startfile1,startfile2,eosfile,opacityfile
      logical throwaway

      real*8 Rform, yscalconst
      common/dust/Rform
      real*8 yscalfactor(numfiltersmax)
      common /yscales/ yscalfactor

C Type of output desired:      
      WRITE (6,*) 'Select output type:  0= no output'
c      WRITE (6,*) '                     1= ppout'
c      WRITE (6,*) '                     22= expansion analysis'
c      write (6,*) '                     40= composition'
c      write (6,*) '                     41= luminosity for V838Mon'
c      write (6,*) '                     50= ascii R,P,RHO,T,MU for SM'
c      write (6,*) '                     66= optical depth'
c      write (6,*) '                     57= optical depth for polytrope'
c      write (6,*) '                     69= trajectories for <=3 stars'
c      write (6,*) '                     70= trajectories for N stars'
      write (6,*) '                     71= more robust trajectories'
c      write (6,*) '                     72= effective pot energies'
c      write (6,*) '                     73= Roche lobe radii'
c      write (6,*) '                    400= auto option 40'
      READ (5,*) IOUT

      WRITE (6,*) 'Enter Nbegin, Nfinal, dN:'
      READ (5,*) NI1,NI2,NI3
      IF(NI3.EQ.0) NI3=NI2-NI1
      IF(NI3.EQ.0) NI3=1

c      if(IOUT.eq.69) then
c         OPEN(24,file='massAndMore.out')
c         write(24,'(a12,99a13)') 'time ','m1',
c     $        'm2',
c     $        'x1','y1','z1',
c     $        'x2','y2','z2',
c     $        'separation','CE+ejecta m',
c     $        'semimajor a','eccentricity',
c     $        'unbound m ','spin1','spin2','orb. period'
c      endif

      if(IOUT.eq.71) then
c         call tabulinit
         OPEN(24,file='massAndMore.out')
         write(24,'(a12,99a13)') 'time ','m1',
     $        'm2',
     $        'x1','y1','z1',
     $        'x2','y2','z2',
     $        'separation','CE+ejecta m',
     $        'semimajor a','eccentricity',
     $        'unbound m ','spin1','spin2','orb. period'
      endif

C     Main loop over binary output files:
      DO INNIT=NI1,NI2,NI3         
         NNIT=INNIT
         IFORM = 4
         CALL READIT(NNIT,IFORM)
         if(IFORM.eq.4)then

            WRITE (6,*) 
            WRITE (6,*) 'NOUT=',NOUT,' NIT=',NIT,' T=',T
            
c            IF (IOUT.EQ.69) CALL trajectories
            IF (IOUT.EQ.71) CALL robusttrajectories
 10      endif
      ENDDO
      
      IF(IOUT.eq.69 .or. iout.eq.71) Close(24)

      END
