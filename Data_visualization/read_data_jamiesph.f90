!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR OUTPUT FROM THE GADGET CODE
!
! NOTE THAT THIS ONLY "OFFICIALLY" WORKS WITH THE PARALLEL CODE AS WE
! REQUIRE KNOWLEDGE OF THE PARTICLE SMOOTHING LENGTHS
!
! SOME CHOICES FOR THIS FORMAT CAN BE SET USING THE FOLLOWING
!  ENVIRONMENT VARIABLES:
!
! GSPLASH_USE_Z if 'YES' uses redshift in the legend instead of time
! GSPLASH_DARKMATTER_HSOFT if given a value > 0.0 will assign a
!  smoothing length to dark matter particles which can then be
!  used in the rendering
!
! the data is stored in the global array dat
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nstepsread  : number of steps read from this file
!
! dat(maxpart,maxplot,maxstep) : main data array
!
! npartoftype(maxstep): number of particles of each type in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step 
!                      (used in calc_quantities for calculating the pressure)
!
! most of these values are stored in global arrays 
! in the module 'particle_data'
!
! Partial data read implemented Nov 2006 means that columns with
! the 'required' flag set to false are not read (read is therefore much faster)
!-------------------------------------------------------------------------

subroutine read_data(rootname,istepstart,ipos,nstepsread)
  use particle_data, only:dat,iamtype,npartoftype,time,gamma,maxpart,maxcol,maxstep
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc,iformat,required,ipartialread
  use settings_page, only:legendtext
  use mem_allocation, only:alloc
  use labels, only:ih,irho
  use system_utils, only:renvironment,lenvironment
  implicit none
  integer, intent(in) :: istepstart,ipos
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  character(len=len(rootname)+10) :: datfile
  integer, dimension(maxparttypes) :: npartoftypei,Nall
  integer, dimension(:), allocatable :: iamtemp
  integer :: i,j,itype,icol,ierr
  integer :: index1,index2,indexstart,indexend,Nmassesdumped
  integer :: ncolstep,npart_max,nstep_max,ntoti
  integer :: iFlagSfr,iFlagFeedback,iFlagCool,nfiles
  logical :: iexist,reallocate
  real(doub_prec) :: timetemp,ztemp, dummy
  real(doub_prec), dimension(6) :: massoftypei
  real, dimension(:), allocatable :: dattemp1
  real :: hsoft
  integer :: ntot, nnopt, nout, nit, nav, ngr, nrelax
  real(doub_prec) :: hmin, hmax, sep0, tf, dtout, alpha, beta, eta2, trelax, dt, omega2
  real(doub_prec) :: dx, dy, dz, dm, dh, drho, dvx, dvy, dvz, dudot
  real(doub_prec) :: duth, dmmu
!  real(doub_prec) :: timeinit
!  save timeinit
!  real(doub_prec) :: timeinit
!  save timeinit
!  integer timeprompt
  integer frame
  data frame/-1/
  save frame
  real(doub_prec) :: theta
  real(doub_prec) :: gram, sec, cm, kelvin, erg, boltz
  parameter(gram=1.d0,sec=1.d0,cm=1.d0,kelvin=1.d0)
  parameter(erg=gram*cm**2/sec**2)
  parameter(boltz=1.380658d-16*erg/kelvin)
  real(doub_prec) :: munit, runit, gravconst
  parameter(munit=1.9891d33,runit=6.9599d10)
  parameter(gravconst = 6.67390d-08)

  nstepsread = 0
  npart_max = maxpart
!   if (maxparttypes.lt.6) then
!      print*,' *** ERROR: not enough particle types for GADGET data read ***'
!      print*,' *** you need to edit splash parameters and recompile ***'
!      stop
!   endif
  
  if (len_trim(rootname).gt.0) then
     datfile = trim(rootname)
  else
     print*,' **** no data read **** ' 
     return
  endif
!
!--check if first data file exists
!
  inquire(file=datfile,exist=iexist)
  if (.not.iexist) then
     print "(a)",' *** error: '//trim(datfile)//': file not found ***'    
     return
  endif
!
!--set parameters which do not vary between timesteps
!
  ndim = 3
  ndimV = 3
!
!--read data from snapshots
!  
  i = istepstart

  write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  !
  !--open data file and read data
  !
  open(11,iostat=ierr,file=datfile,status='old',form='unformatted')
  if (ierr /= 0) then
     print "(a)", '*** ERROR OPENING FILE ***'
     return
  endif
  !
  !--read header for this timestep
  !
!  read(11,iostat=ierr) npartoftypei(1:6),massoftypei,timetemp,ztemp, &
!      iFlagSfr,iFlagFeedback,Nall(1:6),iFlagCool,nfiles

  read(11, iostat = ierr) &
       ntot, nnopt, hmin, hmax, sep0, tf, dtout, nout, nit, timetemp, &
       nav, alpha, beta, eta2, ngr, nrelax, trelax, dt, omega2 

  if(omega2.ne.0.d0 .and. frame.eq.-1) then
     ! frame can equal -1 only the first time through, so this question
     ! will get asked (at most) only once
     print *, 'Period of rotating frame=',8*atan(1.d0)/omega2**0.5d0
     print *, 'Do you want movie in the inertial frame? (1=yes)'
     read(5,*) frame
  endif
  if(frame.eq.-1) frame=0

  if (ierr /= 0) then
     print "(a)", '*** ERROR READING TIMESTEP HEADER ***'
     return
  endif

!  iformat = 0
!  if (iFlagCool.gt.0) then
!     iformat = 1
!     ncolstep = 12 ! 3 x pos, 3 x vel, pmass, utherm, rho, Ne, Nh, h
!     ncolumns = ncolstep
!  else
     iformat = 0
     ncolstep = 13 ! 3 x pos, 3 x vel, pmass, rho, utherm, mean_mu, h, du/dt, temperature
     ncolumns = ncolstep  
!  endif

  irho = 8
  ih   = ncolstep

  ntoti = ntot    !int(sum(npartoftypei(1:6)))
  print*,'time             : ',timetemp
  print*,'N_total          : ',ntoti
  print*,'N data columns   : ',ncolstep

!  if (nfiles.gt.1) then
!     print*,' nfiles = ',nfiles
!     print*,'*** ERROR: read from > 1 files not implemented'
!     return
!  endif
  !
  !--if successfully read header, increment the nstepsread counter
  !
  nstepsread = nstepsread + 1
  !
  !--now read data
  !
  reallocate = .false.
  npart_max = maxpart
  nstep_max = max(maxstep,1)

  if (ntoti.gt.maxpart) then
     reallocate = .true.
     if (maxpart.gt.0) then
        ! if we are reallocating, try not to do it again
        npart_max = int(1.1*ntoti)
     else
        ! if first time, save on memory
        npart_max = int(ntoti)

!        if(int(timetemp).gt.0) then
!           print *,'Should the first data file reset to t=0? (1=yes)'
!           read(5,*) timeprompt
!           if(timeprompt.eq.1) then
!              timeinit=timetemp
!           else
!              timeinit=0.d0
!           endif
!        else
!           timeinit=0.d0
!        endif

     endif
  endif
  if (i.ge.maxstep .and. i.ne.1) then
     nstep_max = i + max(10,INT(0.1*nstep_max))
     reallocate = .true.
  endif

  !
  !--reallocate memory for main data array
  !
  if (reallocate .or. .not.(allocated(dat))) then
     call alloc(npart_max,nstep_max,max(ncolstep+ncalc,maxcol),mixedtypes=.true.)
  endif
  !
  !--copy header into header arrays
  !
  npartoftype(:,i) = 0
!  npartoftype(1,i) = ntoti
  !
  !--set time to be used in the legend
  !

  !--use this line for code time
!  time(i) = real(timetemp-timeinit) 
  time(i) = real(timetemp) 

  !
  !--read particle data
  !
  if (ntoti.gt.0) then
     !
     !--read particles' data
     !
     do j=1,ntot

        read(11, iostat = ierr) & 
             dx, dy, dz,                           &   ! positions
             dm,                                   &   ! mass
             dh, drho,                             &   ! h & rho
             dvx, dvy, dvz,                        &   ! velocities
             dummy, dummy, dummy,                  &   ! vxdot, vydot, vzdot
             duth,                                 &   ! uthermal
             dudot,                                &   ! udot
             dummy, & ! dummy, dummy, dummy,           &   ! gx, gy, gz, gpot
             dmmu,                                 &   ! mean_mu
             dummy, dummy!, dummy, dummy                ! aa, bb, cc, divv
        
        if(frame.eq.1) then
           theta=omega2**0.5d0*timetemp
! rotation is counterclockwise:
           dat(j,1,i)=dx*cos(theta)-dy*sin(theta)
           dat(j,2,i)=dx*sin(theta)+dy*cos(theta)
           dat(j,4,i)=dvx*cos(theta)-dvy*sin(theta)
           dat(j,5,i)=dvx*sin(theta)+dvy*cos(theta)
        else
           dat(j,1,i)  = dx
           dat(j,2,i)  = dy
           dat(j,4,i)  = dvx
           dat(j,5,i)  = dvy
        endif

        dat(j,3,i)  = dz
        dat(j,6,i)  = dvz

        dat(j,7,i)  = dm

        dat(j,8,i)  = drho 
        dat(j,9,i)  = duth
        dat(j,10,i) = dmmu
        dat(j,11,i) = dh
        dat(j,12,i) = dudot
        dat(j,13,i) = duth*gravconst*munit/runit/(1.5d0*boltz/dmmu)

        if(duth.ne.0.d0) then
           iamtype(j,i)=1
           npartoftype(1,i)=npartoftype(1,i)+1
        else
           iamtype(j,i)=2
           npartoftype(2,i)=npartoftype(2,i)+1
        endif
        
     end do
!    pos = 1,2,3
!    vel = 4,5,6
!    mass = 7
!    rho, u, mean_mu, h = 8, 9, 10, 11

  else
     ntoti = 1
     npartoftype(1,i) = 1
     dat(:,:,i) = 0.
  endif

!
!--now memory has been allocated, set arrays which are constant for all time
!
  gamma = 5./3.
!
!--set flag to indicate that only part of this file has been read 
!
  if (.not.all(required(1:ncolstep))) ipartialread = .true.
!
!--close data file and return
!                    
  close(unit=11)

  if (nstepsread.gt.0) then
     print*,'>> last step ntot =',sum(npartoftype(:,istepstart+nstepsread-1))
  endif
  return

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels, only:label,iamvec,labelvec,labeltype,ix,ivx,ipmass,ih,irho,ipr,iutherm
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ntypes,UseTypeInRenderings
  use geometry, only:labelcoord
  use system_utils, only:renvironment
  implicit none
  integer :: i
  real :: hsoft

  if (ndim.le.0 .or. ndim.gt.3) then
     print*,'*** ERROR: ndim = ',ndim,' in set_labels ***'
     return
  endif
  if (ndimV.le.0 .or. ndimV.gt.3) then
     print*,'*** ERROR: ndimV = ',ndimV,' in set_labels ***'
     return
  endif

  do i=1,ndim
     ix(i) = i
  enddo
  ivx = 4
  ipmass = 7
  irho = 8        ! location of rho in data array
  ipr = 0
  iutherm = 9     !  thermal energy
  ih  = 11
  
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  label(irho) = 'density'
  label(iutherm) = 'specific internal energy u'
  label(ih) = 'h'
  label(10) = 'mean_mu'
  label(ipmass) = 'particle mass'
  label(12) = 'du/dt'
  label(13) = 'temperature'
  !
  !--set labels for vector quantities
  !
  iamvec(ivx:ivx+ndimV-1) = ivx
  labelvec(ivx:ivx+ndimV-1) = 'v'
  do i=1,ndimV
     label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
  enddo
  !
  !--set labels for each particle type
  !
  ntypes = 2
  labeltype(1) = 'gas'
  labeltype(2) = 'compact object'
  UseTypeInRenderings(1) = .true.
  UseTypeInRenderings(2) = .false.



!   !
!   !--set labels of the quantities read in
!   !
!   label(ix(1:ndim)) = labelcoord(1:ndim,1)
!   label(irho) = 'density'
!   label(iutherm) = 'u'
!   label(ipmass) = 'particle mass'
  
!   if (ncolumns.gt.10) then
!      label(10) = 'Ne'
!      label(11) = 'Nh'
!      ih = 12        !  smoothing length
!   else
!      ih = 10
!   endif
!   label(ih) = 'h'
!   !
!   !--set labels for vector quantities
!   !
!   iamvec(ivx:ivx+ndimV-1) = ivx
!   labelvec(ivx:ivx+ndimV-1) = 'v'
!   do i=1,ndimV
!      label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
!   enddo
  
!   !--set labels for each particle type
!   !
!   ntypes = 5
!   labeltype(1) = 'gas'
!   labeltype(2) = 'dark matter'
!   labeltype(5) = 'star'
!   UseTypeInRenderings(1) = .true.
!   !
!   !--dark matter particles are of non-SPH type (ie. cannot be used in renderings)
!   !  unless they have had a smoothing length defined
!   !
!   hsoft = renvironment('GSPLASH_DARKMATTER_HSOFT')
!   if (hsoft.gt.tiny(hsoft)) then
!      UseTypeInRenderings(2) = .true.
!   else
!      UseTypeInRenderings(2) = .false.  
!   endif
!   UseTypeInRenderings(3:5) = .false.

!-----------------------------------------------------------
  return
end subroutine set_labels

