! Define the fortran interface to starsmasher here
! 
! Contains subroutine to set all values from Python for now

! Subroutine responsible for setting all of the values that the user puts in python
subroutine PythonSetValues( Pndisplace, Pdisplacex, Pdisplacey, Pdisplacez, Psemimajoraxis&
        &,Pbimpact,Pe0,Pvinf2, Ptf, Pdtout, Pn, Pgflag, Pnnopt, Pnav, Palpha, Pbeta, Pngr,& 
        &Phco, Phfloor,Pnrelax,Ptrelax,Psep0,Pequalmass,Ptreloff,Ptresplintmuoff,Pnitpot,&
        &Ptscanon,Psepfinal,Pnintvar,Pngravprocs,Pqthreads,Pmbh,Prunit,Pmunit,Pcn1,Pcn2,Pcn3,&
        &Pcn4,Pcn5,Pcn6,Pcn7,Pcomputeexclusivemode,Pomega_spin,Pppn,Pneos,Pnselfgravity,Pgam,&
        &Preat, Pstarmass,Pstarradius, Pncooling, Pnkernel, Pteq, Pjumpahead, Pstartfile1,Pstartfile2,&
        &Peosfile,Popacityfile,Pprofilefile, Pthrowaway,Pstellarevolutioncodetype,Psimulationtype)
    
     include 'starsmasher.h'


!    common/displace/displacex,displacey,displacez,ndisplace
!    namelist/input/ tf,dtout,n,nnopt,nav,alpha,beta,ngr,hco,hfloor,& 
!    & 


end subroutine
