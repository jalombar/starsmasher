!***********************************************************************
!     Whether the groups of hgroup.f90 can run on this processor.
!
!     hgroup.f90 is compiled for AVX2 on x86-64 (HGFLAGS in the Makefile,
!     which then also defines HG_AVX2 for this file).  A processor without
!     it, or a virtual machine that hides it, stops on the first AVX2
!     instruction, so this file, which is compiled for any x86-64, reads the
!     flags that /proc/cpuinfo lists the first time it is called and answers
!     .false. when avx2 is not among them; rank 0 then says that the groups
!     are not used.  Where there is no /proc/cpuinfo the flags the file was
!     compiled with are trusted.
!***********************************************************************
      logical function hg_cpu_ok(rank)
      implicit none
      integer, intent(in) :: rank
      logical, save :: checked = .false., ok = .true.
#ifdef HG_AVX2
      integer :: u, ios
      character(len=16384) :: line
      if (.not. checked) then
         checked = .true.
         open(newunit=u, file='/proc/cpuinfo', status='old', action='read', iostat=ios)
         if (ios == 0) then
            do
               read(u, '(a)', iostat=ios) line
               if (ios /= 0) exit
               if (line(1:5) == 'flags') then
                  ok = index(line, ' avx2') > 0
                  exit
               endif
            enddo
            close(u)
         endif
         if (.not. ok .and. rank == 0) then
            write(*,'(a)') 'hgroup: this processor has no AVX2, which hgroup.f90 is compiled for: '// &
                 'the groups are not used (hgroup=0)'
         endif
      endif
#endif
      hg_cpu_ok = ok
      end function hg_cpu_ok
