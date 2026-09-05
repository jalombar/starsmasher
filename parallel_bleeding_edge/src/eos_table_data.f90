module eos_table_data
  implicit none

  integer :: numrho = 0
  integer :: numu = 0
  integer :: numx = 0
  real*8 :: zzz = 0.d0
  real*8 :: rhotable1 = 0.d0
  real*8 :: utable1 = 0.d0
  real*8 :: xtable1 = 0.d0
  real*8 :: steprho = 0.d0
  real*8 :: stepu = 0.d0
  real*8 :: stepx = 0.d0
  real*8, allocatable :: eostable(:,:,:,:)
end module eos_table_data
