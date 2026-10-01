!
! SPDX-License-Identifier: MIT
!
module mod_common_mpi
  !
  use mpi
  use mod_types, only: rp
  !
  implicit none
  !
  integer :: myid
  integer :: comm_cart,ierr
  integer :: status(MPI_STATUS_SIZE)
#if defined(_OPENACC)
  real(rp), allocatable, dimension(:,:) :: xsl_buf, xrl_buf, xsr_buf, xrr_buf, &
                                           ysr_buf, yrr_buf, ysl_buf, yrl_buf, &
                                           zsr_buf, zrr_buf, zsl_buf, zrl_buf
  !@cuf integer :: mydev
  !
#if defined(_GPU_MPI)
  attributes(device)  :: xsl_buf, xrl_buf, xsr_buf, xrr_buf, &
                         ysr_buf, yrr_buf, ysl_buf, yrl_buf, &
                         zsr_buf, zrr_buf, zsl_buf, zrl_buf
#else
  attributes(managed) :: xsl_buf, xrl_buf, xsr_buf, xrr_buf, &
                         ysr_buf, yrr_buf, ysl_buf, yrl_buf, &
                         zsr_buf, zrr_buf, zsl_buf, zrl_buf
#endif
#endif
  !
  integer, dimension(3) :: ijk_start,ijk_start_x,ijk_start_y,ijk_start_z,n_x,n_y,n_z
  integer :: left,right,front,back,top,bottom,ipencil
  !
#if defined(_IBM)
	integer                 :: xhalo,yhalo,xhalocol,yhalocol,xhalocolint,yhalocolint!,restartp,rankcw
	!logical                :: reorder
	integer                 :: request             !v assgn zredistribute
	integer                 :: rightfront, leftfront, leftback, rightback !v assgn initmpi
	integer, dimension(0:8) :: neighbor            !v assgn initmpi
  integer, dimension(1:2) :: coords              !v assgn initmpi 
  real                    :: boundleftmyid,boundfrontmyid !v assgn initmpi
	!real                   :: theta_Spagetti
#endif
!
end module mod_common_mpi
