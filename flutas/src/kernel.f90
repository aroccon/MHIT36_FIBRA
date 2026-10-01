module mod_kernel
implicit none
private
public kernel,kernelchk,kerneltest
contains
!
real function kernel(r)
!
! Discretized delta function. (N.B.: multiplied by the grid spacing!)
!
implicit none
real :: r
!
if (abs(r) .gt. 1.5) then
  kernel=0.
else
  if (abs(r) .gt. 0.5) then
    kernel = (1./6.)*(5.-3.*abs(r)-sqrt(-3.*((1.-abs(r))**2)+1.))
  else
    kernel = (1./3.)*(1.+sqrt(-3.*(r**2)+1.))
  endif
endif
!
return
end function kernel
!
integer function kernelchk(r)
!
! checks if r is within the kernel stencil. if so kernelchk(r) = 1
!
implicit none
real :: r
!
if (abs(r).gt.1.5) then
  kernelchk = 0
else
  kernelChk = 1
endif
!
return
end function kernelchk
!
subroutine kerneltest(sum)
integer :: i,j,k
real :: coorx,coory,coorz
real, intent(out) :: sum
!
! tests if function kernel is working properly (sum ~ 1)
!
sum = 0.
do k=-25,25
  coorz=k/12.
  do j=-25,25
    coory=j/12.
    do i=-25,25
      coorx=i/12.
      sum = sum + (1./12.)*(1./12.)*(1./12.)*kernel(coorx)*kernel(coory)*kernel(coorz)
    enddo
  enddo
enddo
!
return
end subroutine kerneltest
!
end module mod_kernel
