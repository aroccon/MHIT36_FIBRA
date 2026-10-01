module mod_forcing
use mod_common_fibm
use mod_common_mpi
!use mod_comp_tension
implicit none
private
public complagrforces,updtintermediatevel1,updtintermediatevel2, updtintermediatevel
contains
!
subroutine complagrforces
implicit none
integer :: l,p
real :: dti
!
dti = 1./dt
!
!$omp parallel default(shared) &
!$omp& private(p,l)
!$omp do
betafilament=-1./dt
do p=1,pmax

  if (ap(p)%mslv>0) then
    do l=1,NL
      !ap(p)%dxdt(l)=(ap(p)%xfpo(l) - ap(p)%xfpold(l)) / dt
      !ap(p)%dydt(l)=(ap(p)%yfpo(l) - ap(p)%yfpold(l)) / dt
      !ap(p)%dzdt(l)=(ap(p)%zfpo(l) - ap(p)%zfpold(l)) / dt


      if (istep==1.and.inertia==0) ap(p)%dydt(l)=ap(p)%zfp(l)-(lz/2.)

        
      ap(p)%integralx = ap(p)%integralx+(alphafilament*dt*(ap(p)%ul(l)-ap(p)%dxdt(l)))
      ap(p)%integraly = ap(p)%integraly+(alphafilament*dt*(ap(p)%vl(l)-ap(p)%dydt(l)))
      ap(p)%integralz = ap(p)%integralz+(alphafilament*dt*(ap(p)%wl(l)-ap(p)%dzdt(l)))

      ap(p)%fxl(l) = ap(p)%integralx+(betafilament*(ap(p)%ul(l)-ap(p)%dxdt(l)))
      ap(p)%fyl(l) = ap(p)%integraly+(betafilament*(ap(p)%vl(l)-ap(p)%dydt(l)))
      ap(p)%fzl(l) = ap(p)%integralz+(betafilament*(ap(p)%wl(l)-ap(p)%dzdt(l)))


      ! ap(p)%fzl(l)=ap(p)%fzl(l)-.00001*(ap(p)%yfp(l)-ap(p)%y)

    enddo
    !    ap(p)%fxl(:)=0.
    !    ap(p)%fyl(:)=0.
    !    ap(p)%fzl(:)=0.
    !  endif!if inertia

  endif
enddo

!$omp end parallel
!
return
end subroutine complagrforces



!
!=================================================================!
!
subroutine updtintermediatevel1
implicit none
integer i,j,k,p
real sumfx,sumfy,sumfz
real sumfx_all,sumfy_all,sumfz_all,bulk_v_periodic
real forcextot_all,forceytot_all,forceztot_all


!$omp parallel default(shared) &
!$omp& private(i,j,k) 
!$omp do

do k=0,k1
  do j=0,j1
    do i=0,i1

      dudt(i,j,k) = dudt(i,j,k)  + forcex(i,j,k)*dt!*filsub
      unew(i,j,k) = unew(i,j,k)  + forcex(i,j,k)*dt!*filsub
      !forceytot subtracted to get net zero forcing
      !-dp/dy needed to balance total drag force = -forceytot
      dvdt(i,j,k) =  dvdt(i,j,k) + (forcey(i,j,k))*dt!*filsub
      vnew(i,j,k) =  vnew(i,j,k) + (forcey(i,j,k))*dt!*filsub
      dwdt(i,j,k) = dwdt(i,j,k)  + forcez(i,j,k)*dt!*filsub
      wnew(i,j,k) = wnew(i,j,k)  + forcez(i,j,k)*dt!*filsub
      !dwdt(i,j,k) = dwdtold(i,j,k) 

    enddo
  enddo
enddo
      
!$omp end parallel
!
end subroutine updtintermediatevel1




!=================================================================!
subroutine updtintermediatevel2
implicit none
integer i,j,k,p
real sumfx,sumfy,sumfz
real sumfx_all,sumfy_all,sumfz_all,bulk_v_periodic
real forcextot_all,forceytot_all,forceztot_all



!forceytot = sumfy + wallshearnew
!
! flux imposed
!
if (solver=='f') then
  forceytot = 0.
elseif (solver=='r'.and.rkiter==1) then
  forceytot = 0.
endif

if(iniu.eq.'poi'.or.iniu.eq.'log') then


  if (flforce==1) then

    if (solver=='f') then
      forceytot = (v_bulk - bulk_v_sup)/dt
    elseif (solver=='r') then
      forceytot = forceytot + ((v_bulk - bulk_v_sup)/dt)
    endif


    elseif (flforce==2) then


    if (solver=='f') then
      forceytot = (forcextot-(visc0*(-dudyt+dudyb)/lz))
    elseif (solver=='r') then
      forceytot = forceytot+((forcextot-(visc0*(-dudyt+dudyb)/lz))*rkcoeffab(rkiter))
    endif

  endif



  if (myid==0.and.mod(istep,10)==0.and.solver=='f') then
    open(11,file=datadir//'pressure_gradient.txt',position='append')
    write(11,'(2E16.8)') time,forceytot
    close(11)
  endif

endif


!$omp parallel default(shared) &
!$omp& private(i,j,k) 
!$omp do

do k=0,k1
  do j=0,j1
    do i=0,i1

      !      dudt(i,j,k) = dudt(i,j,k)  + forcex(i,j,k)*dt!*filsub
      !      unew(i,j,k) = unew(i,j,k)  + forcex(i,j,k)*dt!*filsub
      !      forceytot subtracted to get net zero forcing
      !     -dp/dy needed to balance total drag force = -forceytot
      dvdt(i,j,k) =  dvdt(i,j,k) -forceytot*dt!*filsub
      vnew(i,j,k) =  vnew(i,j,k) -forceytot*dt!*filsub
      !      dwdt(i,j,k) = dwdt(i,j,k)  + forcez(i,j,k)*dt!*filsub
      !      wnew(i,j,k) = wnew(i,j,k)  + forcez(i,j,k)*dt!*filsub
      !     dwdt(i,j,k) = dwdtold(i,j,k) 

    enddo
  enddo
enddo

      
!$omp end parallel
!
end subroutine updtintermediatevel2


subroutine updtintermediatevel
implicit none
integer i,j,k,p
real sumfx,sumfy,sumfz
real sumfx_all,sumfy_all,sumfz_all,bulk_v_periodic
real forcextot_all,forceytot_all,forceztot_all



!forceytot = sumfy + wallshearnew
!
! flux imposed
!
forcextot = 0.
forceytot = 0.
forceytot = 0.


if(iniu.eq.'poi'.or.iniu.eq.'log') then

  bulk_v_periodic=bulk_v_sup*cos(2.*pi*(1./60.)*time)
  forceytot = -1.0*(bulk_v_sup-v_bulk)/dt
endif

!print*,"--------------forcetot"
!print*, real(forceytot,4)

!$omp parallel default(shared) &
!$omp& private(i,j,k) 
!$omp do

do k=0,k1
  do j=0,j1
    do i=0,i1

      dudt(i,j,k) = dudtold(i,j,k)
      !forceytot subtracted to get net zero forcing
      !-dp/dy needed to balance total drag force = -forceytot
      dvdt(i,j,k) =  dvdtold(i,j,k) - forceytot*dt
      dwdt(i,j,k) = dwdtold(i,j,k)

    enddo
  enddo
enddo
!$omp end parallel
!
end subroutine updtintermediatevel

!
end module mod_forcing
