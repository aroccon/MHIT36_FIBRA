module mod_concentration
use mod_param_fibm
use mod_common_fibm
use mod_common_mpi
use mod_kernel
use mod_bound
use mod_output
use mod_output
use decomp_2d
implicit none
private
public concentration
contains


!v has only this subroutine. Purpose??
!v 
subroutine concentration
implicit none
integer :: i,j,k,l,p
integer :: ilow,ihigh,jlow,jhigh,klow,khigh
real,dimension(itot,ktot)::vf2,vf2t
real :: coorx,coory,coorz,fctot
real :: coorxs,coorys,coorzs
real :: coorxfp,cooryfp,coorzfp
real :: kernelx,kernely,kernelz
real :: kernelxs,kernelys,kernelzs
real :: forcex_sc,forcey_sc,forcez_sc
type pneighbor
  real :: x,y,z
  real, dimension(1:nno) :: xfp,yfp,zfp, qn1, qn2, qn3, qn4
!  real, dimension(1:nl) :: sumu
end type pneighbor
type(pneighbor), dimension(0:8,1:npmax) :: anb
integer :: nb,nbsend,nbrecv
integer :: nrrequests
integer :: arrayrequests(1:3) !3=3*1 (master might have 3 slaves)
integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
real :: boundleftnb,boundrightnb,boundfrontnb,boundbacknb
real :: dVlagrdVeul
character(7) filenumber
!real, dimension(1:nl,1:npmax) :: sumu
integer :: idp,tag
!

dVeul = dx*dy*dz

dVlagrdVeul = dVlagr/dVeul
!
forcexold=forcex
forceyold=forcey
forcezold=forcez
!$omp workshare
anb(0,1:pmax)%x = ap(1:pmax)%x
anb(0,1:pmax)%y = ap(1:pmax)%y
anb(0,1:pmax)%z = ap(1:pmax)%z
anb(1:8,1:pmax)%x = 0.
anb(1:8,1:pmax)%y = 0.
anb(1:8,1:pmax)%z = 0.
forall(l=1:nno,p=1:pmax)
  anb(0,p)%xfp(l) = ap(p)%xfp(l)
  anb(0,p)%yfp(l) = ap(p)%yfp(l)
  anb(0,p)%zfp(l) = ap(p)%zfp(l)
  anb(0,p)%qn1(l) = ap(p)%qn1(l)
  anb(0,p)%qn2(l) = ap(p)%qn2(l)
  anb(0,p)%qn3(l) = ap(p)%qn3(l)
  anb(0,p)%qn4(l) = ap(p)%qn4(l)
end forall

forall(l=1:nl,i=1:8,p=1:pmax)
  anb(i,p)%xfp(l) = 0.
  anb(i,p)%yfp(l) = 0.
  anb(i,p)%zfp(l) = 0.
  anb(0,p)%qn1(l) = 1.
  anb(0,p)%qn2(l) = 0.
  anb(0,p)%qn3(l) = 0.
  anb(0,p)%qn4(l) = 0.
end forall


vf=0.


!$omp end workshare
!
do p=1,pmax
  nrrequests = 0
  do nb=1,8
    nbsend = nb     ! rank of process which receives data ('sends data to neighbor nbsend')
    idp = abs(ap(p)%mslv)
    tag = idp*10+nbsend
    nbrecv  = nb+4  ! rank of process which sends data ('receives data from neighbor nbrecv')
    if (nbrecv .gt. 8) nbrecv = nbrecv - 8
    if (ap(p)%mslv .gt. 0) then
      ! myid is master of particle ap(p)%mslv
      if (ap(p)%nb(nbsend) .eq. 1) then
        ! neighbor(nbsend) is rank of slave for particle ap(p)%mslv
        if ( neighbor(nbsend) .eq. myid ) then
          ! process might be both master and slave of same particle due to periodic b.c.'s
          anb(nbrecv,p)%x = anb(0,p)%x
          anb(nbrecv,p)%y = anb(0,p)%y
          anb(nbrecv,p)%z = anb(0,p)%z
          !$omp parallel default(shared) &
          !$omp&private(l)
          !$omp do

          do l=1,nl
            anb(nbrecv,p)%xfp(l) = anb(0,p)%xfp(l)
            anb(nbrecv,p)%yfp(l) = anb(0,p)%yfp(l)
            anb(nbrecv,p)%zfp(l) = anb(0,p)%zfp(l)
            anb(nbrecv,p)%qn1(l) = anb(0,p)%qn1(l)
            anb(nbrecv,p)%qn2(l) = anb(0,p)%qn2(l)
            anb(nbrecv,p)%qn3(l) = anb(0,p)%qn3(l)
            anb(nbrecv,p)%qn4(l) = anb(0,p)%qn4(l)

          enddo
          if (nbrecv .eq. 1) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
            !$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
            enddo
          endif
          if (nbrecv .eq. 2) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
            !$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
            enddo
          endif
          if (nbrecv .eq. 3) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
            !$omp do
            do l=1,nl
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
            enddo
          endif
          if (nbrecv .eq. 4) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
            !$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
            enddo
          endif
          if (nbrecv .eq. 5) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
            !$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
            enddo
          endif
          if (nbrecv .eq. 6) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
            !$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
            enddo
          endif
          if (nbrecv .eq. 7) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
            !$omp do
            do l=1,nl
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
            enddo
          endif
          if (nbrecv .eq. 8) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
            !$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
            enddo
          endif
          !$omp end parallel
        else

          nrrequests = nrrequests + 1
          call MPI_ISEND(anb(0,p)%x,3*(1+1*nl),MPI_REAL8,neighbor(nbsend), &
                         tag,comm_cart,arrayrequests((nrrequests-1)+1),ierr)
          ! send x,y,z,xfp,yfp,zfp,fxl,fyl,fzl -> 3*(1+2*nl) contiguous info
          ! (see definition of type pneighbor in the begining of the subroutine)
        endif
      endif
    endif
    if (ap(p)%mslv .lt. 0) then
      ! myid is slave of particle -ap(p)%mslv
      if (ap(p)%nb(nbrecv) .eq. 1) then
        ! ap(p)%nb(nbrecv) is rank of master of particle -ap(p)%mslv
        nrrequests = nrrequests + 1
        call MPI_IRECV(anb(nbrecv,p)%x,3*(1+1*nl),MPI_REAL8,neighbor(nbrecv), &
                       tag,comm_cart,arrayrequests((nrrequests-1)+1),ierr)
        ! recv x,y,z,xfp,yfp,zfp,fxl,fyl,fzl -> 3*(1+2*nl) contiguous info
        ! (see definition of type pneighbor in the begining of the subroutine)
      endif
    endif
  enddo ! do nb=
  call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,ierr)
enddo

!
! second step: recompute particle positions for slaves due to periodic b.c.'s.
! required: (part of) particle within domain bounds of slave process.
!
!$omp parallel default(shared) &
!$omp&private(p,l,nbrecv,boundleftnb,boundbacknb,boundrightnb,boundfrontnb)  
!!$omp do schedule(dynamic)
!$omp do
do p=1,pmax
  if (ap(p)%mslv .lt. 0) then
    ! myid is slave of particle -ap(p)%mslv
    nbrecv=1
    if (ap(p)%nb(nbrecv) .eq. 1) then
      ! ap(p)%nb(nbrecv) is rank of master of particle -ap(p)%mslv
      boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left boundary of neighbor nb
      if (anb(nbrecv,p)%x .lt. boundleftnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
        enddo
      endif
    endif
    nbrecv=2
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left boundary of neighbor nb
      boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back boundary of neighbor nb
      if (anb(nbrecv,p)%x .lt. boundleftnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
        enddo
      endif
      if (anb(nbrecv,p)%y .gt. boundbacknb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
        enddo
      endif
    endif
    nbrecv=3
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back boundary of neighbor nb
      if (anb(nbrecv,p)%y .gt. boundbacknb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
        enddo
      endif
    endif
    nbrecv=4
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
      boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
      if (anb(nbrecv,p)%x .gt. boundrightnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
        enddo
      endif
      if (anb(nbrecv,p)%y .gt. boundbacknb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
        enddo
      endif
    endif
    nbrecv=5
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
      if (anb(nbrecv,p)%x .gt. boundrightnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
        enddo
      endif
    endif
    nbrecv=6
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
      boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
      if (anb(nbrecv,p)%x .gt. boundrightnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
        enddo
      endif
      if (anb(nbrecv,p)%y .lt. boundfrontnb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
        enddo
      endif
    endif
    nbrecv=7
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
      if (anb(nbrecv,p)%y .lt. boundfrontnb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
        enddo
      endif
    endif
    nbrecv=8
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
      boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
      if (anb(nbrecv,p)%x .lt. boundleftnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
        enddo
      endif
      if (anb(nbrecv,p)%y .lt. boundfrontnb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
        enddo
      endif
    endif
  endif
enddo
!$omp end parallel
!
! third step: perform partial integration.
!
!$omp workshare
!forcex(:,:,:) = 0.
!forcey(:,:,:) = 0.
!forcez(:,:,:) = 0.
!!!dudt(:,:,:) = 0.
!!!dvdt(:,:,:) = 0.
!!!dwdt(:,:,:) = 0.
!$omp end workshare
!
!$omp parallel default(shared) &
!$omp&private(p,nbrecv,l,coorxfp,cooryfp,coorzfp,ilow,ihigh,jlow,jhigh,klow,khigh) &
!$omp&private(i,j,k,nbsend,nb) &
!$omp&private(coorzs,coorz,kernelzs,kernelz,coorys,coory,kernelys,kernely,coorxs,coorx,kernelxs,kernelx,forcex_sc,forcey_sc,forcez_sc)
!!$omp do schedule(dynamic)
!$omp do 
do p=1,pmax
  if (ap(p)%mslv .ne. 0) then
    ! myid is master or slave of particle abs(ap(p)%mslv)
    do nb=0,8
      if((nb.gt.0.and.ap(p)%mslv.lt.0.and.ap(p)%nb(nb).eq.1) .or. & !slave
         (nb.eq.0.and.ap(p)%mslv.gt.0) .or. & !pure master
         (nb.gt.0.and.ap(p)%mslv.gt.0.and.ap(p)%nb(nb).eq.1.and.neighbor(nb) .eq. myid)) then
        nbrecv = nb
        ! ap(p)%nb(nbrecv) is rank of master of particle -ap(p)%mslv
        if ( neighbor(nb) .eq. myid .and. nb.gt.0.and.ap(p)%mslv.gt.0) then
          nbrecv = nb + 4
          if (nbrecv .gt. 8) nbrecv = nbrecv-8
        endif
        do l=1,nl
          coorxfp = (anb(nbrecv,p)%xfp(l)-boundleftmyid)*dxi
          cooryfp = (anb(nbrecv,p)%yfp(l)-boundfrontmyid)*dyi
          coorzfp =  anb(nbrecv,p)%zfp(l)*dzi
          ilow  = nint( coorxfp - 2. )
          ihigh = nint( coorxfp + 2. )
          jlow  = nint( cooryfp - 2. )
          jhigh = nint( cooryfp + 2. )
          klow  = nint( coorzfp - 2. )
          khigh = nint( coorzfp + 2. )
          if (ilow .lt. 1) ilow = 1
          if (jlow .lt. 1) jlow = 1
          if (klow .lt. 1) klow = 1
          if (ihigh .gt. imax) ihigh = imax
          if (jhigh .gt. jmax) jhigh = jmax
          if (khigh .gt. kmax) khigh = kmax
          do k=klow,khigh
            coorzs   = (1.*k)-coorzfp
            coorz    = coorzs          !vert.    distance in grid points
            kernelzs = kernel(coorzs)
            kernelz  = kernel(coorz)
            do j=jlow,jhigh
              coorys   = (1.*j)-cooryfp
              coory    = coorys         !spanw.   distance in grid points
              kernelys = kernel(coorys)
              kernely  = kernel(coory)        
              do i=ilow,ihigh
                coorxs   = (1.*i)-coorxfp
                coorx    = coorxs       !streamw. distance in grid points
                kernelxs = kernel(coorxs)
                kernelx  = kernel(coorx)

                if (l==1.or.l==NL) then

                  !$omp atomic
                  vf(i,j,k) = vf(i,j,k) + kernelx*kernely*kernelz*((aspectratio**2.)*pi*.25)* &
                  (ds*.5/dVeul)

                else
                  !$omp atomic
                  vf(i,j,k) = vf(i,j,k) + kernelx*kernely*kernelz*((aspectratio**2.)*pi*.25)* &
                  (ds/dVeul)

                endif

              enddo
            enddo
          enddo
        enddo !do l=
      endif
    enddo !do nbrecv=
  endif
enddo
!$omp end parallel
!

! communicate data in x direction (periodic b.c.'s incorporated)
!
!v call updthalos(vf,1)
!call updthalos(forcey,1)
!call updthalos(forcez,1)
!
! communicate data in y direction (periodic b.c.'s incorporated)
!
!v call updthalos(vf,2)
!call updthalos(forcey,2)
!call updthalos(forcez,2)
!
vf2=0.
vf2t=0.
do i=1,imax
  do j=1,jmax
    do k=1,kmax
      vf2t(i+zstart(1)-1,k)=vf2t(i+zstart(1)-1,k)+vf(i,j,k)
    enddo
  enddo
enddo

do i=1,itot
  do k=0,ktot
    call mpi_allreduce(vf2t(i,k),vf2(i,k),1,mpi_real8,mpi_sum,comm_cart,ierr)
  enddo
enddo

vf2=vf2/jtot


!!if (myid==0) then
!!  write(filenumber,'(i7.7)') istep
!!  open(40,file=datadir//'vf'//filenumber//'.txt')
!!  do i=1,itot
!!    do k=1,ktot
!!      write(40,'(3E15.7)')  (i-.5)*dx,(k-.5)*dz,vf2(i,k)
!!    enddo
!!  enddo
!!  close(40)
!!endif


return
end subroutine concentration
!
end module mod_concentration
