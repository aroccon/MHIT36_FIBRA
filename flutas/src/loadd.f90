module mod_loadpart
use mpi
use mod_param_fibm, only: np, nl, ndims, lx, ly, lz, dims, &
											offset, datadir, nno
use mod_common_fibm,   only: ap, npmstr, pmax, pmax, l, npmax
use mod_common_mpi!,  only: MPI_REAL8, comm_cart, ierr, mpi_finalize, &
									!			MPI_INTEGER, MPI_SUM, myid
implicit none
private
public loadpart
contains
!
subroutine loadpart(in,nr)
implicit none

!integer, intent(IN) :: in,nr
integer, intent(in) :: in, nr
type particle_restart
 real :: x,y,z, integralx, integraly,integralz
 real, dimension(nno) :: xll, yll, zll, qn1, qn2, qn3, qn4, & 
                         xfpo,yfpo,zfpo, qno1, qno2, qno3, qno4, &
                         xfpold,yfpold,zfpold, qnold1, qnold2, qnold3, qnold4, & 
                         dxdtl,dydtl,dzdtl, omg1, omg2, omg3,  & 
												 dxdto,dydto,dzdto, omgo1, omgo2, omgo3,  & 
												 ua,va,wa, omgd1, omgd2, omgd3,   &        
												 uao,vao,wao, omgdo1, omgdo2, omgdo3, &
												 fxll,fyll,fzll, mxl, myl, mzl,  &         
												 ull,vll,wll, omgl1, omgl2, omgl3!,  &      
!                      dxds,dyds,dzds, &                        
!                      fcx,fcy,fcz, &              
!                      xt1,yt1,zt1, &             
!                      xt2,yt2,zt2, &
!                      xt3,yt3,zt3, &
!                      xt4,yt4,zt4, &
!                      xtc,ytc,ztc, &             
!                      fx1,fy1,fz1, &            
!                      fx2,fy2,fz2, &
!                      fx3,fy3,fz3, &
!                      fx4,fy4,fz4!, &
                        !cp                        
! amount to be communicated: 6+57*nno	
! real, dimension (nl) :: xfp, yfp, zfp, &
!												dxdt, dydt, dzdt, &
!												ul, vl, wl, &
!                       fxl,fyl,fzl
end type particle_restart
!type(particle_restart), dimension(np/5) :: glob, glob_all
type(particle_restart), dimension(np) :: glob, glob_all

integer, parameter :: skip = 6+57*nno
integer p,i,idp,pp
integer :: fh,rcc,pcc
character(1):: fl
integer :: proccoords(1:ndims),procrank
real :: leftbound,rightbound,frontbound,backbound
real :: xdum
integer :: lenr
real :: dist,angle
real :: xp,yp
real :: ax
real :: ay
integer :: count_mstr_all, icp
character(len=7)  :: istepchar
character(len=4)  :: rankpr
character(len=11) :: tempstr2
integer :: counter,count_slve_loc
integer(kind=MPI_OFFSET_KIND) :: filesize,disp
!
inquire (iolength=lenr) xdum
write(istepchar,'(i7.7)') nr
!
!
!
if (in.eq.0) then
!
ap(1:npmax)%mslv = 0
ap(1:npmax)%x = 0.
ap(1:npmax)%y = 0.
ap(1:npmax)%z = 0.
ap(1:npmax)%integralx = 0.
ap(1:npmax)%integraly = 0.
ap(1:npmax)%integralz = 0.

forall (p=1:npmax)
	ap(p)%xfp(:)  = 0.
  ap(p)%yfp(:)  = 0.
  ap(p)%zfp(:)  = 0.

  ap(p)%fxl(:)  = 0.
  ap(p)%fyl(:)  = 0.
  ap(p)%fzl(:)  = 0.
  
	ap(p)%dxdt(:)  = 0.
  ap(p)%dydt(:)  = 0.
  ap(p)%dzdt(:)  = 0.
  
	ap(p)%ul(:)    = 0.
  ap(p)%vl(:)    = 0.
  ap(p)%wl(:)    = 0.

	ap(p)%xll(:)  = 0.
  ap(p)%yll(:)  = 0.
  ap(p)%zll(:)  = 0.
  ap(p)%qn1(:)  = 0.
  ap(p)%qn2(:)  = 0.
  ap(p)%qn3(:)  = 0.
  ap(p)%qn4(:)  = 0.

  ap(p)%xfpo(:) = 0.
  ap(p)%yfpo(:) = 0.
  ap(p)%zfpo(:) = 0.
  ap(p)%qno1(:) = 0.
  ap(p)%qno2(:) = 0.
  ap(p)%qno3(:) = 0.
  ap(p)%qno4(:) = 0.

  ap(p)%xfpold(:) = 0.
  ap(p)%yfpold(:) = 0.
  ap(p)%zfpold(:) = 0.
  ap(p)%qnold1(:) = 0.
  ap(p)%qnold2(:) = 0.
  ap(p)%qnold3(:) = 0.
  ap(p)%qnold4(:) = 0.

  ap(p)%dxdtl(:)  = 0.
  ap(p)%dydtl(:)  = 0.
  ap(p)%dzdtl(:)  = 0.
  ap(p)%omgd1(:) = 0.
  ap(p)%omgd2(:) = 0.
  ap(p)%omgd3(:) = 0.

  ap(p)%ua(:)    = 0.
  ap(p)%va(:)    = 0.
  ap(p)%wa(:)    = 0.
  ap(p)%omgd1(:) = 0.
  ap(p)%omgd2(:) = 0.
  ap(p)%omgd3(:) = 0.

  ap(p)%ull(:)    = 0.
  ap(p)%vll(:)    = 0.
  ap(p)%wll(:)    = 0.
  ap(p)%omgl1(:) = 0.
  ap(p)%omgl2(:) = 0.
  ap(p)%omgl3(:) = 0.

  ap(p)%fxll(:)  = 0.
  ap(p)%fyll(:)  = 0.
  ap(p)%fzll(:)  = 0.
  ap(p)%mxl(:)  = 0.
  ap(p)%myl(:)  = 0.
  ap(p)%mzl(:)  = 0.

  ap(p)%dxdto(:)  = 0.  !v I added
  ap(p)%dydto(:)  = 0.
  ap(p)%dzdto(:)  = 0.
  ap(p)%omgdo1(:) = 0.
  ap(p)%omgdo2(:) = 0.
  ap(p)%omgdo3(:) = 0.

  ap(p)%uao(:)    = 0.  !v I added
  ap(p)%vao(:)    = 0.
  ap(p)%wao(:)    = 0.
  ap(p)%omgdo1(:) = 0.
  ap(p)%omgdo2(:) = 0.
  ap(p)%omgdo3(:) = 0.

  ap(p)%nb(1:8) = 0
end forall




i = 0
npmstr = 0

do rcc=1,1
write(fl,'(i1.1)') rcc

		
  if(myid.eq.0) then
    !v open(20,file=datadir//'allpart'//fl//'data'//istepchar,access='direct',recl=np/5*skip*lenr)
    !open(20,file=datadir//'allpart'//fl//'data'//istepchar,access='direct',recl=np*skip*lenr)
		!print*, "----np, skip, lenr, np*skip*lenr, nr, istepchar:   ", np, skip, lenr, np*skip*lenr, nr, istepchar

		open(20,file=datadir//'allpart'//fl//'data'//istepchar,access='direct',recl=np*skip*lenr)
    read(20,rec=rcc) glob
		!close(20)
  endif
  !
  ! rank 0 broadcasts the particle data to all the others
  !
  !call MPI_BCAST(glob(1)%x,np/5*skip,MPI_REAL8,0,comm_cart,ierr)
  call MPI_BCAST(glob(1)%x,np*skip,MPI_REAL8,0,comm_cart,ierr)
  ! 
  ! Determine master and slave processes for each particle.
  !
  ! initialisation
!  !
!ap(p(1:npmax)%x = 0.
!ap(1:npmax)%y = 0.
!ap(1:npmax)%z = 0.
!ap(1:npmax)%integralx = 0.
!ap(1:npmax)%integraly = 0.
!ap(1:npmax)%integralz = 0.
!
!forall (p=1:npmax)
!  ap(p)%xfp(:)  = 0.
!  ap(p)%yfp(:)  = 0.
!  ap(p)%zfp(:)  = 0.
!  ap(p)%ul(:)   = 0.
!  ap(p)%vl(:)   = 0.
!  ap(p)%wl(:)   = 0.
!  ap(p)%fxl(:)  = 0.
!  ap(p)%fyl(:)  = 0.
!  ap(p)%fzl(:)  = 0.
!  ap(p)%dxdt(:) = 0.
!  ap(p)%dydt(:) = 0.
!  ap(p)%dzdt(:) = 0.
!  ap(p)%xfpo(:)  = 0.
!  ap(p)%yfpo(:)  = 0.
!  ap(p)%zfpo(:)  = 0.
!  ap(p)%xfpold(:)  = 0.
!  ap(p)%yfpold(:)  = 0.
!  ap(p)%zfpold(:)  = 0.
!  ap(p)%ua(:)   = 0.
!  ap(p)%va(:)   = 0.
!  ap(p)%wa(:)   = 0.
!  ap(p)%nb(1:8)   = 0
!end forall1:npmax)%mslv = 0
  


  ax = 0.5
  ay = 0.5

  !do p=1,np/5
  do p=1,np

    if (glob(p)%x.lt.0..or.glob(p)%x.gt.lx .or. &
        glob(p)%y.lt.0..or.glob(p)%y.gt.ly .or. &
        glob(p)%z.lt.0..or.glob(p)%z.gt.lz) then
      if (myid.eq.0) then
        write(6,*) 'Fatal error in initialisation of particle positions - '
        write(6,*) 'particle outside the domain!'
        write(6,*) 'Program aborted...'
      endif
      call mpi_finalize(ierr)
      stop
    endif
    if (glob(p)%x.eq.lx) ax = 0.51
    if (glob(p)%x.eq.0) ax = 0.49
    if (glob(p)%y.eq.ly) ay = 0.51
    if (glob(p)%y.eq.0) ay = 0.49

    proccoords(1) = nint(glob(p)%x*dims(1)/lx - ax)
    proccoords(2) = nint(glob(p)%y*dims(2)/ly - ay)
    leftbound     = (proccoords(1)  )*lx/(1.*dims(1)) ! left  boundary of particle's master
    rightbound    = (proccoords(1)+1)*lx/(1.*dims(1)) ! right boundary of particle's master
    frontbound    = (proccoords(2)  )*ly/(1.*dims(2)) ! front boundary of particle's master
    backbound     = (proccoords(2)+1)*ly/(1.*dims(2)) ! back  boundary of particle's master
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid.eq.procrank) then
      npmstr = npmstr + 1
      i = i + 1

      !ap(i)%mslv = ((rcc-1)*(np/5))+p
      ap(i)%mslv = ((rcc-1)*(np))+p

      ap(i)%x = glob(p)%x
      ap(i)%y = glob(p)%y
      ap(i)%z = glob(p)%z

      ap(i)%integralx = glob(p)%integralx
      ap(i)%integraly = glob(p)%integraly
      ap(i)%integralz = glob(p)%integralz

			ap(i)%xll(1:nno)   = glob(p)%xll(1:nno)
      ap(i)%yll(1:nno)   = glob(p)%yll(1:nno)
      ap(i)%zll(1:nno)   = glob(p)%zll(1:nno)
      ap(i)%qn1(1:nno)   = glob(p)%qn1(1:nno)
      ap(i)%qn2(1:nno)   = glob(p)%qn2(1:nno)
      ap(i)%qn3(1:nno)   = glob(p)%qn3(1:nno)
      ap(i)%qn4(1:nno)   = glob(p)%qn4(1:nno)
      
      ap(i)%yfpo(1:nno)  = glob(p)%yfpo(1:nno)
      ap(i)%zfpo(1:nno)  = glob(p)%zfpo(1:nno)
      ap(i)%qno1(1:nno)  = glob(p)%qno1(1:nno)
      ap(i)%qno2(1:nno)  = glob(p)%qno2(1:nno)
      ap(i)%qno3(1:nno)  = glob(p)%qno3(1:nno)
      ap(i)%qno4(1:nno)  = glob(p)%qno4(1:nno)

      ap(i)%xfpold(1:nno)= glob(p)%xfpold(1:nno)
      ap(i)%yfpold(1:nno)= glob(p)%yfpold(1:nno)
      ap(i)%zfpold(1:nno)= glob(p)%zfpold(1:nno)
      ap(i)%qnold1(1:nno)= glob(p)%qnold1(1:nno)
      ap(i)%qnold2(1:nno)= glob(p)%qnold2(1:nno)
      ap(i)%qnold3(1:nno)= glob(p)%qnold3(1:nno)
      ap(i)%qnold4(1:nno)= glob(p)%qnold4(1:nno)
			
			!print*, "---value of glob%qn1 from saved data"
      !do icp =1,nno
			!	!print*, "---glob(p)%qN1(1:nno)", real(glob(p)%qn1(icp),4)
			!	print*, "---ap(i)%qN1(1:nno)", real(ap(i)%qn1(icp),8)
			!enddo
      !do icp =1,nno
			!	!print*, "---glob(p)%qN1(1:nno)", real(glob(p)%qn1(icp),4)
			!	print*, "---ap(i)%qNo1(1:nno)", real(ap(i)%qno1(icp),8)
			!enddo
      !do icp =1,nno
			!	!print*, "---glob(p)%qN1(1:nno)", real(glob(p)%qn1(icp),4)
			!	print*, "---ap(i)%qNold1(1:nno)", real(ap(i)%qnold1(icp),8)
			!enddo
			!print*, "---value of glob%qn1 from saved data"

      ap(i)%xfpo(1:nno)  = glob(p)%xfpo(1:nno)

      ap(i)%dxdtl(1:nno) = glob(p)%dxdtl(1:nno)
      ap(i)%dydtl(1:nno) = glob(p)%dydtl(1:nno)
      ap(i)%dzdtl(1:nno) = glob(p)%dzdtl(1:nno)
      ap(i)%omg1(1:nno)  = glob(p)%omg1(1:nno)
      ap(i)%omg2(1:nno)  = glob(p)%omg2(1:nno)
      ap(i)%omg3(1:nno)  = glob(p)%omg3(1:nno)

      ap(i)%ua(1:nno)    = glob(p)%ua(1:nno)
      ap(i)%va(1:nno)    = glob(p)%va(1:nno)
      ap(i)%wa(1:nno)    = glob(p)%wa(1:nno)
      ap(i)%omgd1(1:nno) = glob(p)%omgd1(1:nno)
      ap(i)%omgd2(1:nno) = glob(p)%omgd2(1:nno)
      ap(i)%omgd3(1:nno) = glob(p)%omgd3(1:nno)

      ap(i)%fxll(1:nno)  = glob(p)%fxll(1:nno)
      ap(i)%fyll(1:nno)  = glob(p)%fyll(1:nno)
      ap(i)%fzll(1:nno)  = glob(p)%fzll(1:nno)
      ap(i)%mxl(1:nno)  = glob(p)%mxl(1:nno)
      ap(i)%myl(1:nno)  = glob(p)%myl(1:nno)
      ap(i)%mzl(1:nno)  = glob(p)%mzl(1:nno)

      ap(i)%ull(1:nno)   = glob(p)%ull(1:nno)
      ap(i)%vll(1:nno)   = glob(p)%vll(1:nno)
      ap(i)%wll(1:nno)   = glob(p)%wll(1:nno)
      ap(i)%omgl1(1:nno) = glob(p)%omgl1(1:nno)
      ap(i)%omgl2(1:nno) = glob(p)%omgl2(1:nno)
      ap(i)%omgl3(1:nno) = glob(p)%omgl3(1:nno)

      ap(i)%dxdto(1:nno)  = glob(p)%dxdto(1:nno)  !v I added
      ap(i)%dydto(1:nno)  = glob(p)%dydto(1:nno)
      ap(i)%dzdto(1:nno)  = glob(p)%dzdto(1:nno)
      ap(i)%omgo1(1:nno)  = glob(p)%omgo1(1:nno)
      ap(i)%omgo2(1:nno)  = glob(p)%omgo2(1:nno)
      ap(i)%omgo3(1:nno)  = glob(p)%omgo3(1:nno)

      ap(i)%uao(1:nno)    = glob(p)%uao(1:nno)    !v I added
      ap(i)%vao(1:nno)    = glob(p)%vao(1:nno)
      ap(i)%wao(1:nno)    = glob(p)%wao(1:nno)
      ap(i)%omgdo1(1:nno) = glob(p)%omgdo1(1:nno)
      ap(i)%omgdo2(1:nno) = glob(p)%omgdo2(1:nno)
      ap(i)%omgdo3(1:nno) = glob(p)%omgdo3(1:nno)
      
			do l=1,nno
        !neighbor 1
        if ( (ap(i)%xll(l)+offset) .gt. rightbound ) then
          if ( ((ap(i)%yll(l)+offset) .ge. frontbound) .and. ((ap(i)%yll(l)-offset) .le. backbound) ) then 
             ap(i)%nb(1) = 1 !neighbor 1 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 2
        if ( (ap(i)%xll(l)+offset) .gt. rightbound ) then
          if  ( (ap(i)%yll(l)-offset) .le. frontbound ) then 
            ap(i)%nb(2) = 1 !neighbor 2 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 3
        if ( (ap(i)%yll(l)-offset) .lt. frontbound ) then
          if ( ((ap(i)%xll(l)+offset) .ge. leftbound) .and. ((ap(i)%xll(l)-offset) .le. rightbound )) then 
             ap(i)%nb(3) = 1 !neighbor 3 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 4
        if ( (ap(i)%yll(l)-offset) .lt. frontbound ) then
          if  ( (ap(i)%xll(l)-offset) .le. leftbound ) then 
            ap(i)%nb(4) = 1 !neighbor 4 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 5
        if ( (ap(i)%xll(l)-offset) .lt. leftbound ) then
          if ( ((ap(i)%yll(l)+offset) .ge. frontbound) .and. ((ap(i)%yll(l)-offset) .le. backbound) ) then 
             ap(i)%nb(5) = 1 !neighbor 5 is slave of particle ap(i)%mslv
          endif
        endif 

      !neighbor 6
        if ( (ap(i)%yll(l)+offset) .gt. backbound ) then
          if  ( (ap(i)%xll(l)-offset) .le. leftbound ) then 
            ap(i)%nb(6) = 1 !neighbor 6 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 7
        if ( (ap(i)%yll(l)+offset) .gt. backbound ) then
          if ( ((ap(i)%xll(l)+offset) .ge. leftbound) .and. ((ap(i)%xll(l)-offset) .le. rightbound )) then 
             ap(i)%nb(7) = 1 !neighbor 7 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 8
        if ( (ap(i)%yll(l)+offset) .gt. backbound ) then
          if  ( (ap(i)%xll(l)+offset) .ge. rightbound ) then 
            ap(i)%nb(8) = 1 !neighbor 8 is slave of particle ap(i)%mslv
          endif
        endif
              
      enddo

    else

    count_slve_loc = 0
    !neighbor 1 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax ) + 1
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay ) 
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid .eq. procrank) then

      do l=1,nno

        if ( (glob(p)%xll(l)+offset) .gt. rightbound ) then
          if ( ((glob(p)%yll(l)+offset) .ge. frontbound) .and. ((glob(p)%yll(l)-offset) .le. backbound) ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            !ap(i)%mslv = -(((rcc-1)*(np/5))+p)     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%mslv = -(((rcc-1)*(np))+p)     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%nb(5) = 1      !neighbor 5 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 201

          endif
        endif  
      enddo

201 continue  
    endif

    !neighbor 2 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax ) + 1
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay ) - 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid .eq. procrank) then

      do l=1,nno

        if ( (glob(p)%xll(l)+offset) .gt. rightbound ) then
          if  ( (glob(p)%yll(l)-offset) .le. frontbound ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            !ap(i)%mslv = -(((rcc-1)*(np/5))+p)     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%mslv = -(((rcc-1)*(np))+p)     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%nb(6) = 1      !neighbor 6 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 202

          endif
        endif  
      enddo

202 continue 
    endif



    !neighbor 3 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax )
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay ) - 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid .eq. procrank) then


      do l=1,nno

        if ( (glob(p)%yll(l)-offset) .lt. frontbound ) then
          if ( ((glob(p)%xll(l)+offset) .ge. leftbound) .and. ((glob(p)%xll(l)-offset) .le. rightbound) ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            !ap(i)%mslv = -(((rcc-1)*(np/5))+p)     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%mslv = -(((rcc-1)*(np))+p)     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%nb(7) = 1      !neighbor 7 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 203

          endif
        endif  
      enddo

203 continue 
    endif


    !neighbor 4 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax ) - 1
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay ) - 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid .eq. procrank) then

      do l=1,nno

        if ( (glob(p)%yll(l)-offset) .lt. frontbound ) then
          if  ( (glob(p)%xll(l)-offset) .le. leftbound ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            !ap(i)%mslv = -(((rcc-1)*(np/5))+p)     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%mslv = -(((rcc-1)*(np))+p)     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%nb(8) = 1      !neighbor 8 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 204

          endif
        endif  
      enddo

204 continue 
    endif


    !neighbor 5 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax ) - 1
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay )
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid .eq. procrank) then

      do l=1,nno

      if ( (glob(p)%xll(l)-offset) .lt. leftbound ) then
        if ( ((glob(p)%yll(l)+offset) .ge. frontbound) .and. ((glob(p)%yll(l)-offset) .le. backbound ) ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            !ap(i)%mslv = -(((rcc-1)*(np/5))+p)     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%mslv = -(((rcc-1)*(np))+p)     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%nb(1) = 1      !neighbor 1 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 205

          endif
        endif  
      enddo

205 continue 
    endif


    !neighbor 6 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax ) - 1
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay ) + 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid .eq. procrank) then

      do l=1,nno

      if ( (glob(p)%yll(l)+offset) .gt. backbound ) then
        if  ( (glob(p)%xll(l)-offset) .le. leftbound ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            !ap(i)%mslv = -(((rcc-1)*(np/5))+p)     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%mslv = -(((rcc-1)*(np))+p)     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%nb(2) = 1      !neighbor 2 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 206

          endif
        endif  
      enddo

206 continue 
    endif


    !neighbor 7 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax )
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay ) + 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid .eq. procrank) then

      do l=1,nno

      if ( (glob(p)%yll(l)+offset) .gt. backbound ) then
        if ( ((glob(p)%xll(l)+offset) .ge. leftbound) .and. ((glob(p)%xll(l)-offset) .le. rightbound) ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            !ap(i)%mslv = -(((rcc-1)*(np/5))+p)     !myid is slave of particle p=ap(i)%mslv
            ap(i)%mslv = -(((rcc-1)*(np))+p)     !myid is slave of particle p=ap(i)%mslv
            ap(i)%nb(3) = 1      !neighbor 3 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 207

          endif
        endif  
      enddo

207 continue 
    endif


    !neighbor 8 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax ) + 1
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay ) + 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid .eq. procrank) then

      do l=1,nno

      if ( (glob(p)%yll(l)+offset) .gt. backbound ) then
        if  ( (glob(p)%xll(l)+offset) .ge. rightbound ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            !ap(i)%mslv = -(((rcc-1)*(np/5))+p)     !myid is slave of particle p=ap(i)%mslv
            ap(i)%mslv = -(((rcc-1)*(np))+p)     !myid is slave of particle p=ap(i)%mslv
            ap(i)%nb(4) = 1      !neighbor 4 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 208

          endif
        endif  
      enddo

208 continue 
    endif

  endif
enddo


if (myid==0) then
 close(20)
endif

enddo!do rcc
  
pmax = i 

  !if(pmax.eq.0) pmax = 1 ! for now
  print*, 'Thread ', myid, ' contains ', pmax, 'paricles.'

  call MPI_ALLREDUCE(npmstr,count_mstr_all,1,MPI_INTEGER,MPI_SUM,comm_cart,ierr)
  print*,count_mstr_all,npmstr,np
  if (count_mstr_all.ne.np) then
    write(6,*) 'Fatal error in initialisation of particle positions!'
    write(6,*) 'Program aborted...'
    call mpi_abort(comm_cart,ierr,ierr)
    stop
  elseif(pmax.gt.npmax) then
    write(6,*) 'Size of local particle array of process ', myid, ' is too small!'
    write(6,*) 'Program aborted... (later I will simply write a warning and allocate more memory)'
    call mpi_abort(comm_cart,ierr,ierr)
    stop
  else
    write(6,*) 'The particles were successfully initialized in thread ', myid, ' !'
  endif

endif ! in.eq.0
!










if (in.eq.1) then
!!
!! write particle related data directly in parallel with MPI-IO


!do rcc=1,5
do rcc=1,1

write(fl,'(i1.1)') rcc


!if (myid==0) then
! !open(20,file=datadir//'allpart'//fl//'data'//istepchar,access='direct',recl=np/5*skip*lenr)
! open(20,file=datadir//'allpart'//fl//'data'//istepchar,access='direct',recl=np/1*skip*lenr)
!endif

  !do i=1,np/5
  do i=1,np
   !pcc=((np/5)*(rcc-1))+i
   pcc=(np*(rcc-1))+i
   counter = 0
    do p=1,pmax
      if (pcc.eq.ap(p)%mslv) then
        counter = counter + 1
        
				glob(i)%x = ap(p)%x !WRITE THIS IN A MORE COMPACT WAY!
        glob(i)%y = ap(p)%y
        glob(i)%z = ap(p)%z

        glob(i)%integralx = ap(p)%integralx
        glob(i)%integraly = ap(p)%integraly
        glob(i)%integralz = ap(p)%integralz

				glob(i)%xll(1:nno) = ap(p)%xll(1:nno)
				glob(i)%yll(1:nno) = ap(p)%yll(1:nno)
      	glob(i)%zll(1:nno) = ap(p)%zll(1:nno)
      	glob(i)%qn1(1:nno) = ap(p)%qn1(1:nno)
      	glob(i)%qn2(1:nno) = ap(p)%qn2(1:nno)
      	glob(i)%qn3(1:nno) = ap(p)%qn3(1:nno)
      	glob(i)%qn4(1:nno) = ap(p)%qn4(1:nno)
			
				!print*, "---value of ap%qn1 from saved data"
				!do icp =1,nno
				!!print*, "---ap(i)%qN1(1:nno)", real(ap(p)%qn1(icp),4)
				!print*, "---glob(i)%qN1(1:nno)", real(glob(i)%qn1(icp),4)
				!enddo
				!print*, "-------pmax", pmax
				!print*, "---value of ap%qn1 from saved data"

      	glob(i)%xfpo(1:nno) = ap(p)%xfpo(1:nno)
      	glob(i)%yfpo(1:nno) = ap(p)%yfpo(1:nno)
      	glob(i)%zfpo(1:nno) = ap(p)%zfpo(1:nno)
      	glob(i)%qno1(1:nno) = ap(p)%qno1(1:nno)
      	glob(i)%qno2(1:nno) = ap(p)%qno2(1:nno)
      	glob(i)%qno3(1:nno) = ap(p)%qno3(1:nno)
      	glob(i)%qno4(1:nno) = ap(p)%qno4(1:nno)

      	glob(i)%xfpold(1:nno) = ap(p)%xfpold(1:nno)
      	glob(i)%yfpold(1:nno) = ap(p)%yfpold(1:nno)
      	glob(i)%zfpold(1:nno) = ap(p)%zfpold(1:nno)
      	glob(i)%qnold1(1:nno) = ap(p)%qnold1(1:nno)
      	glob(i)%qnold2(1:nno) = ap(p)%qnold2(1:nno)
      	glob(i)%qnold3(1:nno) = ap(p)%qnold3(1:nno)
      	glob(i)%qnold4(1:nno) = ap(p)%qnold4(1:nno)

      	glob(i)%dxdtl(1:nno) = ap(p)%dxdtl(1:nno)
      	glob(i)%dydtl(1:nno) = ap(p)%dydtl(1:nno)
      	glob(i)%dzdtl(1:nno) = ap(p)%dzdtl(1:nno)
      	glob(i)%omg1(1:nno) = ap(p)%omg1(1:nno)
      	glob(i)%omg2(1:nno) = ap(p)%omg2(1:nno)
      	glob(i)%omg3(1:nno) = ap(p)%omg3(1:nno)

      	glob(i)%ua(1:nno)   = ap(p)%ua(1:nno)
      	glob(i)%va(1:nno)   = ap(p)%va(1:nno)
      	glob(i)%wa(1:nno)   = ap(p)%wa(1:nno)
      	glob(i)%omgd1(1:nno)= ap(p)%omgd1(1:nno)
      	glob(i)%omgd2(1:nno)= ap(p)%omgd2(1:nno)
      	glob(i)%omgd3(1:nno)= ap(p)%omgd3(1:nno)

				glob(i)%fxll(1:nno)  = ap(p)%fxll(1:nno)
      	glob(i)%fyll(1:nno)  = ap(p)%fyll(1:nno)
      	glob(i)%fzll(1:nno)  = ap(p)%fzll(1:nno)
      	glob(i)%mxl(1:nno)  = ap(p)%mxl(1:nno)
      	glob(i)%myl(1:nno)  = ap(p)%myl(1:nno)
      	glob(i)%mzl(1:nno)  = ap(p)%mzl(1:nno)

      	glob(i)%ull(1:nno)  = ap(p)%ull(1:nno)
      	glob(i)%vll(1:nno)  = ap(p)%vll(1:nno)
      	glob(i)%wll(1:nno)  = ap(p)%wll(1:nno)
      	glob(i)%omgl1(1:nno)= ap(p)%omgl1(1:nno)
      	glob(i)%omgl2(1:nno)= ap(p)%omgl2(1:nno)
      	glob(i)%omgl3(1:nno)= ap(p)%omgl3(1:nno)

      	glob(i)%dxdto(1:nno) = ap(p)%dxdto(1:nno)  !v I added
      	glob(i)%dydto(1:nno) = ap(p)%dydto(1:nno)
      	glob(i)%dzdto(1:nno) = ap(p)%dzdto(1:nno)
      	glob(i)%omgo1(1:nno) = ap(p)%omgo1(1:nno)
      	glob(i)%omgo2(1:nno) = ap(p)%omgo2(1:nno)
      	glob(i)%omgo3(1:nno) = ap(p)%omgo3(1:nno)

      	glob(i)%uao(1:nno)   = ap(p)%uao(1:nno)  !v I added
      	glob(i)%vao(1:nno)   = ap(p)%vao(1:nno)
      	glob(i)%wao(1:nno)   = ap(p)%wao(1:nno)
      	glob(i)%omgdo1(1:nno)= ap(p)%omgdo1(1:nno)
      	glob(i)%omgdo2(1:nno)= ap(p)%omgdo2(1:nno)
      	glob(i)%omgdo3(1:nno)= ap(p)%omgdo3(1:nno)
      	
      endif
    enddo

    if(counter.eq.0) then
      glob(i)%x = 0. !WRITE THIS IN A MORE COMPACT WAY!
      glob(i)%y = 0.
      glob(i)%z = 0.

      glob(i)%integralx = 0.
      glob(i)%integraly = 0.
      glob(i)%integralz = 0.

			glob(i)%xll(1:nno) = 0.
      glob(i)%yll(1:nno) = 0.
      glob(i)%zll(1:nno) = 0.
      glob(i)%qn1(1:nno) = 0.
      glob(i)%qn2(1:nno) = 0.
      glob(i)%qn3(1:nno) = 0.
      glob(i)%qn4(1:nno) = 0.

      glob(i)%xfpo(1:nno) = 0.
      glob(i)%yfpo(1:nno) = 0.
      glob(i)%zfpo(1:nno) = 0.
      glob(i)%qno1(1:nno) = 0.
      glob(i)%qno2(1:nno) = 0.
      glob(i)%qno3(1:nno) = 0.
      glob(i)%qno4(1:nno) = 0.

      glob(i)%xfpold(1:nno) = 0.
      glob(i)%yfpold(1:nno) = 0.
      glob(i)%zfpold(1:nno) = 0.
      glob(i)%qnold1(1:nno) = 0.
      glob(i)%qnold2(1:nno) = 0.
      glob(i)%qnold3(1:nno) = 0.
      glob(i)%qnold4(1:nno) = 0.

      glob(i)%dxdtl(1:nno) = 0.
      glob(i)%dydtl(1:nno) = 0.
      glob(i)%dzdtl(1:nno) = 0.
      glob(i)%omg1(1:nno) = 0.
      glob(i)%omg2(1:nno) = 0.
      glob(i)%omg3(1:nno) = 0.

      glob(i)%ua(1:nno)   = 0.
      glob(i)%va(1:nno)   = 0.
      glob(i)%wa(1:nno)   = 0.
      glob(i)%omgd1(1:nno)= 0.
      glob(i)%omgd2(1:nno)= 0.
      glob(i)%omgd3(1:nno)= 0.

      glob(i)%fxll(1:nno)  = 0.
      glob(i)%fyll(1:nno)  = 0.
      glob(i)%fzll(1:nno)  = 0.
      glob(i)%mxl(1:nno)  = 0.
      glob(i)%myl(1:nno)  = 0.
      glob(i)%mzl(1:nno)  = 0.

      glob(i)%ull(1:nno)    = 0.
      glob(i)%vll(1:nno)    = 0.
      glob(i)%wll(1:nno)    = 0.
      glob(i)%omgl1(1:nno) = 0.
      glob(i)%omgl2(1:nno) = 0.
      glob(i)%omgl3(1:nno) = 0.

      glob(i)%dxdto(1:nno) = 0.  !v I added
      glob(i)%dydto(1:nno) = 0.
      glob(i)%dzdto(1:nno) = 0.
      glob(i)%omgo1(1:nno) = 0.
      glob(i)%omgo2(1:nno) = 0.
      glob(i)%omgo3(1:nno) = 0.

      glob(i)%uao(1:nno)   = 0.  !v I added
      glob(i)%vao(1:nno)   = 0.
      glob(i)%wao(1:nno)   = 0.
      glob(i)%omgdo1(1:nno)= 0.
      glob(i)%omgdo2(1:nno)= 0.
      glob(i)%omgdo3(1:nno)= 0.

    endif
  enddo
  !call mpi_reduce(glob(1)%x,glob_all(1)%x,np/5*skip,mpi_real8,mpi_sum,0,MPI_COMM_WORLD,ierr)
  call mpi_reduce(glob(1)%x,glob_all(1)%x,np*skip,mpi_real8,mpi_sum,0,MPI_COMM_WORLD,ierr)
  if (myid.eq.0) then ! myid = 0 writes the data into a single file
		 open(20,file=datadir//'allpart'//fl//'data'//istepchar,access='direct',recl=np*skip*lenr)
     write(20,rec=rcc) glob_all
		!close(20)

		!print*, "---value of glob_all%qn1 after saving data"
    !do icp =1,nno
    !!print*, "---ap(i)%qN1(1:nno)", real(ap(p)%qn1(icp),4)
    !print*, "---glob(i)%qN1(1:nno)", real(glob_all%qn1(icp),4)
    !!print*, "---glob(i)%qN1(1:nno)", real(glob%qn1(icp),4)
    !enddo
		!print*, "---value of glob_all%qn1 after saving data"

  endif

if (myid==0) then
 close(20)
endif

enddo !do rcc
endif ! if (in .eq. 1)



return
end subroutine loadpart
!
end module mod_loadpart
