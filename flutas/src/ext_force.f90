module mod_ext_force
use mod_param_fibm,    only: nl, i1, j1, k1, ds, aspectratio, filsub, &
													alphafilament, npmax, nxi_vec, nno
use mod_common_fibm,   only: betafilament, forcex, forcey, forcez, lcounter, &
													pmax, ap !forcex, forcey, forcez
use mod_interp_spread, only: eulr2lagr, lagr2eulr
use mod_initialSetup,  only: ap_ugp, ap_fgp
!use mod_forcing
!use mod_common_mpi,   only: MPI_REAL8, &
!                            MPI_SUM, status, xhalo, yhalo, MPI_STATUS_SIZE, &
!                            neighbor
implicit none
private
public ext_force
contains
!
!subroutine ext_force(unew, vnew, wnew, dt, istep, time)
subroutine ext_force(unew, vnew, wnew, dt, istep, time, forcex, forcey, forcez)
implicit none
real(8), intent(in), dimension(0:i1, 0:j1, 0:k1)  :: unew, vnew, wnew
real(8), intent(in)                               :: dt,time
integer, intent(in)                               :: istep
real, intent(inout), dimension(0:i1, 0:j1, 0:k1)  :: forcex, forcey, forcez

integer :: l,p,i,j,k,ii
real :: dti
real,dimension(NL)::dxdsf,dydsf,dzdsf,tx,ty,tz
real,dimension(NL)::txs,tys,tzs
real,dimension(NL-1)::xs,ys,zs
real::xn1,xn2,xn3,xn4,yn1,yn2,yn3,yn4,zn1,zn2,zn3,zn4
real::ftx1,fty1,ftz1,ftx2,fty2,ftz2
real::ftx3,fty3,ftz3,ftx4,fty4,ftz4
real::nx,ny,nz,kr,frx,fry,frz
!
dti = 1./dt

!print*,"=-=-=-=-=-=-=-=-= dt inside ext_force:     ", dt

betafilament= -1./dt
!betafilament= -1000


!betafilament=betafilament/filsub
!
!$omp parallel default(shared) &
!$omp& private(p,l)
!$omp do


!v initialized in interp_spread. 
!v forcex(:,:,:) = 0.  !v size 0:i1,0:j1,0:k1, i1 -> max number grid points per proc
!v forcey(:,:,:) = 0.
!v forcez(:,:,:) = 0.


!--------- interp velocity at specified parameteric points
call ap_ugp ! out: ap(pp)%xfp, ap(pp)%dxdt


if (lcounter==1) then  !v lcounter==1 herein     

  do p=1,pmax               !v  herein, computing coords of 4 ghost p (xn1-xn4) and center point (xtc) upto L114
                            !v pmax: ! maximum number of particles in a thread is equal to the number of particles 
                            !          'mastered' and 'enslaved' by it

    if (ap(p)%mslv>0) then

      !do l=1,nno            !v dxdsf etc are used to compute xn1-xn4 etc
      !  if (l==1) then      !v dxdsf are tangents at each node
      !    dxdsf(l)=(ap(p)%xll(l+1)-ap(p)%xll(l))/ds
      !    dydsf(l)=(ap(p)%yll(l+1)-ap(p)%yll(l))/ds
      !    dzdsf(l)=(ap(p)%zll(l+1)-ap(p)%zll(l))/ds
      !  elseif (l==NL) then
      !    dxdsf(l)=(ap(p)%xll(l)-ap(p)%xll(l-1))/ds
      !    dydsf(l)=(ap(p)%yll(l)-ap(p)%yll(l-1))/ds
      !    dzdsf(l)=(ap(p)%zll(l)-ap(p)%zll(l-1))/ds
      !  else
      !    dxdsf(l)=(ap(p)%xll(l+1)-ap(p)%xll(l-1))/(2.*ds) 
      !    dydsf(l)=(ap(p)%yll(l+1)-ap(p)%yll(l-1))/(2.*ds)
      !    dzdsf(l)=(ap(p)%zll(l+1)-ap(p)%zll(l-1))/(2.*ds)
      !  endif
      !enddo

      do l=1,nno
       ! if (dxdsf(l)==0.) then
       !   xn1=1.
       !   yn1=0.
       !   zn1=0.
       ! else
       !   yn1=sqrt(1./(1.+((dydsf(l)**2.)/(dxdsf(l)**2.))))
       !   xn1=-(yn1*dydsf(l))/dxdsf(l)
       !   zn1=0.
       ! endif

       ! xn2=(zn1*dydsf(l))-(yn1*dzdsf(l))
       ! yn2=(xn1*dzdsf(l))-(zn1*dxdsf(l))
       ! zn2=(yn1*dxdsf(l))-(xn1*dydsf(l))

       ! xn3=(zn2*dydsf(l))-(yn2*dzdsf(l))
       ! yn3=(xn2*dzdsf(l))-(zn2*dxdsf(l))
       ! zn3=(yn2*dxdsf(l))-(xn2*dydsf(l))

       ! xn4=(zn3*dydsf(l))-(yn3*dzdsf(l))
       ! yn4=(xn3*dzdsf(l))-(zn3*dxdsf(l))
       ! zn4=(yn3*dxdsf(l))-(xn3*dydsf(l))


       ! ap(p)%xt1(l)=ap(p)%xll(l)+(aspectratio*.5*xn1) !v these are finally the coordinates of four ghost points, from xn1=xn4 etc
       ! ap(p)%yt1(l)=ap(p)%yll(l)+(aspectratio*.5*yn1)
       ! ap(p)%zt1(l)=ap(p)%zll(l)+(aspectratio*.5*zn1)

       ! ap(p)%xt2(l)=ap(p)%xll(l)+(aspectratio*.5*xn2)
       ! ap(p)%yt2(l)=ap(p)%yll(l)+(aspectratio*.5*yn2)
       ! ap(p)%zt2(l)=ap(p)%zll(l)+(aspectratio*.5*zn2)

       ! ap(p)%xt3(l)=ap(p)%xll(l)+(aspectratio*.5*xn3)
       ! ap(p)%yt3(l)=ap(p)%yll(l)+(aspectratio*.5*yn3)
       ! ap(p)%zt3(l)=ap(p)%zll(l)+(aspectratio*.5*zn3)

       ! ap(p)%xt4(l)=ap(p)%xll(l)+(aspectratio*.5*xn4)
       ! ap(p)%yt4(l)=ap(p)%yll(l)+(aspectratio*.5*yn4)
       ! ap(p)%zt4(l)=ap(p)%zll(l)+(aspectratio*.5*zn4)

        ap(p)%xtc(l)=ap(p)%xll(l)    
        ap(p)%ytc(l)=ap(p)%yll(l)
        ap(p)%ztc(l)=ap(p)%zll(l)

      enddo  !do l=1,nl
    endif    !if (ap(p)%mslv>0)
  enddo      !do p=1,pmax



  !v initialized: Lagrangian forces and forces at 4 GPs, for max number of particles per processor
  if (lcounter==1) then
    do p=1,npmax
			ap(p)%fxll(:)=0.
      ap(p)%fyll(:)=0.
      ap(p)%fzll(:)=0.

      ap(p)%fxl(:)=0.
      ap(p)%fyl(:)=0.
      ap(p)%fzl(:)=0.
      ap(p)%mxl(:)=0.
      ap(p)%myl(:)=0.
      ap(p)%mzl(:)=0.

      ap(p)%fx1(:)=0. !v maybe, remove this part
      ap(p)%fy1(:)=0.
      ap(p)%fz1(:)=0.
      ap(p)%fx2(:)=0.
      ap(p)%fy2(:)=0.
      ap(p)%fz2(:)=0.
      ap(p)%fx3(:)=0.
      ap(p)%fy3(:)=0.
      ap(p)%fz3(:)=0.
      ap(p)%fx4(:)=0.
      ap(p)%fy4(:)=0.
      ap(p)%fz4(:)=0.
    enddo
  endif




  !do p=1,pmax
  !if (ap(p)%mslv>0) then
  !  do l=1,NL
  !    ap(p)%dxdt(l)=(ap(p)%xfpo(l) - ap(p)%xfpold(l)) / dt
  !    ap(p)%dydt(l)=(ap(p)%yfpo(l) - ap(p)%yfpold(l)) / dt
  !    ap(p)%dzdt(l)=(ap(p)%zfpo(l) - ap(p)%zfpold(l)) / dt
  !  enddo
  !endif
  !enddo


  !do p=1,pmax
  !  if (ap(p)%mslv>0) then
  !    ap(p)%xfp(:)=ap(p)%xt1(:)
  !    ap(p)%yfp(:)=ap(p)%yt1(:)
  !    ap(p)%zfp(:)=ap(p)%zt1(:)
  !  endif
  !enddo
  !!v interpolating fluid vel from Eulerian to Lagrangian points
  !!v here, ap(p)%xfp etc are used to compute U_ib eventually
  !call eulr2lagr(unew,vnew,wnew)
  !  
  !!v here, Computed Lagrangian forces at fiber points using Eq. (3.32)
  !do p=1,pmax
  !  if (ap(p)%mslv>0) then
  !    do l=1,NL        
  !      ap(p)%fx1(l)=((betafilament)*( ap(p)%ul(l) - ap(p)%dxdt(l) ))*.25
  !      ap(p)%fy1(l)=((betafilament)*( ap(p)%vl(l) - ap(p)%dydt(l) ))*.25
  !      ap(p)%fz1(l)=((betafilament)*( ap(p)%wl(l) - ap(p)%dzdt(l) ))*.25

  !      !ap(p)%fxl(l)=ap(p)%fx1(l)
  !      !ap(p)%fyl(l)=ap(p)%fy1(l)
  !      !ap(p)%fzl(l)=ap(p)%fz1(l)
  !    enddo
  !  endif
  !enddo
  !!call lagr2eulr

  
	!do p=1,pmax
  !  if (ap(p)%mslv>0) then
  !    ap(p)%xfp(:)=ap(p)%xt2(:)
  !    ap(p)%yfp(:)=ap(p)%yt2(:)
  !    ap(p)%zfp(:)=ap(p)%zt2(:)
  !  endif
  !enddo
  !call eulr2lagr(unew,vnew,wnew)
  !do p=1,pmax
  !  if (ap(p)%mslv>0) then
  !    do l=1,NL
  !      ap(p)%fx2(l)=((betafilament)*(ap(p)%ul(l) - ap(p)%dxdt(l)))*.25
  !      ap(p)%fy2(l)=((betafilament)*(ap(p)%vl(l) - ap(p)%dydt(l)))*.25
  !      ap(p)%fz2(l)=((betafilament)*(ap(p)%wl(l) - ap(p)%dzdt(l)))*.25
  !      !ap(p)%fxl(l)=ap(p)%fx2(l)
  !      !ap(p)%fyl(l)=ap(p)%fy2(l)
  !      !ap(p)%fzl(l)=ap(p)%fz2(l)
  !    enddo
  !  endif
  !enddo


  !do p=1,pmax
  !  if (ap(p)%mslv>0) then
  !    ap(p)%xfp(:)=ap(p)%xt3(:)
  !    ap(p)%yfp(:)=ap(p)%yt3(:)
  !    ap(p)%zfp(:)=ap(p)%zt3(:)
  !  endif
  !enddo
  !call eulr2lagr(unew,vnew,wnew)
  !do p=1,pmax
  !  if (ap(p)%mslv>0) then
  !    do l=1,NL
  !      ap(p)%fx3(l)=((betafilament)*(ap(p)%ul(l) - ap(p)%dxdt(l)))*.25
  !      ap(p)%fy3(l)=((betafilament)*(ap(p)%vl(l) - ap(p)%dydt(l)))*.25
  !      ap(p)%fz3(l)=((betafilament)*(ap(p)%wl(l) - ap(p)%dzdt(l)))*.25
  !      !ap(p)%fxl(l)=ap(p)%fx3(l)
  !      !ap(p)%fyl(l)=ap(p)%fy3(l)
  !      !ap(p)%fzl(l)=ap(p)%fz3(l)
  !    enddo
  !  endif
  !enddo



  !do p=1,pmax
  !  if (ap(p)%mslv>0) then
  !    ap(p)%xfp(:)=ap(p)%xt4(:)
  !    ap(p)%yfp(:)=ap(p)%yt4(:)
  !    ap(p)%zfp(:)=ap(p)%zt4(:)
  !  endif
  !enddo
  !call eulr2lagr(unew,vnew,wnew)
  !do p=1,pmax
  !  if (ap(p)%mslv>0) then
  !    do l=1,NL
  !      ap(p)%fx4(l)=((betafilament)*(ap(p)%ul(l) - ap(p)%dxdt(l)))*.25
  !      ap(p)%fy4(l)=((betafilament)*(ap(p)%vl(l) - ap(p)%dydt(l)))*.25
  !      ap(p)%fz4(l)=((betafilament)*(ap(p)%wl(l) - ap(p)%dzdt(l)))*.25
  !      !ap(p)%fxl(l)=ap(p)%fx4(l)
  !      !ap(p)%fyl(l)=ap(p)%fy4(l)
  !      !ap(p)%fzl(l)=ap(p)%fz4(l)
  !    enddo      
  !  endif
  !enddo




  do p=1,pmax
    if (ap(p)%mslv>0) then
      do l=1,nno
        ap(p)%xll(l)=ap(p)%xtc(l)
        ap(p)%yll(l)=ap(p)%ytc(l)
        ap(p)%zll(l)=ap(p)%ztc(l)
      enddo
    endif
  enddo

	!v I: ap(p)%xfp, O: ap(p)%ul
  call eulr2lagr(unew,vnew,wnew)
endif 





ap(1:npmax)%integralx = 0.
ap(1:npmax)%integraly = 0.
ap(1:npmax)%integralz = 0.

!do p=1,pmax
!  if (ap(p)%mslv>0) then
!
!    do l=1,NL
!      !ap(p)%fxl(l)= ap(p)%fxl(l)+((betafilament)*(ap(p)%ul(l)-ap(p)%dxdt(l)))!-ap(p)%fcx(l)
!      !ap(p)%fyl(l)= ap(p)%fyl(l)+((betafilament)*(ap(p)%vl(l)-ap(p)%dydt(l)))!-ap(p)%fcy(l)
!      !ap(p)%fzl(l)= ap(p)%fzl(l)+((betafilament)*(ap(p)%wl(l)-ap(p)%dzdt(l)))!-ap(p)%fcz(l)
!      
!			!v ===== commented to try feedback law based ibm expression for force
!			ap(p)%fxl(l)= ((betafilament)*(ap(p)%ul(l) - ap(p)%dxdt(l))) !-ap(p)%fcx(l)      !v -1* Eq. (3.32)
!      ap(p)%fyl(l)= ((betafilament)*(ap(p)%vl(l) - ap(p)%dydt(l))) !-ap(p)%fcy(l)
!      ap(p)%fzl(l)= ((betafilament)*(ap(p)%wl(l) - ap(p)%dzdt(l))) !-ap(p)%fcz(l)
!
!			!v F_ib based on feedback law 
!			!ap(p)%integralx=ap(p)%integralx+(alphafilament*dt*(ap(p)%ul(l)-ap(p)%dxdt(l)))
!			!ap(p)%integraly=ap(p)%integraly+(alphafilament*dt*(ap(p)%vl(l)-ap(p)%dydt(l)))
!			!ap(p)%integralz=ap(p)%integralz+(alphafilament*dt*(ap(p)%wl(l)-ap(p)%dzdt(l)))
!
!			!ap(p)%fxl(l)= ap(p)%integralx+(betafilament*(ap(p)%ul(l)-ap(p)%dxdt(l))) 
!			!ap(p)%fyl(l)= ap(p)%integraly+(betafilament*(ap(p)%vl(l)-ap(p)%dydt(l))) 
!			!ap(p)%fzl(l)= ap(p)%integralz+(betafilament*(ap(p)%wl(l)-ap(p)%dzdt(l))) 
!
!      !print '(3E12.4)', ap(p)%ul(l), ap(p)%vl(l), ap(p)%wl(l)
!      !print '(3E12.4)', ap(p)%integralx, ap(p)%integraly, ap(p)%integralz
!      !print '(3E12.4)', ap(p)%fxl(l), ap(p)%fyl(l), ap(p)%fzl(l)
!    enddo
!
!    !  do l=1,NL  
!    !    ap(p)%fxl(l)=ap(p)%fx1(l)+ap(p)%fx2(l)+ap(p)%fx3(l)+ap(p)%fx4(l)
!    !    ap(p)%fyl(l)=ap(p)%fy1(l)+ap(p)%fy2(l)+ap(p)%fy3(l)+ap(p)%fy4(l)
!    !    ap(p)%fzl(l)=ap(p)%fz1(l)+ap(p)%fz2(l)+ap(p)%fz3(l)+ap(p)%fz4(l)
!    !  enddo
!		!print*,'=================Forces=========================='
!    !do ii = 1,NL
!    !  print '(3E12.4)', ap(p)%fxl(ii), ap(p)%fyl(ii), ap(p)%fzl(ii)!, ap(p)%mxl(ii), ap(p)%myl(ii), ap(p)%mzl(ii)
!    !enddo
!    ! print*,'=================Fluid vel=========================='
!    ! do ii = 1,NL
!    !   print '(3E12.4)', ap(p)%ul(ii), ap(p)%vl(ii), ap(p)%wl(ii)!, ap(p)%mxl(ii), ap(p)%myl(ii), ap(p)%mzl(ii)
!    ! enddo
!  endif
!enddo



do p=1,pmax
  if (ap(p)%mslv>0) then

    do l=1,nl
			ap(p)%fxl(l)= ((betafilament)*(ap(p)%ul(l) - ap(p)%dxdt(l)))
      ap(p)%fyl(l)= ((betafilament)*(ap(p)%vl(l) - ap(p)%dydt(l)))
      ap(p)%fzl(l)= ((betafilament)*(ap(p)%wl(l) - ap(p)%dzdt(l)))
		enddo	
		!print*,'=================Forces at xiP=========================='
    !do ii = 1,nl
    !  print '(3E12.4)', ap(p)%fxl(ii), ap(p)%fyl(ii), ap(p)%fzl(ii)!, ap(p)%mxl(ii), ap(p)%myl(ii), ap(p)%mzl(ii)
    !enddo
    !print*,'=================Fluid vel=========================='
    !do ii = 1,nl
    !  print '(3E12.4)', ap(p)%ul(ii), ap(p)%vl(ii), ap(p)%wl(ii)!, ap(p)%mxl(ii), ap(p)%myl(ii), ap(p)%mzl(ii)
    !enddo
    !print*,'=================LagP vel=========================='
    !do ii = 1,nl
    !  print '(3E12.4)', ap(p)%dxdt(ii), ap(p)%dydt(ii), ap(p)%dzdt(ii)!, ap(p)%mxl(ii), ap(p)%myl(ii), ap(p)%mzl(ii)
    !enddo

	endif !mslv
enddo	!pmax
		




if (lcounter==filsub) then
   call lagr2eulr(istep, time, forcex, forcey, forcez)  
endif

do p=1,pmax
  if (ap(p)%mslv>0) then
    do l=1,nl
      ap(p)%fxl(l)= ((betafilament)*(ap(p)%ul(l) - ap(p)%dxdt(l)))
      ap(p)%fyl(l)= ((betafilament)*(ap(p)%vl(l) - ap(p)%dydt(l)))
      ap(p)%fzl(l)= ((betafilament)*(ap(p)%wl(l) - ap(p)%dzdt(l)))
    enddo
  endif !mslv
enddo !pmax


!I: ap(p)%fxl, O: ap(p)%fxll
call ap_fgp 


!do p=1,pmax
	!if (ap(p)%mslv>0) then
		!do l=1,NL
			!ap(p)%fxl(l)=ap(p)%fxl(l)+ap(p)%fcx(l)
			!ap(p)%fyl(l)=ap(p)%fyl(l)+ap(p)%fcy(l)
			!ap(p)%fzl(l)=ap(p)%fzl(l)+ap(p)%fcz(l)
		!enddo
	!endif
!enddo




!v computing torque/moment exerted by the hydrodynamic forces on fiber given by Eq. (3.36-3.37) ?
!do p=1,pmax
!  if (ap(p)%mslv>0) then
!
!    do l=1,nl-1
!      xs(l)=(ap(p)%xfp(l+1)-ap(p)%xfp(l))/ds    !v not being used anywhere, maybe implemented for dubugging. 
!      ys(l)=(ap(p)%yfp(l+1)-ap(p)%yfp(l))/ds    !v perhaps, checking if xs==1 or not
!      zs(l)=(ap(p)%zfp(l+1)-ap(p)%zfp(l))/ds
!    enddo
!
!    do l=1,NL   
!      ftx1=(-ap(p)%fx1(l)*abs(ap(p)%dxds(l)))   !v compute forces * tangents along each direction 
!      fty1=(-ap(p)%fy1(l)*abs(ap(p)%dyds(l)))
!      ftz1=(-ap(p)%fz1(l)*abs(ap(p)%dzds(l)))
!      ftx2=(-ap(p)%fx2(l)*abs(ap(p)%dxds(l)))
!      fty2=(-ap(p)%fy2(l)*abs(ap(p)%dyds(l)))
!      ftz2=(-ap(p)%fz2(l)*abs(ap(p)%dzds(l)))
!      ftx3=(-ap(p)%fx3(l)*abs(ap(p)%dxds(l)))
!      fty3=(-ap(p)%fy3(l)*abs(ap(p)%dyds(l)))
!      ftz3=(-ap(p)%fz3(l)*abs(ap(p)%dzds(l)))
!      ftx4=(-ap(p)%fx4(l)*abs(ap(p)%dxds(l)))
!      fty4=(-ap(p)%fy4(l)*abs(ap(p)%dyds(l)))
!      ftz4=(-ap(p)%fz4(l)*abs(ap(p)%dzds(l)))
!
!      ! if (l==1.or.l==nl) then
!
!      !   !v this one setting zero moments at the end points, M = -r x F
!      !   tx(l)=-( (fty1*.25*(ap(p)%ztc(l)-ap(p)%zt1(l))) + (fty2*.25*(ap(p)%ztc(l)-ap(p)%zt2(l))) + &
!      !   (fty3*.25*(ap(p)%ztc(l)-ap(p)%zt3(l)))+(fty4*.25*(ap(p)%ztc(l)-ap(p)%zt4(l))) + &
!      !   (ftz1*.25*(ap(p)%yt1(l)-ap(p)%ytc(l)))+(ftz2*.25*(ap(p)%yt2(l)-ap(p)%ytc(l))) + &
!      !   (ftz3*.25*(ap(p)%yt3(l)-ap(p)%ytc(l)))+(ftz4*.25*(ap(p)%yt4(l)-ap(p)%ytc(l))) ) !v -1*(fty*zt - ftz*yt) for all 4 points
!
!      !   ty(l)=-((ftx1*.25*(ap(p)%zt1(l)-ap(p)%ztc(l)))+(ftx2*.25*(ap(p)%zt2(l)-ap(p)%ztc(l))) + &
!      !   (ftx3*.25*(ap(p)%zt3(l)-ap(p)%ztc(l)))+(ftx4*.25*(ap(p)%zt4(l)-ap(p)%ztc(l))) + &
!      !   (ftz1*.25*(ap(p)%xtc(l)-ap(p)%xt1(l)))+(ftz2*.25*(ap(p)%xtc(l)-ap(p)%xt2(l))) + &
!      !   (ftz3*.25*(ap(p)%xtc(l)-ap(p)%xt3(l)))+(ftz4*.25*(ap(p)%xtc(l)-ap(p)%xt4(l)))) !v -1*(ftx*zt - ftz*xt) for all 4 points
!
!      !   tz(l)=-((ftx1*.25*(ap(p)%ytc(l)-ap(p)%yt1(l)))+(ftx2*.25*(ap(p)%ytc(l)-ap(p)%yt2(l)))+ &
!      !   (ftx3*.25*(ap(p)%ytc(l)-ap(p)%yt3(l)))+(ftx4*.25*(ap(p)%ytc(l)-ap(p)%yt4(l))) + &
!      !   (fty1*.25*(ap(p)%xt1(l)-ap(p)%xtc(l)))+(fty2*.25*(ap(p)%xt2(l)-ap(p)%xtc(l))) + &
!      !   (fty3*.25*(ap(p)%xt3(l)-ap(p)%xtc(l)))+(fty4*.25*(ap(p)%xt4(l)-ap(p)%xtc(l)))) !v -1*(ftx*zt - ftz*xt) for all 4 points
!
!      ! else
!
!       !v tx(l)=(fty1*.25*(ap(p)%ztc(l)-ap(p)%zt1(l)))+(fty2*.25*(ap(p)%ztc(l)-ap(p)%zt2(l))) + &
!       !v (fty3*.25*(ap(p)%ztc(l)-ap(p)%zt3(l)))+(fty4*.25*(ap(p)%ztc(l)-ap(p)%zt4(l))) + &
!       !v (ftz1*.25*(ap(p)%yt1(l)-ap(p)%ytc(l)))+(ftz2*.25*(ap(p)%yt2(l)-ap(p)%ytc(l))) + &
!       !v (ftz3*.25*(ap(p)%yt3(l)-ap(p)%ytc(l)))+(ftz4*.25*(ap(p)%yt4(l)-ap(p)%ytc(l)))
!
!       !v ty(l)=(ftx1*.25*(ap(p)%zt1(l)-ap(p)%ztc(l)))+(ftx2*.25*(ap(p)%zt2(l)-ap(p)%ztc(l)))+ &
!       !v (ftx3*.25*(ap(p)%zt3(l)-ap(p)%ztc(l)))+(ftx4*.25*(ap(p)%zt4(l)-ap(p)%ztc(l))) + &
!       !v (ftz1*.25*(ap(p)%xtc(l)-ap(p)%xt1(l)))+(ftz2*.25*(ap(p)%xtc(l)-ap(p)%xt2(l))) + &
!       !v (ftz3*.25*(ap(p)%xtc(l)-ap(p)%xt3(l)))+(ftz4*.25*(ap(p)%xtc(l)-ap(p)%xt4(l)))
!
!
!       !v tz(l)=(ftx1*.25*(ap(p)%ytc(l)-ap(p)%yt1(l)))+(ftx2*.25*(ap(p)%ytc(l)-ap(p)%yt2(l))) + &
!       !v (ftx3*.25*(ap(p)%ytc(l)-ap(p)%yt3(l)))+(ftx4*.25*(ap(p)%ytc(l)-ap(p)%yt4(l))) + &
!       !v (fty1*.25*(ap(p)%xt1(l)-ap(p)%xtc(l)))+(fty2*.25*(ap(p)%xt2(l)-ap(p)%xtc(l))) + &
!       !v (fty3*.25*(ap(p)%xt3(l)-ap(p)%xtc(l)))+(fty4*.25*(ap(p)%xt4(l)-ap(p)%xtc(l)))
!
!        
!        ap(p)%mxl(l)=0. !tx(l)
!        ap(p)%myl(l)=0. !ty(l)
!        ap(p)%mzl(l)=0. !tz(l)
!
!      ! endif  !if l==1...
!
!    enddo    !do l=1,nl    !v I guess, at the end the tension is computed here
!
!
!
!    ! !v this is d/ds (D(M)) as per Eq. (3.36), not needed in my case, need to find distributed force in other way,
!    ! do l=1,nl                               
!    !   if (l==1) then
!    !     txs(l)=(tx(l+1)-tx(l))/ds       
!    !     tys(l)=(ty(l+1)-ty(l))/ds
!    !     tzs(l)=(tz(l+1)-tz(l))/ds
!    !   elseif (l==nl) then
!    !     txs(l)=(tx(l)-tx(l-1))/ds
!    !     tys(l)=(ty(l)-ty(l-1))/ds
!    !     tzs(l)=(tz(l)-tz(l-1))/ds
!    !   else
!    !     txs(l)=(tx(l+1)-tx(l-1))/(2.*ds)
!    !     tys(l)=(ty(l+1)-ty(l-1))/(2.*ds)
!    !     tzs(l)=(tz(l+1)-tz(l-1))/(2.*ds)
!    !   endif
!    ! enddo
!!v this is d/ds (D(M)) as per Eq. (3.36), not needed in my case, need to find distributed force in other way,
!    ! do l=1,nl                               
!    !   if (l==1) then
!    !     txs(l)=(tx(l+1)-tx(l))/ds       
!    !     tys(l)=(ty(l+1)-ty(l))/ds
!    !     tzs(l)=(tz(l+1)-tz(l))/ds
!    !   elseif (l==nl) then
!    !     txs(l)=(tx(l)-tx(l-1))/ds
!    !     tys(l)=(ty(l)-ty(l-1))/ds
!    !     tzs(l)=(tz(l)-tz(l-1))/ds
!    !   else
!    !     txs(l)=(tx(l+1)-tx(l-1))/(2.*ds)
!    !     tys(l)=(ty(l+1)-ty(l-1))/(2.*ds)
!    !     tzs(l)=(tz(l+1)-tz(l-1))/(2.*ds)
!    !   endif
!    ! enddo
!
!  
!    ! do l=1,nl   !v computed fyl and fzl from fry and frz
!
!    !   if (ap(p)%dyds(l)==0. .and. ap(p)%dzds(l)==0.) then
!    !     fry=0.
!    !     frz=0.
!    !   elseif (ap(p)%dyds(l)/=0. .and. ap(p)%dzds(l)==0.) then
!    !     fry= (txs(l)*ap(p)%dzds(l))/((ap(p)%dyds(l)**2.)+(ap(p)%dzds(l)**2.))
!    !     frz= -txs(l)/ap(p)%dyds(l)
!    !   else
!    !     fry= (txs(l)*ap(p)%dzds(l))/((ap(p)%dyds(l)**2.)+(ap(p)%dzds(l)**2.))
!    !     frz=-(fry*ap(p)%dyds(l))/ap(p)%dzds(l)
!    !   endif
!
!    !   ap(p)%fyl(l)=ap(p)%fyl(l)-fry
!    !   ap(p)%fzl(l)=ap(p)%fzl(l)-frz
!    ! enddo
!
!
!    ! do l=1,nl  !v computed fxl and fyl from frx and fry
!
!    !   if (ap(p)%dxds(l)==0..and.ap(p)%dyds(l)==0.) then
!    !     frx=0.
!    !     fry=0.
!    !   elseif (ap(p)%dxds(l)/=0..and.ap(p)%dyds(l)==0.) then
!    !     frx= (tzs(l)*ap(p)%dyds(l))/((ap(p)%dxds(l)**2.)+(ap(p)%dyds(l)**2.))
!    !     fry=-tzs(l)/ap(p)%dxds(l)
!    !   else
!    !     frx= (tzs(l)*ap(p)%dyds(l))/((ap(p)%dxds(l)**2.)+(ap(p)%dyds(l)**2.))
!    !     fry=-(frx*ap(p)%dxds(l))/ap(p)%dyds(l)
!    !   endif
!      
!    !   ap(p)%fxl(l)=ap(p)%fxl(l)-frx
!    !   ap(p)%fyl(l)=ap(p)%fyl(l)-fry
!
!    ! enddo
!
!
!
!    ! do l=1,nl
!
!    !   if (ap(p)%dzds(l)==0..and.ap(p)%dxds(l)==0.) then
!    !     frz=0.
!    !     frx=0.
!    !   elseif (ap(p)%dzds(l)/=0..and.ap(p)%dxds(l)==0.) then
!    !     frz= (tys(l)*ap(p)%dxds(l))/((ap(p)%dzds(l)**2.)+(ap(p)%dxds(l)**2.))
!    !     frx=-tys(l)/ap(p)%dzds(l)
!    !   else
!    !     frz= (tys(l)*ap(p)%dxds(l))/((ap(p)%dzds(l)**2.)+(ap(p)%dxds(l)**2.))
!    !     frx=-(frz*ap(p)%dzds(l))/ap(p)%dxds(l)
!    !   endif
!      
!    !   ap(p)%fzl(l)=ap(p)%fzl(l)-frz
!    !   ap(p)%fxl(l)=ap(p)%fxl(l)-frx
!    ! enddo
!  
!    ! do l=1,nl   !v computed fyl and fzl from fry and frz
!
!    !   if (ap(p)%dyds(l)==0. .and. ap(p)%dzds(l)==0.) then
!    !     fry=0.
!    !     frz=0.
!    !   elseif (ap(p)%dyds(l)/=0. .and. ap(p)%dzds(l)==0.) then
!    !     fry= (txs(l)*ap(p)%dzds(l))/((ap(p)%dyds(l)**2.)+(ap(p)%dzds(l)**2.))
!    !     frz= -txs(l)/ap(p)%dyds(l)
!    !   else
!    !     fry= (txs(l)*ap(p)%dzds(l))/((ap(p)%dyds(l)**2.)+(ap(p)%dzds(l)**2.))
!    !     frz=-(fry*ap(p)%dyds(l))/ap(p)%dzds(l)
!    !   endif
!
!    !   ap(p)%fyl(l)=ap(p)%fyl(l)-fry
!    !   ap(p)%fzl(l)=ap(p)%fzl(l)-frz
!    ! enddo
!
!
!    ! do l=1,nl  !v computed fxl and fyl from frx and fry
!
!    !   if (ap(p)%dxds(l)==0..and.ap(p)%dyds(l)==0.) then
!    !     frx=0.
!    !     fry=0.
!    !   elseif (ap(p)%dxds(l)/=0..and.ap(p)%dyds(l)==0.) then
!    !     frx= (tzs(l)*ap(p)%dyds(l))/((ap(p)%dxds(l)**2.)+(ap(p)%dyds(l)**2.))
!    !     fry=-tzs(l)/ap(p)%dxds(l)
!    !   else
!    !     frx= (tzs(l)*ap(p)%dyds(l))/((ap(p)%dxds(l)**2.)+(ap(p)%dyds(l)**2.))
!    !     fry=-(frx*ap(p)%dxds(l))/ap(p)%dyds(l)
!    !   endif
!      
!    !   ap(p)%fxl(l)=ap(p)%fxl(l)-frx
!    !   ap(p)%fyl(l)=ap(p)%fyl(l)-fry
!
!    ! enddo
!
!
!
!    ! do l=1,nl
!
!    !   if (ap(p)%dzds(l)==0..and.ap(p)%dxds(l)==0.) then
!    !     frz=0.
!    !     frx=0.
!    !   elseif (ap(p)%dzds(l)/=0..and.ap(p)%dxds(l)==0.) then
!    !     frz= (tys(l)*ap(p)%dxds(l))/((ap(p)%dzds(l)**2.)+(ap(p)%dxds(l)**2.))
!    !     frx=-tys(l)/ap(p)%dzds(l)
!    !   else
!    !     frz= (tys(l)*ap(p)%dxds(l))/((ap(p)%dzds(l)**2.)+(ap(p)%dxds(l)**2.))
!    !     frx=-(frz*ap(p)%dzds(l))/ap(p)%dxds(l)
!    !   endif
!      
!    !   ap(p)%fzl(l)=ap(p)%fzl(l)-frz
!    !   ap(p)%fxl(l)=ap(p)%fxl(l)-frx
!    ! enddo
!
!
!
!  endif
!enddo




!$omp end parallel
!
return
end subroutine ext_force
end module mod_ext_force
