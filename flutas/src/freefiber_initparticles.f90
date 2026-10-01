module initparticles_fibm
use decomp_2d
use mod_param_fibm,  only : dims, ndims, lx, ly, lz, ds, &
                        dx, dy, dz, nderiv, &
                        L_fibr, bc0_dof, bcL_dof, ndof, &
                        np, pi, fboundary, &
                        ifattached, ifcurved, nrow, &
                        ncolumn, deg_ele, add_ele, noCP,  &
                        uKnoti, nel, nl, offset, &
                        nno, p_ordr, p1_fordr, ngps, &
                        ngpb, rhof, radius, datadir, npi, npj, &
												npk, Eflag, rhoA, rhoI_yy, rhoI_zz, diametr
use mod_common_fibm, only : all_dof, ir_dof, uKnotp, coefsp, &
                        coefsp_temp, coefsp_np, uKnoth, coefsh, &
												coefsh_temp, coefsh_np, uKnotad, coefsad, &
												coefsad_temp, vec_dummy, kntins, &
												uKnot, coefs, coefs_temp, coefs_np, &
                        ap, npmax, pmax, npmstr, Xi, l, & !nb
                        CONS, unqKnot, elRangeU, xn, &
                        xigs, xigb, wgs, wgb, J_s, J_b, &
                        Arf, I_yy, I_zz, J_xx, NdN_array, &
                        NdN_arrayb, mInerCP, theta_Spagetti, &
                        phi_Spagetti, mInerCP2, strain_k0

use mod_common_mpi,  only : ierr, comm_cart, myid, MPI_INTEGER, MPI_REAL8, &
												MPI_SUM, status
use mod_initialSetup, only: setupInitial
use mod_linspace
use bspline



implicit none
private
public initparticles
contains
subroutine initparticles
implicit none
!integer,parameter extrapoints = nint(2*(radius+0.1)*dxi)! for visualisation purposes
real, dimension(np) :: xcglob,ycglob,zcglob,thetacglob,phicglob
integer :: i,j,k,p,pp,rk
integer :: proccoords(1:ndims),procrank
integer, dimension(2) :: sbuf,rbuf
real :: leftbound,rightbound,frontbound,backbound
real :: dist,distx,disty,distz,distzw,angle
real :: xp,yp!v (not there in gen ver),distzw2,xf,zf
real :: ax
real :: ay
real :: rn
integer :: counter,crys
integer :: idp
character(len=5) rankpr
character(len=7) :: tempstr
character(len=11) :: tempstr2
integer :: count_mstr,count_slve,count_mstr_all,count_slve_loc
real, dimension(NL) :: xfp_temp,yfp_temp,zfp_temp,qn1_temp,qn2_temp,qn3_temp,qn4_temp,theta
integer:: ifil,jfil,kfil
!v-----------------------------------------------------------------------------------------------
REAL(8)  :: coefsi(4,2)
REAL(8)  :: coefsi_np(4,2,np), coefsi_temp(4,2)														 
!REAL(8)     :: coefsi(4,2) = RESHAPE( (/ 0.0, 0.0, 0.0, & 
!                1.0, L_fibr, 0.0, 0.0, 1.0/), (/ 4, 2 /) )!redefined below
coefs = RESHAPE( (/ 0.0, 0.0, 0.0, 1.0, L_fibr, 0.0, 0.0, 1.0/), (/ 4, 2 /) )

do p=1,np
  do i=1,4
    do j=1,2
      coefsi_np(i,j,p) = coefsi(i,j)
      !print*,coefsi_np(i,j,np)
    enddo
  enddo
enddo


!-- compute reduced DOFs
if (fboundary=="h" .or. fboundary=="c") then
  counter=1
  do i = (bc0_dof(ubound(bc0_dof,1))+1), ndof !bc0_dof(ubound(bc0_dof,1)), bcL_dof(1)
    ir_dof(counter) = all_dof(i)
    counter = counter+1
  end do  
elseif (fboundary=="f") then
  do i = 1, ndof !bc0_dof(ubound(bc0_dof,1)), bcL_dof(1)
    ir_dof(i) = all_dof(i)
  end do
endif
! counter=0
!print*, "----------------------------dof",ir_dof
!stop

!
! position of spheres: global initialisation by root (myid = 0)
!
crys=0
theta=-0.0*pi


do p=1,np
  call random_number(rn)
  ! theta_Spagetti(p) = rn*pi!(pi/2.)-((pi/18.)*(-1.**int(10.*rn)))!rn*pi
  theta_Spagetti(p) = (pi/2.)-((pi/18.)*(-1.**int(10.*rn)))!rn*pi
  call random_number(rn)
  !phi_Spagetti(p) =rn*pi
  phi_Spagetti(p) =rn*2.*pi
enddo



if (myid .eq. 0) then
  !!!call random_seed( put = (/16112006/))
  write(6,*) 'part as crystal 1, random 0: ',crys
  open(23,file=datadir//"position_filaments.txt")
  if(crys.eq.1) then
    !    p=0
    !    ! pseudo-crystal
    !    do k=1,4
    !      do j=1,16
    !        do i=1,8
    !          p=p+1
    !          thetacglob(p) = 0.
    !          phicglob(p)   = 0.
    !          call random_number(rn)
    !          xcglob(p)     = (lx/4.)*((1.*i)-0.5) + 0.25*(rn-0.5)*radius
    !          call random_number(rn)
    !          ycglob(p)     = (ly/8.)*((1.*j)-0.5) + 0.25*(rn-0.5)*radius
    !          call random_number(rn)
    !          zcglob(p)     = (lz/2. )*((1.*k)-0.5) + (rn-0.5)*radius
    !          write(6,*) 'Location of sphere #',p
    !          write(6,*) 'x,y,z = ',xcglob(p),ycglob(p),zcglob(p)
    !          write(23,'(I3,5E32.16)') p,thetacglob(p),phicglob(p), &
    !                                  xcglob(p),ycglob(p),zcglob(p)
    !        enddo
    !      enddo
    !    enddo
  else

    ! pseudo-random
    counter = 0
    do p=1,NP
      if(NP<=5) then

        if (ifattached==0) then

          jfil=int((p-1)/(npi*npk))+1
          ifil=int(((p-((jfil-1)*npi*npk))-1)/real(npk))+1
          kfil=(p-((jfil-1)*npi*npk))-((ifil-1)*npk)

          !            xcglob(p)     = ifil*(lx/(npi+1))
          !            ycglob(p)     = jfil*(ly/(npj+1))
          !            zcglob(p)     = kfil*2.*(lz/((2*npk)+2))+2.

          !            xcglob(p)     = lx/2.
          !            ycglob(p)     = ly/2.
          !            zcglob(p)     = lz/2.

          !if (p==1) then
            xcglob(p)     = lx/2.
            ycglob(p)     = 0.5*ly !+ 0.5*L_fibr*cos(-0.1*pi) !ly/2.
            zcglob(p)     = lz/2. !- 0.5*L_fibr*sin(-0.1*pi) !lz/2.
            
						!----
						!ycglob(p)     = 2.0 - 0.5*diametr !ly/2.
            !zcglob(p)     = abs(0.5*L_fibr - lz*(p-1)) ! 2 filament case
						!----
            !ycglob(p)     = ly/2. 
            !ycglob(p)     = 4.0 + 0.5*(p) 
            !ycglob(p)     = 3.5 + 0.5*(p) 
            !zcglob(p)     = 0.5*L_fibr !abs(0.5*L_fibr - lz*(p-1)) ! 2 filament case
            !zcglob(p)     = abs(0.5*L_fibr - lz)

						!print*,"------xcglob, ycglob, zcglob", xcglob(p), ycglob(p), zcglob(p)
            
            ! ! based on ceter point start and end point of spline curve is created in 3D
            ! coefsi(1,1:2) = xcglob(p)
            ! coefsi(2,1:2) = (/ ycglob(p)-0.5*cos(0.1*pi), ycglob(p)+0.5*cos(0.1*pi) /)
            ! coefsi(3,1:2) = (/ zcglob(p)-0.5*sin(0.1*pi), zcglob(p)+0.5*sin(0.1*pi) /)

          !else if (p==2) then
            !xcglob(p)     = lx/2.
            !ycglob(p)     = ly/2.
            !zcglob(p)     = 5.+(1.1*aspectratio)
          !endif

        elseif (ifattached==1) then

          xcglob(p)     = ((lx/nrow)*(p-(floor(real(p-1)/nrow)*nrow)-1))+(lx/(2.*nrow))
          ycglob(p)     = ((ly/ncolumn)*floor(real(p-1)/nrow))+(ly/(2.*ncolumn))-.5
          zcglob(p)     = 4.*dz

        endif

      else      
        111  continue
        call random_number(rn)
        xcglob(p)     = lx*rn
        call random_number(rn)
        ycglob(p)     = ly*rn
        call random_number(rn)
        zcglob(p)     = lz-(lz*rn)
        !v below section not there in genVer
				! do l=1,nl
        !   xf = xcglob(p)+(l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p))
        !   zf = zcglob(p) + (l-int((NL+1)/2.))*ds*cos(theta_Spagetti(p))
        !   distzw = min(abs(lz-zf),abs(zf))
        !   distzw2 = min(abs(lx-xf),abs(xf))
        !   distz = abs(zcglob(p)- zcglob(pp))
        !   if (myid .eq. 0) write(*,*) l,"done 1.1"
        !   if(distzw.lt.0.05*radius.or.distzw2.lt.0.05*radius) goto 111
        !   if (myid .eq. 0) write(*,*) l,"done 1.2"
        ! enddo
        !do pp=1,p
        !  do j=-1,1
        !    disty = abs(ycglob(p)- (ycglob(pp)+j*ly))
        !    do i=-1,1
        !      distx = abs(xcglob(p)- (xcglob(pp)+i*lx))
        !      if(distx.gt.dx.or. &
        !        disty.gt.dy.or. &
        !        distz.gt.dz.or.p.eq.pp) then
        !        ! good particle
        !      else
        !        dist = distx**2+disty**2.+distz**2.
        !        if((dist.lt.(4.2*radius**2.))) then 
        !          !write(*,*)'RANDOM DEVIATION'
        !          !write(*,*)dist,distw
        !          counter=counter+1
        !          goto 111
        !        endif
        !      endif
        !    enddo
        !  enddo
        !222 continue
        !enddo
        !endif
				!      write(6,*) 'Location of filament #',p
      	!      write(6,*) 'x,y,z = ',xcglob(p),ycglob(p),zcglob(p)
      	!      write(23,'(I3,3E32.16)') p,xcglob(p),ycglob(p),zcglob(p) 
			do pp=1,p
          distzw = min(abs(lz-zcglob(p)),abs(zcglob(p)))
          distz = abs(zcglob(p)- zcglob(pp))
          if(distzw.lt.1.05*radius) goto 111
          do j=-1,1
            disty = abs(ycglob(p)- (ycglob(pp)+j*ly))
            do i=-1,1
              distx = abs(xcglob(p)- (xcglob(pp)+i*lx))
              if(distx.gt.dx.or. &
                 disty.gt.dy.or. &
                 distz.gt.dz.or.p.eq.pp) then
                ! good particle
              else
                dist = distx**2+disty**2.+distz**2.
                if((dist.lt.(4.2*radius**2.))) then
                  !write(*,*)'RANDOM DEVIATION'
                  !write(*,*)dist,distw
                  counter=counter+1
                  goto 111
                endif
              endif
            enddo
          enddo
        222 continue
        enddo
      endif
      write(6,*) 'Location of filament #',p
      write(6,*) 'x,y,z = ',xcglob(p),ycglob(p),zcglob(p)
      write(23,'(I3,3E32.16)') p,xcglob(p),ycglob(p),zcglob(p)
                               
    enddo
    write(*,*)'RANDOM DEVIATIONS: ',counter
    p=p-1
    !    p = np
  endif
  close(23)

  if ( (p.ne.np) ) then
    print*,counter,np
    write(6,*) 'Fatal error in initialisation of particle positions!'
    write(6,*) 'Program aborted...'
    call mpi_finalize(ierr) !error
    stop
  endif

   
  do rk=1,Nproc-1
    call MPI_SSEND(xcglob,np,MPI_REAL8,rk,rk+0*(Nproc-1),comm_cart,ierr) !error
    call MPI_SSEND(ycglob,np,MPI_REAL8,rk,rk+1*(Nproc-1),comm_cart,ierr)
    call MPI_SSEND(zcglob,np,MPI_REAL8,rk,rk+2*(Nproc-1),comm_cart,ierr)
  enddo

else ! if myid is not 0:
  call MPI_RECV(xcglob,np,MPI_REAL8,0,myid+0*(Nproc-1),comm_cart,status,ierr) !error
  call MPI_RECV(ycglob,np,MPI_REAL8,0,myid+1*(Nproc-1),comm_cart,status,ierr)
  call MPI_RECV(zcglob,np,MPI_REAL8,0,myid+2*(Nproc-1),comm_cart,status,ierr)
endif


!
! Refine the coarse version of a (straight/curved) filament
!
do p=1,NP
   !based on ceter point start and end point of spline curve is created in 3D
   coefsi_np(1,1:2,p) = xcglob(p) !(/ xcglob(p)-0.5, xcglob(p)+0.5 /)
   !coefsi_np(2,1:2,p) = (/ ycglob(p)-0.5*L_fibr*cos(-0.0*pi), ycglob(p)+0.5*L_fibr*cos(-0.0*pi) /)
   !coefsi_np(3,1:2,p) = (/ zcglob(p)-0.5*L_fibr*sin(-0.0*pi), zcglob(p)+0.5*L_fibr*sin(-0.0*pi) /)
	 !coefsi_np(2,1:2,p) = (/ ycglob(p)-0.5*L_fibr*cos(-0.1*pi), ycglob(p)+0.5*L_fibr*cos(-0.1*pi) /)
   !coefsi_np(3,1:2,p) = (/ zcglob(p)-0.5*L_fibr*sin(-0.1*pi), zcglob(p)+0.5*L_fibr*sin(-0.1*pi) /)
	 coefsi_np(2,1:2,p) = (/ ycglob(p)-0.5*L_fibr*cos(-1.0*pi/180), ycglob(p)+0.5*L_fibr*cos(-1.0*pi/180) /)
   coefsi_np(3,1:2,p) = (/ zcglob(p)-0.5*L_fibr*sin(-1.0*pi/180), zcglob(p)+0.5*L_fibr*sin(-1.0*pi/180) /)


  
   !----
   !coefsi_np(1,1:2,p) = xcglob(p) !(/ xcglob(p)-0.5, xcglob(p)+0.5 /)
   !coefsi_np(2,1:2,p) = ycglob(p) !(/ ycglob(p)-0.5, ycglob(p)+0.5 /)
   !coefsi_np(3,1:2,p) = (/ zcglob(p)+0.5*L_fibr*(-1)**p, zcglob(p)-0.5*L_fibr*(-1)**p /)

	 !----
   !coefsi_np(1,1:2,p) = xcglob(p) 
   !coefsi_np(2,1:2,p) = ycglob(p)
   !coefsi_np(3,1:2,p) = (/ zcglob(p)-0.5*L_fibr, zcglob(p)+0.5*L_fibr /)

	 !do i = 1,4
	 ! 	print '(2E10.3)', (coefsi_np(i,j,p),j=1,2)
	 !enddo
	 !print*,"------------ coefsi_np ----------"

	 !print*,"------size coefsp_np", size(coefsi_np,2)


	 if (ifcurved==0) then


	   !---- Order Elevatation
	   if (deg_ele .GT. 0) then
			 
			 if (p .eq. 1) then
	 	 	   ALLOCATE( uKnotp( noCP + p_ordr + 1 + deg_ele*2) )
	 	 	   ALLOCATE( coefsp(4, noCP-1 + deg_ele + 1) ) 
	 	 	   ALLOCATE( coefsp_np(4, noCP-1 + deg_ele + 1, np) )
	 	 	   ALLOCATE( coefsp_temp(size(coefsp_np,1), size(coefsp_np,2)) )
	 	 	 endif


			 !--- initialize
     	 coefsi_temp = 0.0D+00
			 coefsp_temp = 0.0D+00
     	 do i=1,4
     	    do j=1,2!(size(coefsi_np,2))
     	       coefsi_temp(i,j)  = coefsi_np(i,j,p)
     	    enddo
    	    !do j=1,(size(coefsp_np,2))
     	    !   coefsp_temp(i,j) = coefsp_np(i,j,p)
     	    !enddo
					!print'(11E14.2)', (coefsp_np(i,j,1), j = 1,size(coefsp_np,2))
					!print'(11E14.2)', (coefsi_temp(i,j), j = 1,size(coefsi_temp,2))
     	 enddo
			 !print*, "-----end of coefsp_temp print----"
       !print*,"----coefsi_temp", shape(coefsi_temp)
			 !print*,"----porder", p_ordr
       !print*,"----uKnoti", uKnoti
       !print*,"----uKnotp", uKnotp
       !print*,"---- deg_ele", deg_ele 
       !print*,"---- noCP", noCP 
			 !print*,"----coefsp_temp", shape(coefsp_temp)
			 !stop


			!do i=1,4
			! print'(11E14.2)', (coefsi_temp(i,j), j = 1,size(coefsi_temp,2))
			!enddo	

			 !--- Degree Elevation
	     !CALL DegreeElevate(4, size(coefsi_np,2)-1, p_ordr, uknoti, & 
	     !                    coefsi_temp, deg_ele, noCP-1, uknotp, coefsp_temp)

			 CALL DegreeElevate(4, size(coefsi_np,2)-1, 1, uknoti, &
                          coefsi_temp, deg_ele, noCP-1, uknotp, coefsp_temp)

			 !print*, "size of coefsp_np:    ", shape(coefsp_np)
			 !--- store 
			 !print'(11E14.2)', deg_ele 
			 !print*, "p_ordr", p_ordr
			 !print*,"deg_ele", deg_ele
			 !print*, uKnoti 
			 do i=1,4
					!print'(11E14.2)', (coefsp_temp(i,j), j = 1,size(coefsp_temp,2))
    	    do j=1,(size(coefsp_temp,2))
    	       coefsp_np(i,j,p) = coefsp_temp(i,j)
    	    enddo
					!print'(11E14.2)', (coefsp_np(i,j,p), j = 1,size(coefsp_np,2))
    	 enddo
			 !print*, "-----end of coefsp_np after deg ele----"
			 !if (p .eq. 2) stop
										 
													 
			 if (p .eq. 1) then
					p_ordr = p_ordr + deg_ele        
	        noCP   = size(coefsp_np,2)
			 endif		

	   else
      
	     uKnotp = uKnoti
			 if (p .eq. 1) then
			   ALLOCATE( coefsp_np(4, noCP, np) )
			   ALLOCATE( coefsp_temp(size(coefsp_np,1), size(coefsp_np,2)) )
			 endif	 

			 !print*, "coefsp_np w/o degree elevate"
  		 do i=1,4
					do j=1,size(coefsp_np,2)
						 coefsp_np(i,j,p) = coefsi_np(i,j,p)
					enddo
					!print'(11E14.2)', (coefsp_np(i,j,p), j = 1,size(coefsp_np,2))
			 enddo	
			 !coefsp_np(i,j,p) = coefsi_np(i,j,p)

	   endif
	 
	   ! print*,"-------------------------coefsi 1 ",coefsi(1,:)
	   ! print*,"-------------------------coefsi 2 ",coefsi(2,:)
	   ! print*,"-------------------------coefsi 3 ",coefsi(3,:)
	   ! print*,""
	   
	   ! print*,"---------------------------coefs1 before refineMesh", coefsp(1,:)
	   ! print*,"---------------------------coefs2 before refineMesh", coefsp(2,:)
	   ! print*,"---------------------------coefs3 before refineMesh", coefsp(3,:)
		 !do i=1,4
     !   print'(11E14.2)', (coefsp_np(i,j,p), j=1,size(coefsp_np,2))
     !enddo
     !print*, "-----end of coefsp_np after deg ele----"


	   !------  Refine Mesh
		 if (p .eq. 1)  then 
				ALLOCATE(vec_dummy(nel+1), kntins(nel-1))
			  !vec_dummy = linspace( REAL(coefsp(2,1),8), REAL(coefsp(2,),8), nel+1)
				!vec_dummy = linspace(0.0_dp, L_fibr*1.0_dp, nel+1)
				vec_dummy = linspace(0.0_dp, 1.0_dp, nel+1)
				kntins    = vec_dummy(2 : nel)    
				!print*, "kntins", kntins
				ALLOCATE( uKnoth(size(uKnotp) + size(kntins) ),  & 
					 coefsh_temp(size(coefsp_np,1), size(coefsp_np,2)+size(kntins) ), &
	         coefsh_np(size(coefsp_np,1), size(coefsp_np,2)+size(kntins), np))
			  !ALLOCATE( coefs_temp(size(coefs_np,1), size(coefs_np,2) ) )
		 endif			 
		
		 !print*,"---shape coefsp_1", 	size(coefsp,1)
		 !print*,"---shape coefsp_temp", shape(coefsp_temp)
		 !print*,"---shape coefs_temp", shape(coefs_temp)
		 !print*, "coefsp_temp --> input to refine knot vector"

		 !--- initialize
		 coefsp_temp = 0.0D+00
		 coefsh_temp = 0.0D+00
		 do i=1,4
        do j=1,(size(coefsp_np,2))
           coefsp_temp(i,j) = coefsp_np(i,j,p)
        enddo
			  !print'(11E14.2)', (coefsp_temp(i,j), j = 1,size(coefsp_temp,2))
        !do j=1,(size(coefsp_np,2))
        !   coefsh_temp(i,j)  = coefsp_np(i,j,p)
        !enddo
			  !print'(11E14.2)', (coefs_temp(i,j), j=1,size(coefs_temp,2))
     enddo
		 !stop
	
		 !print*,'-----end of coefsp_temp-----'
		 !print*,'-----end of coefs-----'

		 !print*,"---p:   ", p
		 !do i=1,4
     !   print '(2E10.2)', (coefs_temp(i,j), j=1,p1_fordr)
     !enddo
     !print*, "-----end pring coefsp_temp"
		 !print*, "p_ordr", p_ordr

	   CALL RefineKnotVector(4, size(coefsp_np,2)-1, p_ordr, &   
         uKnotp, coefsp_temp, size(kntins,1)-1, kntins, uKnoth, coefsh_temp)
		 !CALL RefineKnotVector(4, size(coefsp,2)-1, p_ordr, & 
	   !    uKnotp, coefsp_np(:,:,np), size(kntins,1)-1, kntins, uKnot, coefs_np(:,:,np))

		 
		 !print*,"---p:   ", p
     !do i=1,4
		 ! 	print '(2E10.2)', (coefsp_temp(i,j), j=1,nno)
		 !enddo
		 !print*, "-----end pring coefsp_temp"
		 !print*, "coefsh_temp after mesh refinement"
		 !--- store 
     do i=1,4
        do j=1,(size(coefsh_np,2))
           coefsh_np(i,j,p) = coefsh_temp(i,j)
        enddo
			  !print '(11E10.2)', (coefsh_temp(i,j), j=1,size(coefsh_temp,2))
				!print '(11E10.2)', (coefsh_np(i,j,p), j=1,size(coefsh_np,2))
     enddo
		 !print*, "----- end coefs_np ---- after mesh_refine"
		 !stop


		 if (p .eq. 1) then
			 noCP = size(coefsh_np,2)
		 endif
	   !print*,"---------------------------coefs1 after refineMesh", coefs(1,:)
	   !print*,"---------------------------coefs2 after refineMesh", coefs(2,:)
	   !print*,"---------------------------coefs3 after refineMesh", coefs(3,:)
	 	 !DEALLOCATE(uKnotp, coefsp)
	 	 !DEALLOCATE(uKnotp, coefsp, uKnot)



	   !------ Additional order elevation 
		 if (add_ele .gt. 0) then
			 
			 if (p .eq. 1)  then 
        unqKnot = uKnoth(p_ordr+1 : ubound(uKnoth,1)-(1+deg_ele))
				ALLOCATE( uKnot(size(uKnoth) + add_ele*size(unqKnot) ),  & 
					 coefs_temp(size(coefsh_np,1), size(coefsh_np,2)+add_ele*nel ), &
	         coefs_np(size(coefsh_np,1), size(coefsh_np,2)+add_ele*nel, np))
					 !coefs_temp(size(coefsh_np,1), size(coefsh_np,2)+add_ele*size(unqKnot)-1 ), &
	         !coefs_np(size(coefsh_np,1), size(coefsh_np,2)+add_ele*size(unqKnot)-1, np))
			 endif


			 coefsh_temp = 0.0D+00
		 	 do i=1,4
         do j=1,(size(coefsh_np,2))
            coefsh_temp(i,j) = coefsh_np(i,j,p)
         enddo
			 enddo	

			 CALL DegreeElevate(4, size(coefsh_np,2)-1, 1+deg_ele, uknoth, &
                        coefsh_temp, add_ele, noCP-1, uknot, coefs_temp)
			  								
			 do i=1,4
     	    do j=1,(size(coefs_np,2))
     	       coefs_np(i,j,p) = coefs_temp(i,j)
     	    enddo
				  !print '(11E10.2)', (coefs_np(i,j,p), j=1,size(coefs_np,2))
		 	 enddo		
			 !print*,"----------------"

			 if (p .eq. 1) then
					!print*,"----------------p_ordr", p_ordr
			 		!print*,"----------------coefs_np", size(coefs_np,2)
			 		!print*,"----------------p_ordr", p_ordr
	        noCP   = size(coefs_np,2)
			  	!p_ordr = 1 + deg_ele + add_ele 
			 endif		

		else
		   !	
	     uKnot = uKnoth
			 if (p .eq. 1) then
			   ALLOCATE( coefs_np(4, nno, np) )
			 endif	 
  		 do i=1,4
					do j=1,size(coefs_np,2)
						 coefs_np(i,j,p) = coefsh_np(i,j,p)
					enddo
				  !print '(11E10.2)', (coefs_np(i,j,p), j=1,size(coefs_np,2))
			 enddo	
       !
			 !print*, "============= Coming to add_ele == 0 case ============"
			 !stop
			 !print*, uKnot


		endif										
	 


	 elseif (ifcurved==1) then
	 
	   print*, "---------------------------------(ifcurved==1) yet to be done"
	   stop
	 endif
enddo

!print*, "--------uKnot----------", real(uKnot,4)
!print*, "--------nno----------", nno
!print*, "--------p1_fordr---------", p1_fordr




!
! Reading coordinates for curved fiber 
!
!open(unit=199,file='curved_fib_coords.out',status='old',action='read',iostat=ierr)
!do j = 1,nno
!		read(199,*) coefs_np(1:3,j,1)
!	  !print '(11E10.2)', (coefs_np(i,j,1), i=1,3)
!enddo







!
! Determine master and slave processes for each particle.
!
! initialisation
!
ap(1:npmax)%x = 0.
ap(1:npmax)%y = 0.
ap(1:npmax)%z = 0.
ap(1:npmax)%mslv = 0.
forall(i=1:npmax) ap(i)%nb(1:8) = 0
!
count_mstr = 0
count_slve = 0
i = 0
ax = 0.5
ay = 0.5
!
pmax = 0
!
!
do p=1,np
  if (xcglob(p).lt.0..or.xcglob(p).gt.lx .or. &
    ycglob(p).lt.0..or.ycglob(p).gt.ly .or. &
    zcglob(p).lt.radius.or.zcglob(p).gt.lz-radius) then
    if (myid.eq.0) then
      write(6,*) 'Fatal error in initialisation of particle positions - '
      write(6,*) 'particle outside the domain!'
      write(6,*) 'Program aborted...'
    endif
    !call mpi_finalize(error)
    !stop
  endif
  if (xcglob(p).eq.lx) ax = 0.51
  if (xcglob(p).eq.0) ax = 0.49
  if (ycglob(p).eq.ly) ay = 0.51
  if (ycglob(p).eq.0) ay = 0.49
  !  
  proccoords(1) = nint(xcglob(p)*dims(1)/lx - ax)
  proccoords(2) = nint(ycglob(p)*dims(2)/ly - ay)
  leftbound     = (proccoords(1)  )*lx/(1.*dims(1)) ! left  boundary of particle's master
  rightbound    = (proccoords(1)+1)*lx/(1.*dims(1)) ! right boundary of particle's master
  frontbound    = (proccoords(2)  )*ly/(1.*dims(2)) ! front boundary of particle's master
  backbound     = (proccoords(2)+1)*ly/(1.*dims(2)) ! back  boundary of particle's master
  call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr) !error


  if (myid.eq.procrank) then


    i = i + 1
    count_mstr = count_mstr + 1
    ap(i)%x = xcglob(p)
    ap(i)%y = ycglob(p)
    ap(i)%z = zcglob(p)
    ap(i)%mslv = p
    
    
    do l=1,NL


      if (ifcurved==0) then

        !ap(i)%xfp(l) = ap(i)%x + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
        !ap(i)%yfp(l) = ap(i)%y + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)) 
        !ap(i)%zfp(l) = ap(i)%z + (l-int((NL+1)/2.))*ds*cos(theta_Spagetti(p)) 
        
        !ap(i)%xfp(l) = ap(i)%x !+ (l-int((NL+1)/2.))*ds 
        !ap(i)%yfp(l) = ap(i)%y + (l-int((NL+1)/2.))*ds * cos(0.1*pi) 
        !ap(i)%zfp(l) = ap(i)%z + (l-int((NL+1)/2.))*ds * sin(0.1*pi)

        !print*,"------------------->", coefs_np(1,l,p), coefs_np(2,l,p), coefs_np(3,l,p)        
        ap(i)%xfp(l) = coefs_np(1,l,p)!ap(i)%x + (l-int((NL+1)/2.))*ds
        ap(i)%yfp(l) = coefs_np(2,l,p)!ap(i)%y + (l-int((NL+1)/2.))*ds * cos(0.1*pi)
        ap(i)%zfp(l) = coefs_np(3,l,p)!ap(i)%z + (l-int((NL+1)/2.))*ds * sin(0.1*pi)
        ap(i)%qn1(l) = 1.0
        ap(i)%qn2(l) = 0.0
        ap(i)%qn3(l) = 0.0
        ap(i)%qn4(l) = 0.0




      elseif (ifcurved==1) then

        ! if (l==1) then
        ! ap(i)%xfp(l) = ap(i)%x  
        ! ap(i)%yfp(l) = ap(i)%y  
        ! ap(i)%zfp(l) = ap(i)%z 
        ! elseif (l==2) then
        ! ap(i)%xfp(l) = ap(i)%x  
        ! ap(i)%yfp(l) = ap(i)%y+ds  
        ! ap(i)%zfp(l) = ap(i)%z 
        ! else
        ! ap(i)%xfp(l) = ap(i)%xfp(l-1) + ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)*(l-2)/(2.*NL)) 
        ! ap(i)%yfp(l) = ap(i)%yfp(l-1) + ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)*(l-2)/(2.*NL)) 
        ! ap(i)%zfp(l) = ap(i)%zfp(l-1) + ds*cos(theta_Spagetti(p))
        ! endif

      endif

      ap(i)%xfpo(l) = ap(i)%xfp(l)
      ap(i)%yfpo(l) = ap(i)%yfp(l)
      ap(i)%zfpo(l) = ap(i)%zfp(l)    
      ap(i)%qno1(l) = ap(i)%qn1(l)
      ap(i)%qno2(l) = ap(i)%qn2(l)
      ap(i)%qno3(l) = ap(i)%qn3(l)
      ap(i)%qno4(l) = ap(i)%qn4(l)


			!neighbor 1
      if ( (ap(i)%xfp(l)+offset) .gt. rightbound ) then
        if ( ((ap(i)%yfp(l)+offset) .ge. frontbound) .and. ((ap(i)%yfp(l)-offset) .le. backbound) ) then
           ap(i)%nb(1) = 1 !neighbor 1 is slave of particle ap(i)%mslv
        endif
      endif

    !neighbor 2
      if ( (ap(i)%xfp(l)+offset) .gt. rightbound ) then
        if  ( (ap(i)%yfp(l)-offset) .le. frontbound ) then
          ap(i)%nb(2) = 1 !neighbor 2 is slave of particle ap(i)%mslv
        endif
      endif

    !neighbor 3
      if ( (ap(i)%yfp(l)-offset) .lt. frontbound ) then
        if ( ((ap(i)%xfp(l)+offset) .ge. leftbound) .and. ((ap(i)%xfp(l)-offset) .le. rightbound )) then
           ap(i)%nb(3) = 1 !neighbor 3 is slave of particle ap(i)%mslv
        endif
      endif

    !neighbor 4
      if ( (ap(i)%yfp(l)-offset) .lt. frontbound ) then
        if  ( (ap(i)%xfp(l)-offset) .le. leftbound ) then
          ap(i)%nb(4) = 1 !neighbor 4 is slave of particle ap(i)%mslv
        endif
      endif

    !neighbor 5
      if ( (ap(i)%xfp(l)-offset) .lt. leftbound ) then
        if ( ((ap(i)%yfp(l)+offset) .ge. frontbound) .and. ((ap(i)%yfp(l)-offset) .le. backbound) ) then
           ap(i)%nb(5) = 1 !neighbor 5 is slave of particle ap(i)%mslv
        endif
      endif

    !neighbor 6
      if ( (ap(i)%yfp(l)+offset) .gt. backbound ) then
        if  ( (ap(i)%xfp(l)-offset) .le. leftbound ) then
          ap(i)%nb(6) = 1 !neighbor 6 is slave of particle ap(i)%mslv
        endif
      endif

    !neighbor 7
      if ( (ap(i)%yfp(l)+offset) .gt. backbound ) then
        if ( ((ap(i)%xfp(l)+offset) .ge. leftbound) .and. ((ap(i)%xfp(l)-offset) .le. rightbound )) then
           ap(i)%nb(7) = 1 !neighbor 7 is slave of particle ap(i)%mslv
        endif
      endif

    !neighbor 8
      if ( (ap(i)%yfp(l)+offset) .gt. backbound ) then
        if  ( (ap(i)%xfp(l)+offset) .ge. rightbound ) then
          ap(i)%nb(8) = 1 !neighbor 8 is slave of particle ap(i)%mslv
        endif
      endif

    !106 continue  
    enddo

    !  print*, "---------------------------------"
		!  if (ap(p)%mslv .ge. 0) then  
		! 	  do l=1,NL
		! 			print*, ap(p)%xfp(l), ap(p)%yfp(l), ap(p)%zfp(l)
		!     enddo
		!  endif
		!  print*, "---------------------------------"

  else

    
    count_slve_loc = 0
    !neighbor 1 of particle's master
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) + 1
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) 
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid .eq. procrank) then

      do l=1,NL


        if (ifcurved==0) then

          ! xfp_temp(l) = xcglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
          ! yfp_temp(l) = ycglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)) 
          ! zfp_temp(l) = zcglob(p) + (l-int((NL+1)/2.))*ds*cos(theta_Spagetti(p)) 

          xfp_temp(l) = coefs_np(1,l,p)
          yfp_temp(l) = coefs_np(2,l,p)
          zfp_temp(l) = coefs_np(3,l,p)
          !qn1_temp(l) = 1.0
          !qn2_temp(l) = 0.0
          !qn3_temp(l) = 0.0
          !qn4_temp(l) = 0.0  ! didn't do it as not required

        elseif (ifcurved==1) then

          ! if (l==1) then
          !   xfp_temp(l) = xcglob(p)  
          !   yfp_temp(l) = ycglob(p)  
          !   zfp_temp(l) = zcglob(p)
          ! elseif (l==2) then
          !   xfp_temp(l) = xcglob(p)  
          !   yfp_temp(l) = ycglob(p)+ds  
          !   zfp_temp(l) = zcglob(p)      
          ! else
          !   xfp_temp(l) = xfp_temp(l-1) + ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)*(l-2)/(2.*NL)) 
          !   yfp_temp(l) = yfp_temp(l-1) + ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)*(l-2)/(2.*NL)) 
          !   zfp_temp(l) = zfp_temp(l-1) + ds*cos(theta_Spagetti(p))
          ! endif

        endif
        !        xfp_temp(l) = xcglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
        !        yfp_temp(l) = ycglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p))
        !        zfp_temp(l) = zcglob(p) + (l-1)*ds *cos(theta_Spagetti(p))   

        if ( (xfp_temp(l)+offset) .gt. rightbound ) then
          if ( ((yfp_temp(l)+offset) .ge. frontbound) .and. ((yfp_temp(l)-offset) .le. backbound) ) then 
           !if (typeflow=="d".and.rightbound==lx) goto 201
            if(count_slve_loc.eq. 0 ) i = i+1
              count_slve = count_slve + 1
              ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
              ap(i)%nb(5) = 1      !neighbor 5 of myid is particle's master
              count_slve_loc = count_slve_loc + 1
              goto 201

          endif
        endif  
      enddo

    201 continue  
    endif


    !neighbor 2 of particle's master
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) + 1
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) - 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)!error
    if (myid .eq. procrank) then

      do l=1,NL

        if (ifcurved==0) then

          !xfp_temp(l) = xcglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
          !yfp_temp(l) = ycglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)) 
          !zfp_temp(l) = zcglob(p) + (l-int((NL+1)/2.))*ds*cos(theta_Spagetti(p)) 

          !xfp_temp(l) = coefs(1,l)
          !yfp_temp(l) = coefs(2,l)
          !zfp_temp(l) = coefs(3,l)
          xfp_temp(l) = coefs_np(1,l,p)
          yfp_temp(l) = coefs_np(2,l,p)
          zfp_temp(l) = coefs_np(3,l,p)

        elseif (ifcurved==1) then

          ! if (l==1) then
          !   xfp_temp(l) = xcglob(p)  
          !   yfp_temp(l) = ycglob(p)  
          !   zfp_temp(l) = zcglob(p)
          ! elseif (l==2) then
          !   xfp_temp(l) = xcglob(p)  
          !   yfp_temp(l) = ycglob(p)+ds  
          !   zfp_temp(l) = zcglob(p)      
          ! else
          !   xfp_temp(l) = xfp_temp(l-1) + ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)*(l-2)/(2.*NL)) 
          !   yfp_temp(l) = yfp_temp(l-1) + ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)*(l-2)/(2.*NL)) 
          !   zfp_temp(l) = zfp_temp(l-1) + ds*cos(theta_Spagetti(p))
          ! endif

        endif

        !        xfp_temp(l) = xcglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
        !        yfp_temp(l) = ycglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p))
        !        zfp_temp(l) = zcglob(p) + (l-1)*ds *cos(theta_Spagetti(p))  

        if ( (xfp_temp(l)+offset) .gt. rightbound ) then
          if  ( (yfp_temp(l)-offset) .le. frontbound ) then 
            !if (typeflow=="d".and.rightbound==lx) goto 202
            if(count_slve_loc.eq. 0 ) i = i+1
            count_slve = count_slve + 1
            ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%nb(6) = 1      !neighbor 6 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 202

          endif
        endif  

      enddo

    202 continue 
    endif



    !neighbor 3 of particle's master
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax )
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) - 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid .eq. procrank) then

      do l=1,NL

        if (ifcurved==0) then

          !xfp_temp(l) = xcglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
          !yfp_temp(l) = ycglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)) 
          !zfp_temp(l) = zcglob(p) + (l-int((NL+1)/2.))*ds*cos(theta_Spagetti(p)) 

          !xfp_temp(l) = coefs(1,l)
          !yfp_temp(l) = coefs(2,l)
          !zfp_temp(l) = coefs(3,l)
          xfp_temp(l) = coefs_np(1,l,p)
          yfp_temp(l) = coefs_np(2,l,p)
          zfp_temp(l) = coefs_np(3,l,p)

        elseif (ifcurved==1) then

          ! if (l==1) then
          !   xfp_temp(l) = xcglob(p)  
          !   yfp_temp(l) = ycglob(p)  
          !   zfp_temp(l) = zcglob(p)
          ! elseif (l==2) then
          !   xfp_temp(l) = xcglob(p)  
          !   yfp_temp(l) = ycglob(p)+ds  
          !   zfp_temp(l) = zcglob(p)      
          ! else
          !   xfp_temp(l) = xfp_temp(l-1) + ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)*(l-2)/(2.*NL)) 
          !   yfp_temp(l) = yfp_temp(l-1) + ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)*(l-2)/(2.*NL)) 
          !   zfp_temp(l) = zfp_temp(l-1) + ds*cos(theta_Spagetti(p))
          ! endif

        endif

        !        xfp_temp(l) = xcglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
        !        yfp_temp(l) = ycglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p))
        !        zfp_temp(l) = zcglob(p) + (l-1)*ds *cos(theta_Spagetti(p)) 

        if ( (yfp_temp(l)-offset) .lt. frontbound ) then
          if ( ((xfp_temp(l)+offset) .ge. leftbound) .and. ((xfp_temp(l)-offset) .le. rightbound) ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
              count_slve = count_slve + 1
              ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
              ap(i)%nb(7) = 1      !neighbor 7 of myid is particle's master
              count_slve_loc = count_slve_loc + 1
              goto 203

          endif
        endif  
      enddo

    203 continue 
    endif


    !neighbor 4 of particle's master
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) - 1
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) - 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid .eq. procrank) then

      do l=1,NL

        if (ifcurved==0) then

          !xfp_temp(l) = xcglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
          !yfp_temp(l) = ycglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)) 
          !zfp_temp(l) = zcglob(p) + (l-int((NL+1)/2.))*ds*cos(theta_Spagetti(p)) 

          !xfp_temp(l) = coefs(1,l)
          !yfp_temp(l) = coefs(2,l)
          !zfp_temp(l) = coefs(3,l)
          xfp_temp(l) = coefs_np(1,l,p)
          yfp_temp(l) = coefs_np(2,l,p)
          zfp_temp(l) = coefs_np(3,l,p)


        elseif (ifcurved==1) then


          ! if (l==1) then
          !   xfp_temp(l) = xcglob(p)  
          !   yfp_temp(l) = ycglob(p)  
          !   zfp_temp(l) = zcglob(p)
          ! elseif (l==2) then
          !   xfp_temp(l) = xcglob(p)  
          !   yfp_temp(l) = ycglob(p)+ds  
          !   zfp_temp(l) = zcglob(p)      
          ! else
          !   xfp_temp(l) = xfp_temp(l-1) + ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)*(l-2)/(2.*NL)) 
          !   yfp_temp(l) = yfp_temp(l-1) + ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)*(l-2)/(2.*NL)) 
          !   zfp_temp(l) = zfp_temp(l-1) + ds*cos(theta_Spagetti(p))
          ! endif

          endif

          !        xfp_temp(l) = xcglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
          !        yfp_temp(l) = ycglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p))
          !        zfp_temp(l) = zcglob(p) + (l-1)*ds *cos(theta_Spagetti(p))   

          if ( (yfp_temp(l)-offset) .lt. frontbound ) then
            if ( (xfp_temp(l)-offset) .le. leftbound ) then 
              !if (typeflow=="d".and.leftbound==0.) goto 204
                if(count_slve_loc.eq. 0 ) i = i+1
                count_slve = count_slve + 1
                ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
                ap(i)%nb(8) = 1      !neighbor 8 of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                goto 204
                
            endif
          endif  
      enddo

    204 continue 
    endif


    ! neighbor 5 of particle's master
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) - 1
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay )
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid .eq. procrank) then

      do l=1,NL

        if (ifcurved==0) then

          !xfp_temp(l) = xcglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
          !yfp_temp(l) = ycglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)) 
          !zfp_temp(l) = zcglob(p) + (l-int((NL+1)/2.))*ds*cos(theta_Spagetti(p)) 

          !xfp_temp(l) = coefs(1,l)
          !yfp_temp(l) = coefs(2,l)
          !zfp_temp(l) = coefs(3,l)
          xfp_temp(l) = coefs_np(1,l,p)
          yfp_temp(l) = coefs_np(2,l,p)
          zfp_temp(l) = coefs_np(3,l,p)


        elseif (ifcurved==1) then

          ! if (l==1) then
          !   xfp_temp(l) = xcglob(p)  
          !   yfp_temp(l) = ycglob(p)  
          !   zfp_temp(l) = zcglob(p)
          ! elseif (l==2) then
          !   xfp_temp(l) = xcglob(p)  
          !   yfp_temp(l) = ycglob(p)+ds  
          !   zfp_temp(l) = zcglob(p)      
          ! else
          !   xfp_temp(l) = xfp_temp(l-1) + ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)*(l-2)/(2.*NL)) 
          !   yfp_temp(l) = yfp_temp(l-1) + ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)*(l-2)/(2.*NL)) 
          !   zfp_temp(l) = zfp_temp(l-1) + ds*cos(theta_Spagetti(p))
          ! endif


        endif

        !        xfp_temp(l) = xcglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
        !        yfp_temp(l) = ycglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p))
        !        zfp_temp(l) = zcglob(p) + (l-1)*ds *cos(theta_Spagetti(p))  

        if ( (xfp_temp(l)-offset) .lt. leftbound ) then
          if ( ((yfp_temp(l)+offset) .ge. frontbound) .and. ((yfp_temp(l)-offset) .le. backbound ) ) then 
            !if (typeflow=="d".and.leftbound==0.) goto 205
            if(count_slve_loc.eq. 0 )   i = i+1
            count_slve = count_slve + 1
            ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%nb(1) = 1      !neighbor 1 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 205

          endif
        endif  
      enddo

    205 continue 
    endif


    !neighbor 6 of particle's master
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) - 1
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) + 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid .eq. procrank) then

      do l=1,NL

        if (ifcurved==0) then

          !xfp_temp(l) = xcglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
          !yfp_temp(l) = ycglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)) 
          !zfp_temp(l) = zcglob(p) + (l-int((NL+1)/2.))*ds*cos(theta_Spagetti(p)) 

          !xfp_temp(l) = coefs(1,l)
          !yfp_temp(l) = coefs(2,l)
          !zfp_temp(l) = coefs(3,l)
          xfp_temp(l) = coefs_np(1,l,p)
          yfp_temp(l) = coefs_np(2,l,p)
          zfp_temp(l) = coefs_np(3,l,p)

        elseif (ifcurved==1) then

          ! if (l==1) then
          !   xfp_temp(l) = xcglob(p)  
          !   yfp_temp(l) = ycglob(p)  
          !   zfp_temp(l) = zcglob(p)
          ! elseif (l==2) then
          !   xfp_temp(l) = xcglob(p)  
          !   yfp_temp(l) = ycglob(p)+ds  
          !   zfp_temp(l) = zcglob(p)      
          ! else
          !   xfp_temp(l) = xfp_temp(l-1) + ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)*(l-2)/(2.*NL)) 
          !   yfp_temp(l) = yfp_temp(l-1) + ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)*(l-2)/(2.*NL)) 
          !   zfp_temp(l) = zfp_temp(l-1) + ds*cos(theta_Spagetti(p))
          ! endif

        endif

        !        xfp_temp(l) = xcglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
        !        yfp_temp(l) = ycglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p))
        !        zfp_temp(l) = zcglob(p) + (l-1)*ds *cos(theta_Spagetti(p))   

        if ( (yfp_temp(l)+offset) .gt. backbound ) then
          if  ( (xfp_temp(l)-offset) .le. leftbound ) then 
            !if (typeflow=="d".and.leftbound==0.) goto 206
              if(count_slve_loc.eq. 0 ) i = i+1
              count_slve = count_slve + 1
              ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
              ap(i)%nb(2) = 1      !neighbor 2 of myid is particle's master
              count_slve_loc = count_slve_loc + 1
              goto 206

          endif
        endif  
    enddo

    206 continue 
    endif


    !neighbor 7 of particle's master
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax )
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) + 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid .eq. procrank) then

      do l=1,NL

        if (ifcurved==0) then

          !xfp_temp(l) = xcglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
          !yfp_temp(l) = ycglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)) 
          !zfp_temp(l) = zcglob(p) + (l-int((NL+1)/2.))*ds*cos(theta_Spagetti(p)) 

          !xfp_temp(l) = coefs(1,l)
          !yfp_temp(l) = coefs(2,l)
          !zfp_temp(l) = coefs(3,l)
          xfp_temp(l) = coefs_np(1,l,p)
          yfp_temp(l) = coefs_np(2,l,p)
          zfp_temp(l) = coefs_np(3,l,p)

        elseif (ifcurved==1) then

          ! if (l==1) then
          !   xfp_temp(l) = xcglob(p)  
          !   yfp_temp(l) = ycglob(p)  
          !   zfp_temp(l) = zcglob(p)
          ! elseif (l==2) then
          !   xfp_temp(l) = xcglob(p)  
          !   yfp_temp(l) = ycglob(p)+ds  
          !   zfp_temp(l) = zcglob(p)      
          ! else
          !   xfp_temp(l) = xfp_temp(l-1) + ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)*(l-2)/(2.*NL)) 
          !   yfp_temp(l) = yfp_temp(l-1) + ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)*(l-2)/(2.*NL)) 
          !   zfp_temp(l) = zfp_temp(l-1) + ds*cos(theta_Spagetti(p))
          ! endif

        endif

        !        xfp_temp(l) = xcglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
        !        yfp_temp(l) = ycglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p))
        !        zfp_temp(l) = zcglob(p) + (l-1)*ds *cos(theta_Spagetti(p))  

        if ( (yfp_temp(l)+offset) .gt. backbound ) then
          if ( ((xfp_temp(l)+offset) .ge. leftbound) .and. ((xfp_temp(l)-offset) .le. rightbound) ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            count_slve = count_slve + 1
            ap(i)%mslv = -p     !myid is slave of particle p=ap(i)%mslv
            ap(i)%nb(3) = 1      !neighbor 3 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 207

          endif
        endif

      enddo

    207 continue 
    endif


    !neighbor 8 of particle's master
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) + 1
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) + 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,ierr)
    if (myid .eq. procrank) then

    do l=1,NL

      if (ifcurved==0) then

        !xfp_temp(l) = xcglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
        !yfp_temp(l) = ycglob(p) + (l-int((NL+1)/2.))*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)) 
        !zfp_temp(l) = zcglob(p) + (l-int((NL+1)/2.))*ds*cos(theta_Spagetti(p)) 

        !xfp_temp(l) = coefs(1,l)
        !yfp_temp(l) = coefs(2,l)
        !zfp_temp(l) = coefs(3,l)
        xfp_temp(l) = coefs_np(1,l,p)
        yfp_temp(l) = coefs_np(2,l,p)
        zfp_temp(l) = coefs_np(3,l,p)

      elseif (ifcurved==1) then

        ! if (l==1) then
        !   xfp_temp(l) = xcglob(p)  
        !   yfp_temp(l) = ycglob(p)  
        !   zfp_temp(l) = zcglob(p)
        ! elseif (l==2) then
        !   xfp_temp(l) = xcglob(p)  
        !   yfp_temp(l) = ycglob(p)+ds  
        !   zfp_temp(l) = zcglob(p)      
        ! else
        !   xfp_temp(l) = xfp_temp(l-1) + ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)*(l-2)/(2.*NL)) 
        !   yfp_temp(l) = yfp_temp(l-1) + ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p)*(l-2)/(2.*NL)) 
        !   zfp_temp(l) = zfp_temp(l-1) + ds*cos(theta_Spagetti(p))
        ! endif

      endif

      !        xfp_temp(l) = xcglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*sin(phi_Spagetti(p)) 
      !        yfp_temp(l) = ycglob(p) + (l-1)*ds*sin(theta_Spagetti(p))*cos(phi_Spagetti(p))
      !        zfp_temp(l) = zcglob(p) + (l-1)*ds *cos(theta_Spagetti(p))   

      if ( (yfp_temp(l)+offset) .gt. backbound ) then
        if  ( (xfp_temp(l)+offset) .ge. rightbound ) then 
            !if (typeflow=="d".and.rightbound==lx) goto 208
            if(count_slve_loc.eq. 0 ) i = i+1
            count_slve = count_slve + 1
            ap(i)%mslv = -p     !myid is slave of particle p=ap(i)%mslv
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

!print*, "---------------------------------"
!do p=1,np
!if (ap(p)%mslv .ge. 0) then  
!print*,"p at which ap(p)%mslv > 0 ", p 
!  do l=1,NL
!	   print*, ap(p)%xfp(l), ap(p)%yfp(l), ap(p)%zfp(l)
!  enddo
!endif
!enddo
!print*, "---------------------------------"
!stop


! maximum number of particles in a thread is equal to the number of particles
! 'mastered' and 'enslaved' by it
!
!v pmax = i!count_mstr+count_slve
pmax = count_mstr+count_slve
npmstr = count_mstr




! print*,""    
! print*, "--------------------------------count_master", count_mstr
! print*, "--------------------------------pmax, npmstr", pmax, npmstr
! print*,""
 


do p=1,pmax   !v for all master and slave
  if (ap(p)%mslv .le. 0) then   !v among them for slave ones only
    ! myid is either slave of particle ap(p)%mslv or does not contain this particle

    do l=1,NL

      ap(p)%xfp(l)  = 0.
      ap(p)%yfp(l)  = 0.
      ap(p)%zfp(l)  = 0.
      ap(p)%qn1(l)  = 1.
      ap(p)%qn2(l)  = 0.
      ap(p)%qn3(l)  = 0.
      ap(p)%qn4(l)  = 0.

      ap(p)%xfpo(l) = 0.
      ap(p)%yfpo(l) = 0.
      ap(p)%zfpo(l) = 0.
      ap(p)%qno1(l) = 1.
      ap(p)%qno2(l) = 0.
      ap(p)%qno3(l) = 0.
      ap(p)%qno4(l) = 0.

      ap(p)%xfpold(l) = 0.  !v this I added additionally
      ap(p)%yfpold(l) = 0.
      ap(p)%zfpold(l) = 0.
      ap(p)%qnold1(l) = 1.
      ap(p)%qnold2(l) = 0.
      ap(p)%qnold3(l) = 0.
      ap(p)%qnold4(l) = 0.      

      ! ap(p)%dxdt(l) = 0.   !v added to check if this information is stored like xfp in new time step
      ! ap(p)%dydt(l) = 0.
      ! ap(p)%dzdt(l) = 0.
      ! ap(p)%omg1(l) = 0.
      ! ap(p)%omg2(l) = 0.
      ! ap(p)%omg3(l) = 0.

      ! ap(p)%ua(l)    = 0. !v added to check if this information is stored like xfp in new time step
      ! ap(p)%va(l)    = 0.
      ! ap(p)%wa(l)    = 0.
      ! ap(p)%omgd1(l) = 0.
      ! ap(p)%omgd2(l) = 0.
      ! ap(p)%omgd3(l) = 0.
      
    enddo

  endif
enddo



ap(1:npmax)%integralx = 0.
ap(1:npmax)%integraly = 0.
ap(1:npmax)%integralz = 0.
forall (p=1:npmax)         !v npmax is for max number of particles 
  ap(p)%dxdt(:) = 0.
  ap(p)%dydt(:) = 0.
  ap(p)%dzdt(:) = 0.
  ap(p)%omg1(:) = 0.
  ap(p)%omg2(:) = 0.
  ap(p)%omg3(:) = 0.

  ap(p)%ua(:)    = 0.
  ap(p)%va(:)    = 0.
  ap(p)%wa(:)    = 0.
  ap(p)%omgd1(:) = 0.
  ap(p)%omgd2(:) = 0.
  ap(p)%omgd3(:) = 0.

  ap(p)%fxl(:)  = 0.
  ap(p)%fyl(:)  = 0.
  ap(p)%fzl(:)  = 0.
  ap(p)%mxl(:)  = 0.
  ap(p)%myl(:)  = 0.
  ap(p)%mzl(:)  = 0.

  ap(p)%ul(:)   = 0.
  ap(p)%vl(:)   = 0.
  ap(p)%wl(:)   = 0.
  ap(p)%omgl1(:) = 0.
  ap(p)%omgl2(:) = 0.
  ap(p)%omgl3(:) = 0.

  ap(p)%dxdto(:) = 0. !v I added
  ap(p)%dydto(:) = 0.
  ap(p)%dzdto(:) = 0.
  ap(p)%omgo1(:) = 0.
  ap(p)%omgo2(:) = 0.
  ap(p)%omgo3(:) = 0.

  ap(p)%uao(:)    = 0. !v I added
  ap(p)%vao(:)    = 0.
  ap(p)%wao(:)    = 0.
  ap(p)%omgdo1(:) = 0.
  ap(p)%omgdo2(:) = 0.
  ap(p)%omgdo3(:) = 0.
  
end forall

forall (p=pmax+1:npmax)  !v for master+slve to npmax (max numer of particles)
  ap(p)%xfp(:)  = 0. 
  ap(p)%yfp(:)  = 0.
  ap(p)%zfp(:)  = 0.
  ap(p)%qn1(:)  = 1.
  ap(p)%qn2(:)  = 0.
  ap(p)%qn3(:)  = 0.
  ap(p)%qn4(:)  = 0.

  ap(p)%xfpo(:) = 0.
  ap(p)%yfpo(:) = 0.
  ap(p)%zfpo(:) = 0.
  ap(p)%qno1(:) = 1.
  ap(p)%qno2(:) = 0.
  ap(p)%qno3(:) = 0.
  ap(p)%qno4(:) = 0.

  ap(p)%xfpold(:) = 0.
  ap(p)%yfpold(:) = 0.
  ap(p)%zfpold(:) = 0.
  ap(p)%qnold1(:) = 1.
  ap(p)%qnold2(:) = 0.
  ap(p)%qnold3(:) = 0.
  ap(p)%qnold4(:) = 0.

  ! ap(p)%dxdt(:) = 0.   !v added to check if it initializes properly in next timestep as xfp
  ! ap(p)%dydt(:) = 0.
  ! ap(p)%dzdt(:) = 0.
  ! ap(p)%omg1(:) = 0.
  ! ap(p)%omg2(:) = 0.
  ! ap(p)%omg3(:) = 0.

  ! ap(p)%ua(:)    = 0. !v added to check if it initializes properly in next timestep as xfp
  ! ap(p)%va(:)    = 0.
  ! ap(p)%wa(:)    = 0.
  ! ap(p)%omgd1(:) = 0.
  ! ap(p)%omgd2(:) = 0.
  ! ap(p)%omgd3(:) = 0.
  
end forall




write(6,'(A7,I5,A8,I5,A18,I5,A11,A8,I5)') 'Thread ', myid, ' masters ', count_mstr, &
' and is slave for ', count_slve, ' particles. ', ' pmax = ', pmax

!
call MPI_ALLREDUCE(count_mstr,count_mstr_all,1,MPI_INTEGER,MPI_SUM,comm_cart,ierr)
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


!--- Coordinate arrays (coords are already stored in ap(p)%xfp etc)
! this is used to later extract that information in xn (simplifies things)
! NOTE: ideally Xi should be done for each particle, as J_s; J_b, etc can 
!       be differnt for other than straght fibers. But, fine for straight
!       fibers.
!if (myid .eq. 0) then
if (myid.eq.procrank) then
  do p=1,np
    do i=1, nno
      Xi(i,1:3,p) = coefs_np(1:3,i,p)
      Xi(i,4:7,p) = (/ 1., 0., 0., 0. /)
      print '(10E10.4)', Xi(i,1:3,p)
    end do   
    !xn=Xi
  enddo
  print*,""
  print*,""
endif


! Control Point  Connectivity Arrays
!if (myid .eq. 0) then
if (myid.eq.procrank) then
do j = 1, p1_fordr
  CONS(1,j) = j 
end do
!print*, CONS(1,:)
do i = 2, nel
  CONS(i,:) = CONS(i-1,:) + 1 + add_ele
	!print*, CONS(i,:)
end do
!print*,"-----------end of CONS------------"
endif
!print*, "nno, nl, noCP, p1_fordr", nno, nl, noCP, p1_fordr
!print*,"ndof:            ", ndof


! Knot connectivity array (note: only for non-repeated internal knot entries)
!if (myid .eq. 0) then
if (myid.eq.procrank) then
unqKnot = linspace(0.0_dp, 1.0_dp, nel+1)
!unqKnot = linspace(0.0_dp, 1.0_dp, nel+1)
!unqKnot = !uKnot(p1_fordr : ubound(uKnot,1)-p1_fordr)
do i = 1, nel
  do j = 1, 2
     elRangeU(i, j) = unqKnot(i+j-1)
  end do
end do
endif

!print*, "unqKnot:    ", real(unqKnot,4)
!print*, "elRangeU:   ", real(elRangeU,4)
!stop

if (myid .eq. procrank) then
	!--- compute quad points & weights
	call quad_rule(ngps, xigs, wgs)
	call quad_rule(ngpb, xigb, wgb)
	
	
	call setupInitial( nno, nel, CONS, Xi, elRangeU, xigs, wgs, xigb, &
			 wgb, ngps, ngpb, J_s, J_b, rhof, Arf, I_yy, I_zz, J_xx, p1_fordr, &
			 NdN_array, NdN_arrayb, mInerCP, uKnot, nderiv, mInerCP2, Eflag, &
			 rhoA, rhoI_yy, rhoI_zz, np, strain_k0)

			 !do i  = 1,nel*ngps
			 ! 	print '(3E12.4)', (strain_k0(i,j,1), j=1,3)
			 !enddo		
endif		 
!Output: J_s(nel*ngps), J_b(nel*ngpb), mInerCP(6x6)
! !        NdN_array, NdN_arrayb(nel*ngps*p1_fordr, 2)

!
!if (myid.eq.procrank) then
!	print*,"------Js", real(J_s,4)
!	print*,"------Jb", real(J_b,4)
!	print '(6E15.4)',(mInerCP(1,1)), (mInerCP(2,2)), (mInerCP(3,3)), (mInerCP(4,4)), (mInerCP(5,5)), (mInerCP(6,6)) 
!endif
!stop




!counter = 1
!do i= 1,nel

!  do gp = 1,ngps
!
!     do icp = 1,p1_fordr
!        do jcp = 1, nderiv
!           NdN_array((i-1)*ngps*p1_fordr+icp,jcp) = dNshpfun(icp,jcp)
!           !NdN_array((counter-1)*p1_fordr+icp,jcp) = dNshpfun(icp,jcp)
!        end do
!        print*, NdN_array(icp,:)
!     end do
!     print*,",.,.,.,.,.,.,.,.,.,.,.,.,."
!
!print '(5E10.3)', strain_k0(counter,:,1)
!counter=counter+1
!enddo
!stop
!do p = 9,16!61_fordr*2
!    print*, NdN_arrayb(p,:)! = dNshpfun(icp,jcp)
!enddo		
!stop

!print*,'  ----ap coords ---'
! do pp=1,pmax
! 	if (ap(pp)%mslv .ge. 0) then  
! 		!if (myid .eq. 0) then
! 		!write(number,'(i3.3)') ap(pp)%mslv
! 		!write(filenumber,'(i7.7)') 
! 		!open(42,file='data/Filament_initial'//'.txt')
! 		do j=1,NL
! 			!write(42,'(3E15.7)')  ap(pp)%xfp(j),ap(pp)%yfp(j),ap(pp)%zfp(j)
!  	  !print '(3E15.7)', ap(pp)%xfp(l),ap(pp)%yfp(l),ap(pp)%zfp(l)
! 	  enddo
! 		close(42)
! 	endif
! enddo
!stop

!print*, "---------------------------------",(size(coefs_np,2))
!do p = 1,np
!	if (ap(p)%mslv .ge. 0) then                                
!		do i=1,nno                                                
!       !do j=1,(size(coefs_np,2))                          
!					!print '(6E12.4)', (coefs_np(i,j,p), j=1,NL)
!					print '(3E12.4)', (Xi(i,j,p), j=1,3)
!       !enddo
!    enddo
!	endif
!	print*, "--------------- looping over p------------------"
!	print*,""
!enddo	
!STOP

!
return
end subroutine initparticles
!
end module initparticles_fibm
