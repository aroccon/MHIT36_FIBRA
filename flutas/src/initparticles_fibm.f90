module initparticles_fibm
use decomp_2d
use mod_param_fibm,  only : dims, ndims, lx, ly, lz, ds, &
                        dx, dy, dz, nderiv, &
                        L_fibr, bc0_dof, bcL_dof, ndof, &
                        np, pi, fboundary, &
                        ifattached, ifcurved, nrow, &
                        ncolumn, deg_ele, add_ele, noCP,  &
                        uKnoti, nel, nl, offset, PI, &
                        nno, p_ordr, p1_fordr, ngps, c_type, &
                        ngpb, rhof, radius, datadir, npi, npj, &
												npk, Eflag, rhoA, rhoI_yy, rhoI_zz, diametr, &
												nxie, nxi_tvec, nxi_vec
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
use mod_initialSetup,  only: setupInitial, ap_ugp, ap_fgp
use mod_linspace
use bspline



implicit none
private
public initparticles
contains
subroutine initparticles
implicit none
real, dimension(np) :: xcglob,ycglob,zcglob,thetacglob,phicglob
integer :: i,j,k,p,pp,rk
integer :: proccoords(1:ndims),procrank
integer, dimension(2) :: sbuf,rbuf
real :: leftbound,rightbound,frontbound,backbound
real :: dist,distx,disty,distz,distzw,angle
real :: xp,yp
real :: ax
real :: ay
real :: rn
integer :: counter,crys
integer :: idp
character(len=5) rankpr
character(len=7) :: tempstr
character(len=11) :: tempstr2
integer :: count_mstr,count_slve,count_mstr_all,count_slve_loc
real, dimension(nno) :: xfp_temp,yfp_temp,zfp_temp,qn1_temp,qn2_temp,qn3_temp,qn4_temp,theta
integer :: ifil,jfil,kfil
!v-----------------------------------------------------------------------------------------------
REAL(8)  :: coefsi(4,2)
REAL(8)  :: coefsi_np(4,2,np), coefsi_temp(4,2)														 

coefs = RESHAPE( (/ 0.0, 0.0, 0.0, 1.0, L_fibr, 0.0, 0.0, 1.0/), (/ 4, 2 /) )
!coefsi = RESHAPE( (/ 0.5*lx, 2.0, 4.0, 1.0, &
!					0.5*lx, 2.0+L_fibr*cos(0.1*pi), 4.0+L_fibr*sin(0.1*pi), 1.0/), (/ 4, 2 /) )

do p=1,np
  do i=1,4
    do j=1,2
      coefsi_np(i,j,p) = coefsi(i,j)
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

! compute area and second moment of area
if (c_type .eq. "c") then
	Arf  = PI * (0.5 * diametr)**2
	I_zz = (PI/4.0) * (0.5*diametr)**4
  I_yy = I_zz
elseif (c_type .eq. "s") then
	Arf  = diametr**2
	I_zz = (diametr**4)/12.0
  I_yy = I_zz
endif
J_xx = I_zz + I_yy

! location of points in parametric system
nxi_tvec = linspace(0.0_dp, L_fibr*1.0_dp, 2+nl)
nxi_vec  = nxi_tvec(2:1+nl)



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
	   jfil = 0
  else

    ! pseudo-random
    counter = 0
    do p=1,NP
      if(NP<=5) then

        if (ifattached==0) then

          jfil=int((p-1)/(npi*npk))+1
          ifil=int(((p-((jfil-1)*npi*npk))-1)/real(npk))+1
          kfil=(p-((jfil-1)*npi*npk))-((ifil-1)*npk)

          xcglob(p)     = lx/2.
          ycglob(p)    = 0.95*ly  
          zcglob(p)    = 0.5*lz 

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
   coefsi_np(1,1:2,p) = xcglob(p)
   coefsi_np(2,1:2,p) = ycglob(p)
   coefsi_np(3,1:2,p) = (/ zcglob(p)-0.5*L_fibr, zcglob(p)+0.5*L_fibr /)

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
     	 enddo

			 CALL DegreeElevate(4, size(coefsi_np,2)-1, 1, uknoti, &
                          coefsi_temp, deg_ele, noCP-1, uknotp, coefsp_temp)

			 do i=1,4
    	    do j=1,(size(coefsp_temp,2))
    	       coefsp_np(i,j,p) = coefsp_temp(i,j)
    	    enddo
    	 enddo
													 
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

  		 do i=1,4
					do j=1,size(coefsp_np,2)
						 coefsp_np(i,j,p) = coefsi_np(i,j,p)
					enddo
			 enddo	

	   endif
	 
	   !------  Refine Mesh
		 if (p .eq. 1)  then 
				ALLOCATE(vec_dummy(nel+1), kntins(nel-1))
				vec_dummy = linspace(0.0_dp, 1.0_dp, nel+1)
				kntins    = vec_dummy(2 : nel)    
				ALLOCATE( uKnoth(size(uKnotp) + size(kntins) ),  & 
					 coefsh_temp(size(coefsp_np,1), size(coefsp_np,2)+size(kntins) ), &
	         coefsh_np(size(coefsp_np,1), size(coefsp_np,2)+size(kntins), np))
		 endif			 
		
		 !--- initialize
		 coefsp_temp = 0.0D+00
		 coefsh_temp = 0.0D+00
		 do i=1,4
        do j=1,(size(coefsp_np,2))
           coefsp_temp(i,j) = coefsp_np(i,j,p)
        enddo
     enddo
	
	   CALL RefineKnotVector(4, size(coefsp_np,2)-1, p_ordr, &   
         uKnotp, coefsp_temp, size(kntins,1)-1, kntins, uKnoth, coefsh_temp)
     
		 do i=1,4
        do j=1,(size(coefsh_np,2))
           coefsh_np(i,j,p) = coefsh_temp(i,j)
        enddo
     enddo


		 if (p .eq. 1) then
			 noCP = size(coefsh_np,2)
		 endif

	   !------ Additional order elevation 
		 if (add_ele .gt. 0) then
			 
			 if (p .eq. 1)  then 
        unqKnot = uKnoth(p_ordr+1 : ubound(uKnoth,1)-(1+deg_ele))
				ALLOCATE( uKnot(size(uKnoth) + add_ele*size(unqKnot) ),  & 
					 coefs_temp(size(coefsh_np,1), size(coefsh_np,2)+add_ele*nel ), &
	         coefs_np(size(coefsh_np,1), size(coefsh_np,2)+add_ele*nel, np))
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
		 	 enddo		

			 if (p .eq. 1) then
	        noCP   = size(coefs_np,2)
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
			 enddo	
		endif										
	 

	 elseif (ifcurved==1) then
	 
	   print*, "---------------------------------(ifcurved==1) "
	   stop
	 endif
enddo




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
i  = 0
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
    
    
    do l=1,nno


      if (ifcurved==0) then

        ap(i)%xll(l) = coefs_np(1,l,p)!ap(i)%x + (l-int((NL+1)/2.))*ds
        ap(i)%yll(l) = coefs_np(2,l,p)!ap(i)%y + (l-int((NL+1)/2.))*ds * cos(0.1*pi)
        ap(i)%zll(l) = coefs_np(3,l,p)!ap(i)%z + (l-int((NL+1)/2.))*ds * sin(0.1*pi)
        ap(i)%qn1(l) = 1.0
        ap(i)%qn2(l) = 0.0
        ap(i)%qn3(l) = 0.0
        ap(i)%qn4(l) = 0.0

      endif

      ap(i)%xfpo(l) = ap(i)%xll(l)
      ap(i)%yfpo(l) = ap(i)%yll(l)
      ap(i)%zfpo(l) = ap(i)%zll(l)    
      ap(i)%qno1(l) = ap(i)%qn1(l)
      ap(i)%qno2(l) = ap(i)%qn2(l)
      ap(i)%qno3(l) = ap(i)%qn3(l)
      ap(i)%qno4(l) = ap(i)%qn4(l)


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

      do l=1,nno


        if (ifcurved==0) then

          xfp_temp(l) = coefs_np(1,l,p)
          yfp_temp(l) = coefs_np(2,l,p)
          zfp_temp(l) = coefs_np(3,l,p)

        endif

        if ( (xfp_temp(l)+offset) .gt. rightbound ) then
          if ( ((yfp_temp(l)+offset) .ge. frontbound) .and. ((yfp_temp(l)-offset) .le. backbound) ) then 
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

      do l=1,nno

        if (ifcurved==0) then

          xfp_temp(l) = coefs_np(1,l,p)
          yfp_temp(l) = coefs_np(2,l,p)
          zfp_temp(l) = coefs_np(3,l,p)

        endif

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

      do l=1,nno

        if (ifcurved==0) then

          xfp_temp(l) = coefs_np(1,l,p)
          yfp_temp(l) = coefs_np(2,l,p)
          zfp_temp(l) = coefs_np(3,l,p)

        endif

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

      do l=1,nno

        if (ifcurved==0) then

          xfp_temp(l) = coefs_np(1,l,p)
          yfp_temp(l) = coefs_np(2,l,p)
          zfp_temp(l) = coefs_np(3,l,p)

          endif

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

      do l=1,nno

        if (ifcurved==0) then

          xfp_temp(l) = coefs_np(1,l,p)
          yfp_temp(l) = coefs_np(2,l,p)
          zfp_temp(l) = coefs_np(3,l,p)

        endif


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

      do l=1,nno

        if (ifcurved==0) then

          xfp_temp(l) = coefs_np(1,l,p)
          yfp_temp(l) = coefs_np(2,l,p)
          zfp_temp(l) = coefs_np(3,l,p)

        endif

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

      do l=1,nno

        if (ifcurved==0) then

          xfp_temp(l) = coefs_np(1,l,p)
          yfp_temp(l) = coefs_np(2,l,p)
          zfp_temp(l) = coefs_np(3,l,p)

        endif

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

    do l=1,nno

      if (ifcurved==0) then

        xfp_temp(l) = coefs_np(1,l,p)
        yfp_temp(l) = coefs_np(2,l,p)
        zfp_temp(l) = coefs_np(3,l,p)

      endif

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


! maximum number of particles in a thread is equal to the number of particles
! 'mastered' and 'enslaved' by it
!
pmax = count_mstr+count_slve
npmstr = count_mstr




 
do p=1,pmax   !v for all master and slave
  if (ap(p)%mslv .le. 0) then   !v among them for slave ones only
    ! myid is either slave of particle ap(p)%mslv or does not contain this particle

    do l=1,nno

      ap(p)%xll(l)  = 0.
      ap(p)%yll(l)  = 0.
      ap(p)%zll(l)  = 0.
      ap(p)%qn1(l)  = 0.
      ap(p)%qn2(l)  = 0.
      ap(p)%qn3(l)  = 0.
      ap(p)%qn4(l)  = 0.

      ap(p)%xfpo(l) = 0.
      ap(p)%yfpo(l) = 0.
      ap(p)%zfpo(l) = 0.
      ap(p)%qno1(l) = 0.
      ap(p)%qno2(l) = 0.
      ap(p)%qno3(l) = 0.
      ap(p)%qno4(l) = 0.

      ap(p)%xfpold(l) = 0.  
      ap(p)%yfpold(l) = 0.
      ap(p)%zfpold(l) = 0.
      ap(p)%qnold1(l) = 0.
      ap(p)%qnold2(l) = 0.
      ap(p)%qnold3(l) = 0.
      ap(p)%qnold4(l) = 0.      

    enddo

  endif
enddo



ap(1:npmax)%integralx = 0.
ap(1:npmax)%integraly = 0.
ap(1:npmax)%integralz = 0.
forall (p=1:npmax)         
  ap(p)%dxdtl(:) = 0.
  ap(p)%dydtl(:) = 0.
  ap(p)%dzdtl(:) = 0.
  ap(p)%omg1(:)  = 0.
  ap(p)%omg2(:)  = 0.
  ap(p)%omg3(:)  = 0.

  ap(p)%ua(:)    = 0.
  ap(p)%va(:)    = 0.
  ap(p)%wa(:)    = 0.
  ap(p)%omgd1(:) = 0.
  ap(p)%omgd2(:) = 0.
  ap(p)%omgd3(:) = 0.

  ap(p)%fxll(:)  = 0.
  ap(p)%fyll(:)  = 0.
  ap(p)%fzll(:)  = 0.
  ap(p)%mxl(:)  = 0.
  ap(p)%myl(:)  = 0.
  ap(p)%mzl(:)  = 0.

  ap(p)%ull(:)   = 0.
  ap(p)%vll(:)   = 0.
  ap(p)%wll(:)   = 0.
  ap(p)%omgl1(:) = 0.
  ap(p)%omgl2(:) = 0.
  ap(p)%omgl3(:) = 0.

  ap(p)%dxdto(:) = 0. 
  ap(p)%dydto(:) = 0.
  ap(p)%dzdto(:) = 0.
  ap(p)%omgo1(:) = 0.
  ap(p)%omgo2(:) = 0.
  ap(p)%omgo3(:) = 0.

  ap(p)%uao(:)    = 0. 
  ap(p)%vao(:)    = 0.
  ap(p)%wao(:)    = 0.
  ap(p)%omgdo1(:) = 0.
  ap(p)%omgdo2(:) = 0.
  ap(p)%omgdo3(:) = 0.
  
end forall

forall (p=pmax+1:npmax)  
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

end forall




write(6,'(A7,I5,A8,I5,A18,I5,A11,A8,I5)') 'Thread ', myid, ' masters ', count_mstr, &
' and is slave for ', count_slve, ' particles. ', ' pmax = ', pmax

!
call MPI_ALLREDUCE(count_mstr,count_mstr_all,1,MPI_INTEGER,MPI_SUM,comm_cart,ierr)

!--- Coordinate arrays (coords are already stored in ap(p)%xfp etc)
! this is used to later extract that information in xn (simplifies things)
! NOTE: ideally Xi should be done for each particle, as J_s; J_b, etc can 
!       be differnt for other than straght fibers. But, fine for straight
!       fibers.
do p=1,np
  do i=1, nno
    Xi(i,1:3,p) = coefs_np(1:3,i,p)
    Xi(i,4:7,p) = (/ 1., 0., 0., 0. /)
  end do   
enddo


! Control Point  Connectivity Arrays
do j = 1, p1_fordr
  CONS(1,j) = j 
end do
do i = 2, nel
  CONS(i,:) = CONS(i-1,:) + 1 + add_ele
end do


! Knot connectivity array (note: only for non-repeated internal knot entries)
unqKnot = linspace(0.0_dp, 1.0_dp, nel+1)
do i = 1, nel
	do j = 1, 2
		 elRangeU(i, j) = unqKnot(i+j-1)
	end do
end do


call quad_rule(ngps, xigs, wgs)
call quad_rule(ngpb, xigb, wgb)

call setupInitial( nno, nel, CONS, Xi, elRangeU, xigs, wgs, xigb, &
		 wgb, ngps, ngpb, J_s, J_b, rhof, Arf, I_yy, I_zz, J_xx, p1_fordr, &
		 NdN_array, NdN_arrayb, mInerCP, uKnot, nderiv, mInerCP2, Eflag, &
		 rhoA, rhoI_yy, rhoI_zz, np, strain_k0)

!
return
end subroutine initparticles
!
end module initparticles_fibm
