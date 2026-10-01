module mod_initialSetup
use bspline

contains

subroutine setupInitial(nno, nel, CONS, Xi, elRangeU, xigs, wgs, xigb, wgb, &
           ngps, ngpb, J_s, J_b, rhof, Arf, I_yy, I_zz,  J_xx, & 
           p1_fordr, NdN_array, NdN_arrayb, mInerCP, uKnot, nderiv, mInerCP2, &
					 Eflag, rhoA, rhoI_yy, rhoI_zz, np, strain_k0)


use mod_param,         only: rho_sp
					 
IMPLICIT NONE

INTEGER, PARAMETER   :: rk=kind(1.0D+00) !no. of derv + shpfun
INTEGER, INTENT(IN)  :: nno, nel, p1_fordr, ngps, ngpb, &
                        CONS(nel, p1_fordr), nderiv, np
!REAL(16),INTENT(IN) :: Xi(nno, 7), xigs(ngps), wgs(ngps),&
REAL(16), INTENT(IN) :: Xi(nno,7,np), xigs(ngps), wgs(ngps),&
                        xigb(ngpb), wgb(ngpb)
REAL, INTENT(IN)     :: elRangeU(nel, 2), rhof
REAL(8), INTENT(IN)  :: uKnot(nno+p1_fordr)
REAL(16),INTENT(IN)  :: Arf, I_yy, I_zz, J_xx
LOGICAL, INTENT(IN)  :: Eflag
REAL(16),INTENT(IN)  :: rhoA, rhoI_yy, rhoI_zz
!REAL(16), INTENT(INOUT) :: J_s(nel*ngps), J_b(nel*ngpb), &
REAL(16), INTENT(INOUT) :: J_s(nel*ngps,np), J_b(nel*ngpb,np), &
                        NdN_array (nel*ngps*p1_fordr, nderiv), &
                        NdN_arrayb(nel*ngpb*p1_fordr, nderiv), &
                        mInerCP(6,6), mInerCP2(6,6)
												!dN_array(nel * ngps * p1_fordr), &  
                        !dN_arrayb(nel * ngps * p1_fordr), &                      
                        !N_array(nel * ngps * p1_fordr),  &
                        !N_arrayb(nel * ngps * p1_fordr), &
REAL(8),INTENT(INOUT):: strain_k0(nel*ngps,3,np)
LOGICAL              :: mflag 
INTEGER(4)           :: cone(p1_fordr), i, j, gp, icp, jcp, ispan, &
                        counter, inx, p
REAL                 :: eta, elU(2), c1, c2, N(p1_fordr), Jsu
REAL(8)              :: dNshpfun(p1_fordr,nderiv), uu
REAL(16)             :: state_x0(p1_fordr,3), dr0(3) 
!---------------------------------------------------------------------!

do p= 1,np

	counter = 1
	do i= 1,nel
	  cone     = CONS(i,1:p1_fordr)
	  state_x0 = Xi(cone,1:3, p)
	  elU      = elRangeU(i,1:2)
	  c1       = 0.5 * (elU(2) - elU(1))
	  c2       = 0.5 * (elU(2) + elU(1))
	  
	  do gp = 1,ngps
	     eta   = xigs(gp)
	     uu    = c1*eta + c2
	     
	     ispan = FindSpan(nno, p1_fordr-1, uu, uKnot) 
	     call DersBasisFuns(ispan, uu, p1_fordr-1, nderiv-1, uKnot, dNshpfun) 
	     !call BasisFuns(ispan, uu, p1_fordr-1, uKnot, Nshpfun) 
	 
		   inx = (i-1)*ngps*p1_fordr + (gp-1)*p1_fordr
			 do icp = 1,p1_fordr
     	    NdN_array(inx+icp,:) = dNshpfun(icp,:)
     	 end do

	     !do icp = 1,p1_fordr
	     !   do jcp = 1, nderiv
	     !      NdN_array((counter-1)*p1_fordr+icp,jcp) = dNshpfun(icp,jcp)
	     !   end do
	     !end do
	
	     dr0 = 0.0 
	     do icp = 1,p1_fordr
	        dr0 = dr0 + state_x0(icp,1:3)*dNshpfun(icp,2)
	     end do
       Jsu = SQRT( dr0(1)*dr0(1) + dr0(2)*dr0(2) + dr0(3)*dr0(3) )
	     J_s(counter,p) = Jsu

			 ! compute initial strain
			 strain_k0(counter,1:3,p) = dr0/Jsu
	     
	     counter = counter + 1     
	  end do  
	end do
enddo


if (ngpb .NE. ngps) then
	
	do p = 1,np
	   counter = 1
	   do i=1,nel
	      cone     = CONS(i,1:p1_fordr)
	      state_x0 = Xi(cone,1:3,p)
	      elU      = elRangeU(i,1:2)
	      c1       = 0.5 * (elU(2) - elU(1))
	      c2       = 0.5 * (elU(2) + elU(1))
	      
	      do gp = 1,ngpb
	         eta   = xigb(gp)
	         uu    = c1*eta + c2
	         
	         ispan = FindSpan(nno, p1_fordr-1, uu, uKnot) 
	         call DersBasisFuns(ispan, uu, p1_fordr-1, nderiv-1, uKnot, dNshpfun)
	         
					 inx = (i-1)*ngpb*p1_fordr + (gp-1)*p1_fordr
         	 do icp = 1,p1_fordr
         	    NdN_arrayb(inx+icp,:) = dNshpfun(icp,:)
         	 end do
	         !do icp = 1,p1_fordr
	         !   do jcp = 1, nderiv
	         !      NdN_arrayb((counter-1)*p1_fordr+icp,jcp) = dNshpfun(icp,jcp)
	         !   end do
	         !end do
	
	         dr0 = 0.0 !(/ 0. 0. 0. /)     
	         do icp = 1,p1_fordr
	            dr0 = dr0 + state_x0(icp,1:3)*dNshpfun(icp,2)
	         end do
	         J_b(counter,p) = SQRT( dr0(1)*dr0(1) + dr0(2)*dr0(2) + dr0(3)*dr0(3) )
	         
	         counter = counter + 1  
	      end do  
	   end do
	enddo	 

else
	 !do p = 1,np
   J_b = J_s
   NdN_arrayb = NdN_array
end if


!!--- compute vec_ex
!do p = 1,np
!      cone     = CONS(1,1:p1_fordr)
!      state_x0 = Xi(cone,1:3,p)
!      elU      = elRangeU(1,1:2)
!      c1       = 0.5 * (elU(2) - elU(1))
!      c2       = 0.5 * (elU(2) + elU(1))
!      
!      do gp = 1,1 
!         eta   = xigb(gp)
!         uu    = c1*eta + c2
!         
!         ispan = FindSpan(nno, p1_fordr-1, uu, uKnot) 
!         call DersBasisFuns(ispan, uu, p1_fordr-1, nderiv-1, uKnot, dNshpfun)
!         
!         do icp = 1,p1_fordr
!            do jcp = 1, nderiv
!               NdN_arrayb((counter-1)*p1_fordr+icp,jcp) = dNshpfun(icp,jcp)
!            end do
!         end do
!
!         dr0 = 0.0 !(/ 0. 0. 0. /)     
!         do icp = 1,p1_fordr
!            dr0 = dr0 + state_x0(icp,1:3)*dNshpfun(icp,2)
!         end do
!         J_b(counter,p) = SQRT( dr0(1)*dr0(1) + dr0(2)*dr0(2) + dr0(3)*dr0(3) )
!         
!         counter = counter + 1  
!      end do  
!enddo	 


!--- Compute Mass Inertia Matrix 
if (Eflag .eqv. .false.) then

	mInerCP = 0.0
	mInerCP(1,1) = rhof * Arf
	mInerCP(2,2) = rhof * Arf
	mInerCP(3,3) = rhof * Arf
	mInerCP(4,4) = rhof * J_xx 
	mInerCP(5,5) = rhof * I_yy
	mInerCP(6,6) = rhof * I_zz
	
	mInerCP2 = 0.0
	mInerCP2(1,1) = (rhof - rho_sp) * Arf
	mInerCP2(2,2) = (rhof - rho_sp) * Arf
	mInerCP2(3,3) = (rhof - rho_sp) * Arf
	mInerCP2(4,4) = rhof * J_xx
	mInerCP2(5,5) = rhof * I_yy
	mInerCP2(6,6) = rhof * I_zz

else

	mInerCP      = 0.0
	mInerCP(1,1) = rhoA
	mInerCP(2,2) = rhoA
	mInerCP(3,3) = rhoA
	mInerCP(4,4) = rhoI_yy + rhoI_zz 
	mInerCP(5,5) = rhoI_yy
	mInerCP(6,6) = rhoI_zz
	
	mInerCP2      = 0.0
	mInerCP2(1,1) = rhoA - rho_sp*Arf
	mInerCP2(2,2) = rhoA - rho_sp*Arf
	mInerCP2(3,3) = rhoA - rho_sp*Arf
	mInerCP2(4,4) = rhof * (I_yy+I_zz) 
	mInerCP2(5,5) = rhoI_yy
	mInerCP2(6,6) = rhoI_zz
end if


!print*,"-------m_CP-----------"
!print "(E10.2)", mInerCP
!
!print*,""
!print*,"-------m_CP2-----------"
!print "(E10.2)", mInerCP2
!stop

!do i =1, nel*ngps*p1_fordr
!	write(*,'(F10.4)') (Ndn_arrayb(i,j), j=1,nderiv)
!	print*,''
!enddo
!print*,'------end of print Ndn_array----'

end subroutine setupInitial




!--------- In: ap(pp)%dxdtl at CPs, Out: ap(pp)%dxdt, ap(pp)%xfp(nxi) at xip
subroutine ap_ugp
use mod_param_fibm,  only: nxie, nl, nxi_vec, nel, nderiv, &
											 p1_fordr, nno
use mod_common_fibm,  only: ap, pmax, gp, ispan, uKnot, CONS, &
											 dNshpfun, uu, icp, jcp, cone

IMPLICIT NONE

INTEGER, PARAMETER   :: rk=kind(1.0D+00) !no. of derv + shpfun
INTEGER              :: pp, ii, jj, nxi
REAL                 :: state_x0(p1_fordr,3)
!-------------------------------------------

do pp= 1,pmax
  if (ap(pp)%mslv .gt. 0) then

	! initializing vel and coords
	ap(pp)%xfp(1:nl)  = 0.d0
	ap(pp)%yfp(1:nl)  = 0.d0
	ap(pp)%zfp(1:nl)  = 0.d0
	ap(pp)%dxdt(1:nl) = 0.d0
	ap(pp)%dydt(1:nl) = 0.d0
	ap(pp)%dzdt(1:nl) = 0.d0

	do ii = 1,nel
	  cone = CONS(ii,1:p1_fordr)
		!do jj = 1,p1_fordr
		!ie_con(jj,1:3)  = (/ 6*cone(jj)-5, 6*cone(jj)-4, 6*cone(jj)-3 /) !define ie_con(p1_fordr,3) 					
		!state_x0(jj,1:3)= (/ ap(pp)%xfp(cone(jj)), ap(pp)%yfp(cone(jj)), ap(pp)%zfp(cone(jj)) /) !state_x0(p1_fordr,3)
		!dxdtE(jj,1:3)   = (/ap(pp)%dxdt(cone(jj)), ap(pp)%dydt(cone(jj)), ap(pp)%dzdt(cone(jj)) /) !define dxdtE(p1_fordr,3)
		!enddo
	  
		do gp = 1,nxie   

			nxi = (ii-1)*nxie+gp
			uu  = nxi_vec(nxi) !real :: xi 
		  
			!note: here basisfun can be used instead of DersBasisFuns
		  ispan = FindSpan(nno, p1_fordr-1, uu, uKnot)
			call DersBasisFuns(ispan, uu, p1_fordr-1, nderiv-1, & 
		                                              uKnot, dNshpfun)
      
			! compute velocity and coordinates at xi points
			do jj = 1,p1_fordr
				 ap(pp)%dxdt(nxi)= ap(pp)%dxdt(nxi)+ dNshpfun(jj,1)*ap(pp)%dxdtl(cone(jj))
				 ap(pp)%dydt(nxi)= ap(pp)%dydt(nxi)+ dNshpfun(jj,1)*ap(pp)%dydtl(cone(jj))
				 ap(pp)%dzdt(nxi)= ap(pp)%dzdt(nxi)+ dNshpfun(jj,1)*ap(pp)%dzdtl(cone(jj))
        
				 ap(pp)%xfp(nxi)= ap(pp)%xfp(nxi)+ dNshpfun(jj,1)*ap(pp)%xll(cone(jj))
				 ap(pp)%yfp(nxi)= ap(pp)%yfp(nxi)+ dNshpfun(jj,1)*ap(pp)%yll(cone(jj))
				 ap(pp)%zfp(nxi)= ap(pp)%zfp(nxi)+ dNshpfun(jj,1)*ap(pp)%zll(cone(jj))
			enddo
		enddo !gp
	enddo !nel	

	!print*, "------------ ap(pp)%xfp,ap(pp)%yfp, ap(pp)%zfp ------------"	
	!do ii = 1,nl
	!	print'(3E12.4)', ap(pp)%xfp(ii)-0.025,ap(pp)%yfp(ii)-2.d0,ap(pp)%zfp(ii)-4.d0
	!enddo
	
	endif !mslv
enddo !pp
end subroutine ap_ugp



!-----subroutine to extrapolate forces from xip to CPs
subroutine ap_fgp
use mod_param_fibm,  only: nxie, nl, nxi_vec, nel, nderiv, &
											 p1_fordr, nno
use mod_common_fibm, only: ap, pmax, gp, ispan, uKnot, CONS, &
											 dNshpfun, uu, icp, jcp, cone

IMPLICIT NONE

INTEGER              :: pp, ii, jj, nxi, info, IPIV(nno)
!REAL                :: !state_x0(p1_fordr,3)
REAL*8               :: Amat(nno, nno), WORK(nno), &
												Bxmat(nno),Bymat(nno),Bzmat(nno)
external                SGETRF, SGETRI


do pp= 1,pmax
  if (ap(pp)%mslv .gt. 0) then

	! initializing forces
	ap(pp)%fxll(1:nno) = 0.0!D+00
	ap(pp)%fyll(1:nno) = 0.0!D+00
	ap(pp)%fzll(1:nno) = 0.0!D+00
	Amat               = 0.0!D+00
	Bxmat              = 0.0
	Bymat              = 0.0
	Bzmat              = 0.0

	do ii = 1,nel
	  cone = CONS(ii,1:p1_fordr)

		do gp = 1,nxie   
			nxi = (ii-1)*nxie+gp
			uu  = nxi_vec(nxi) 

		  ispan = FindSpan(nno, p1_fordr-1, uu, uKnot)
			call DersBasisFuns(ispan, uu, p1_fordr-1, nderiv-1, & 
		                                              uKnot, dNshpfun)
			! compute A and b arrays
			do icp = 1,p1_fordr
			   do jcp = 1,p1_fordr
						Amat(cone(icp),cone(jcp))=Amat(cone(icp),cone(jcp))+dNshpfun(icp,1)*dNshpfun(jcp,1)
				 enddo
				 Bxmat(cone(icp)) = Bxmat(cone(icp))+ ap(pp)%fxl(nxi)*dNshpfun(icp,1)
				 Bymat(cone(icp)) = Bymat(cone(icp))+ ap(pp)%fyl(nxi)*dNshpfun(icp,1)
				 Bzmat(cone(icp)) = Bzmat(cone(icp))+ ap(pp)%fzl(nxi)*dNshpfun(icp,1)
			enddo

		enddo !gp	
	enddo !nel	
	
	!print*,"------------ Amat raw --------------"
  !do ii = 1,nno
	!	print '(7E12.4)', (Amat(ii,jcp), jcp=1,nno)
	!enddo

	call DGETRF(nno,nno,Amat,nno,IPIV,info)
	
	!print*,"------------ Amat after TRF --------------"
  !do ii = 1,nno
	!	print '(7E12.4)', (Amat(ii,jcp), jcp=1,nno)
	!enddo
  
	call DGETRI(nno,Amat,nno,IPIV,WORK,nno,info)

	!print*,"------------ Amat after inverse --------------"
  !do ii = 1,nno
	!	print '(7E12.4)', (Amat(ii,jcp), jcp=1,nno)
	!enddo
	!
	!print*,"------------ Bx By Bz mat --------------"
  !do ii = 1,nno
	!	print '(3E12.4)', Bxmat(ii), Bymat(ii), Bzmat(ii) 
	!enddo

	ap(pp)%fxll= matmul(Amat,Bxmat)
	ap(pp)%fyll= matmul(Amat,Bymat)
	ap(pp)%fzll= matmul(Amat,Bzmat)

	!print*, "------------ ap(pp)%fxll,ap(pp)%fyll, ap(pp)%fzll ------------"	
	!do ii = 1,nno
	!	print'(3E12.4)', ap(pp)%fxll(ii),ap(pp)%fyll(ii),ap(pp)%fzll(ii)
	!enddo
	

	endif !mslv
enddo !pp
endsubroutine ap_fgp


END MODULE mod_initialSetup
