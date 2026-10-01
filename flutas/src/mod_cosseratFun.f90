MODULE mod_cosseratFun
use bspline
CONTAINS


SUBROUTINE stateIncrement(state_x, state_delta, p1_fordr, state_x_inc)
IMPLICIT NONE

INTEGER                 :: i
INTEGER, PARAMETER      :: rk = kind (1.0D+00) 
INTEGER, INTENT(IN)     :: p1_fordr
REAL(16), INTENT(INOUT) :: state_x_inc(p1_fordr,7)
REAL(16), INTENT(IN)    :: state_x(p1_fordr,7), state_delta(6*p1_fordr)
REAL(16)                :: rel_q(4), qN_old(4), rel_q_temp(3)

do i = 1, p1_fordr
    ! increment to pos 'xyz' coordinates
    state_x_inc(i,1) = state_x(i,1) + state_delta((i-1)*6 + 1)
    state_x_inc(i,2) = state_x(i,2) + state_delta((i-1)*6 + 2)
    state_x_inc(i,3) = state_x(i,3) + state_delta((i-1)*6 + 3)
    
    ! increment to rot coordinates in terms of quaternions
    qN_old              = state_x(i, 4:7)
    rel_q_temp          = state_delta((i-1)*6 + 4 : i*6)
    call q_from_Rotv(rel_q_temp, rel_q)
    call quad_prod(qN_old, rel_q, state_x_inc(i, 4:7))
end do
END SUBROUTINE stateIncrement



SUBROUTINE q_from_Rotv(angle_axis, q)
IMPLICIT NONE

INTEGER, PARAMETER      :: rk = kind (1.0D+00) 
REAL(16), INTENT(IN)    :: angle_axis(3) 
REAL(16), INTENT(INOUT) :: q(4) 
REAL(16)                :: theta, theta_squared, half_theta, k

theta_squared = angle_axis(1)*angle_axis(1) + &
    angle_axis(2)*angle_axis(2) + angle_axis(3)*angle_axis(3)

if (theta_squared > 1e-30)  then  ! non-zero rotation
    theta      = sqrt(theta_squared)
    half_theta = theta/2
    k          = sin(half_theta) / theta
    q(1)       = cos(half_theta)
    q(2)       = angle_axis(1) * k
    q(3)       = angle_axis(2) * k
    q(4)       = angle_axis(3) * k
else
    ! For almost zero rotation
    q(1)       = 1.0  !as cos(0) = 1
    q(2)       = 0.0  !angle_axis(1) * k
    q(3)       = 0.0  !angle_axis(2) * k
    q(4)       = 0.0  !angle_axis(3) * k
end if
END SUBROUTINE q_from_Rotv 


SUBROUTINE quad_prod(qa, qb, qp)
IMPLICIT NONE

INTEGER, PARAMETER      :: rk = kind (1.0D+00) 
REAL(16), INTENT(IN)    :: qa(4), qb(4)
REAL(16), INTENT(INOUT) :: qp(4) 

qp(1) = qa(1)*qb(1) - qa(2)*qb(2) - qa(3)*qb(3) - qa(4)*qb(4)
qp(2) = qa(1)*qb(2) + qa(2)*qb(1) + qa(3)*qb(4) - qa(4)*qb(3)
qp(3) = qa(1)*qb(3) - qa(2)*qb(4) + qa(3)*qb(1) + qa(4)*qb(2)
qp(4) = qa(1)*qb(4) + qa(2)*qb(3) - qa(3)*qb(2) + qa(4)*qb(1) 

!if (qp(1) .LE. 0) then
!    qp(1) = -1*qp(1);
!    qp(2) = -1*qp(2);
!    qp(3) = -1*qp(3);
!    qp(4) = -1*qp(4);
!    print*, 'stopped inside quad_prod as -ve qN is obtained'
!    STOP
!endif

END SUBROUTINE quad_prod



SUBROUTINE getConjugate(q, q_astrk)
IMPLICIT NONE

INTEGER, PARAMETER      :: rk = kind (1.0D+00) 
REAL(16), INTENT(IN)    :: q(4)
REAL(16), INTENT(INOUT) :: q_astrk(4)

q_astrk = (/ q(1), -q(2), -q(3), -q(4) /)

END SUBROUTINE getConjugate


!INTEGER, PARAMETER      :: rk = kind(1.0D+00) 
!REAL(rk), INTENT(IN)    :: angle_axis(3) 
!REAL(rk), INTENT(INOUT) :: q(4) 




SUBROUTINE q_to_AngAxis(q, a_angle, a_axis)
IMPLICIT NONE

INTEGER, PARAMETER      :: rk = kind(1.0D+00) 
REAL(16), INTENT(IN)    :: q(4)
REAL(16), INTENT(INOUT) :: a_angle, a_axis(3)
REAL(16)                :: sin_squared, sin_theta, k 

sin_squared = q(2)*q(2) + q(3)*q(3) + q(4)*q(4) !mag. of rot axis vector
if (sin_squared .GT. 1.0D-30) then
	sin_theta  = sqrt(sin_squared)
    	a_angle    = 2.0 * atan2 (sin_theta, q(1))
    	k          = 1.0 / sin_theta
    	a_axis(1)  = q(2) * k
    	a_axis(2)  = q(3) * k
    	a_axis(3)  = q(4) * k
        a_axis     = a_axis / norm2(a_axis)
else
	a_angle   = 0.0d00
    	a_axis(1) = 1.0d00
    	a_axis(2) = 0.0d00
    	a_axis(3) = 0.0d00
end if
END SUBROUTINE q_to_AngAxis



SUBROUTINE Rotatef(q, A, A_rotated)
IMPLICIT NONE

INTEGER, PARAMETER     :: rk=kind(1.0D+00) !no. of derv + shpfun
REAL(16), INTENT(IN)    :: q(4), A(3)
REAL(16), INTENT(INOUT) :: A_rotated(3) 
REAL(16)                :: q0q0, q1q1, q2q2, q3q3, q0q1, q0q2, q0q3, &
			  q1q2, q1q3, q2q3

q0q0 = q(1)*q(1)
q1q1 = q(2)*q(2)
q2q2 = q(3)*q(3)
q3q3 = q(4)*q(4)
q0q1 = q(1)*q(2)
q0q2 = q(1)*q(3)
q0q3 = q(1)*q(4)
q1q2 = q(2)*q(3)
q1q3 = q(2)*q(4)
q2q3 = q(3)*q(4)

A_rotated(1) = ((q0q0 + q1q1) * 2. - 1.) * A(1) + &
        	((q1q2 - q0q3) * 2.) * A(2) + &
		((q1q3 + q0q2) * 2.) * A(3)
A_rotated(2) = ((q1q2 + q0q3) * 2.) * A(1) + &
		((q0q0 + q2q2) * 2. - 1.) * A(2) + &
		((q2q3 - q0q1) * 2.) * A(3)
A_rotated(3) = ((q1q3 - q0q2) * 2.) * A(1) + &
		((q2q3 + q0q1) * 2.) * A(2) + &
		((q0q0 + q3q3) * 2. - 1.) * A(3)
END SUBROUTINE Rotatef



SUBROUTINE RotateBack(q, A, A_rotated)
IMPLICIT NONE

INTEGER, PARAMETER      :: rk=kind(1.0D+00) !no. of derv + shpfun
REAL(16), INTENT(IN)    :: q(4), A(3)
REAL(16), INTENT(INOUT) :: A_rotated(3) 
REAL(16)                :: q0q0, q1q1, q2q2, q3q3, q0q1, q0q2, q0q3, &
			  q1q2, q1q3, q2q3

q0q0 = +q(1)*q(1)
q1q1 = +q(2)*q(2)
q2q2 = +q(3)*q(3)
q3q3 = +q(4)*q(4)
q0q1 = -q(1)*q(2)
q0q2 = -q(1)*q(3)
q0q3 = -q(1)*q(4)
q1q2 = +q(2)*q(3)
q1q3 = +q(2)*q(4)
q2q3 = +q(3)*q(4)

A_rotated(1) = ((q0q0 + q1q1) * 2. - 1.) * A(1) + & 
    ((q1q2 - q0q3) * 2.) * A(2) + ((q1q3 + q0q2) * 2.) * A(3)
A_rotated(2) = ((q1q2 + q0q3) * 2.) * A(1) + &
    ((q0q0 + q2q2) * 2. - 1.) * A(2) + ((q2q3 - q0q1) * 2.) * A(3)
A_rotated(3) = ((q1q3 - q0q2) * 2.) * A(1) + &
    ((q2q3 + q0q1) * 2.) * A(2) + ((q0q0 + q3q3) * 2. - 1.) * A(3)

END SUBROUTINE RotateBack


SUBROUTINE rotMat_from_q(qR, R)
IMPLICIT NONE

INTEGER, PARAMETER     :: rk=kind(1.0D+00) !no. of derv + shpfun
REAL(16), INTENT(IN)    :: qR(4)
REAL(16), INTENT(INOUT) :: R(3,3)
REAL(16)                :: q0q0, q1q1, q2q2, q3q3, q0q1, q0q2, q0q3, &
			  q1q2, q1q3, q2q3

q0q0 = qR(1)*qR(1)
q1q1 = qR(2)*qR(2)
q2q2 = qR(3)*qR(3)
q3q3 = qR(4)*qR(4)
q0q1 = qR(1)*qR(2)
q0q2 = qR(1)*qR(3)
q0q3 = qR(1)*qR(4)
q1q2 = qR(2)*qR(3)
q1q3 = qR(2)*qR(4)
q2q3 = qR(3)*qR(4)

R(1,1) = (q0q0 + q1q1) * 2. - 1.
R(1,2) = (q1q2 - q0q3) * 2.
R(1,3) = (q1q3 + q0q2) * 2.
R(2,1) = (q1q2 + q0q3) * 2.
R(2,2) = (q0q0 + q2q2) * 2. - 1.
R(2,3) = (q2q3 - q0q1) * 2.
R(3,1) = (q1q3 - q0q2) * 2.
R(3,2) = (q2q3 + q0q1) * 2.
R(3,3) = (q0q0 + q3q3) * 2. - 1.

END SUBROUTINE rotMat_from_q 




SUBROUTINE vCross(va, vb, v_cross)
IMPLICIT NONE

INTEGER, PARAMETER     :: rk=kind(1.0D+00) !no. of derv + shpfun
REAL(16), INTENT(IN)    :: va(3), vb(3)
REAL(16), INTENT(INOUT) :: v_cross(3) 

v_cross(1) = va(2)*vb(3) - va(3)*vb(2)
v_cross(2) = va(3)*vb(1) - va(1)*vb(3)
v_cross(3) = va(1)*vb(2) - va(2)*vb(1)

END SUBROUTINE



SUBROUTINE compute_fdise(fd_CP, state_x, elU, nderiv, ngpb, xigb, & 
              wgb, NdNe_array, p1_fordr, nno, uKnot, Arf, fdise)
IMPLICIT NONE

INTEGER, PARAMETER     :: rk=kind(1.0D+00) !no. of derv + shpfun
INTEGER, INTENT(IN)    :: nderiv, ngpb, p1_fordr, nno
REAL(8), INTENT(IN)    :: uKnot(nno+p1_fordr)

REAL(16), INTENT(IN)   :: state_x(p1_fordr,3), xigb(ngpb), wgb(ngpb), &
                          NdNe_array(ngpb*p1_fordr, nderiv), &
													fd_CP(6*p1_fordr), Arf
REAL, INTENT(IN)       :: elU(2) 

REAL(16),INTENT(INOUT) :: fdise(6*p1_fordr)
REAL                   :: eta, c1, c2, state_x0(p1_fordr,3), &
                          N(p1_fordr)
REAL(8)                :: dNshpfun(p1_fordr,nderiv), uu, dr(3), &
                          Jue, Jsu, mfactor, wg
INTEGER                :: i, gp, icp, ispan, inx

c1 = 0.5 * (elU(2) - elU(1))
c2 = 0.5 * (elU(2) + elU(1))

do gp = 1,ngpb
	eta = xigb(gp)
	wg  = wgb(gp)
	uu  = c1*eta + c2
	Jue = c1

	!ispan = FindSpan(nno, p1_fordr-1, uu, uKnot)
	!call DersBasisFuns(ispan, uu, p1_fordr-1, nderiv-1, & 
  !                                            uKnot, dNshpfun)

	inx = (gp-1)*p1_fordr
  do icp=1,p1_fordr
    dNshpfun(icp,:)= NdNe_array(inx+icp,:)
  enddo
 
	dr = 0.0 !(/ 0. 0. 0. /)     
  do icp = 1,p1_fordr
  		dr = dr + state_x(icp,1:3)*dNshpfun(icp,2)
  end do
  Jsu = SQRT( dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3) )
	
	mfactor = Jue * Jsu * wg * Arf
	do icp = 1, p1_fordr
		  fdise((icp-1)*6+1 : icp*6) = fdise((icp-1)*6+1 : icp*6) + &
							dNshpfun(icp,1) * fd_CP(6*(icp-1)+1: 6*icp) * mfactor 
  end do


end do ! quadrature loop

END SUBROUTINE compute_fdise



SUBROUTINE compute_fbody(g_vec, rhof, rho_sp, state_x, elU, nderiv, ngpb, xigb, & 
              wgb, NdNe_array, p1_fordr, nno, uKnot, Arf, fbody)
IMPLICIT NONE

INTEGER, PARAMETER     :: rk=kind(1.0D+00) !no. of derv + shpfun
INTEGER, INTENT(IN)    :: nderiv, ngpb, p1_fordr, nno
REAL(8), INTENT(IN)    :: uKnot(nno+p1_fordr), rhof, rho_sp

REAL(16), INTENT(IN)   :: state_x(p1_fordr,3), xigb(ngpb), wgb(ngpb), &
                          NdNe_array(ngpb*p1_fordr, nderiv), &
													g_vec(3), Arf 
REAL, INTENT(IN)       :: elU(2) 

REAL(16),INTENT(INOUT) :: fbody(6*p1_fordr)
REAL                   :: eta, c1, c2, state_x0(p1_fordr,3), &
                          N(p1_fordr)
REAL(8)                :: dNshpfun(p1_fordr,nderiv), uu, dr(3), &
                          Jue, Jsu, mfactor, wg
INTEGER                :: i, gp, icp, ispan, inx

c1 = 0.5 * (elU(2) - elU(1))
c2 = 0.5 * (elU(2) + elU(1))

do gp = 1,ngpb
	eta = xigb(gp)
	wg  = wgb(gp)
	uu  = c1*eta + c2
	Jue = c1

	!ispan = FindSpan(nno, p1_fordr-1, uu, uKnot)
	!call DersBasisFuns(ispan, uu, p1_fordr-1, nderiv-1, & 
  !                                            uKnot, dNshpfun)

	inx = (gp-1)*p1_fordr
  do icp=1,p1_fordr
    dNshpfun(icp,:)= NdNe_array(inx+icp,:)
  enddo
 
	dr = 0.0 !(/ 0. 0. 0. /)     
  do icp = 1,p1_fordr
  		dr = dr + state_x(icp,1:3)*dNshpfun(icp,2)
  end do
  Jsu = SQRT( dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3) )
	
	mfactor = Jue * Jsu * wg * Arf
	do icp = 1, p1_fordr
		  fbody((icp-1)*6+1 : icp*6-3) = fbody((icp-1)*6+1 : icp*6-3) + &
							dNshpfun(icp,1) * (rhof-rho_sp)  * g_vec * mfactor 
  end do


end do ! quadrature loop

END SUBROUTINE compute_fbody



SUBROUTINE compute_finer( mInerCP, elU, state_x, state_w, nderiv, &
              ngpb, xigb, wgb, rhof, I_yy, I_zz, J, uae, NdNe_array, &
	            p1_fordr, nno, uKnot, finert, Eflag, rhoA, rhoI_yy, rhoI_zz)
IMPLICIT NONE

INTEGER, PARAMETER     :: rk=kind(1.0D+00)
INTEGER, INTENT(IN)    :: nderiv, ngpb, p1_fordr, nno
REAL(8), INTENT(IN)    :: uKnot(nno+p1_fordr)
REAL(16), INTENT(IN)   :: state_x(p1_fordr,7), state_w(p1_fordr,6), &
												  xigb(ngpb), wgb(ngpb), mInerCP(6,6), &
                          uae(p1_fordr,6), NdNe_array(ngpb*p1_fordr, nderiv) 
REAL, INTENT(IN)       :: elU(2), rhof
REAL(16), INTENT(IN)   :: I_zz, I_yy, J
LOGICAL, INTENT(IN)    :: Eflag
REAL(16),INTENT(IN)    :: rhoA, rhoI_yy, rhoI_zz
REAL(16),INTENT(INOUT) :: finert(6*p1_fordr)
REAL                   :: eta, c1, c2, state_x0(p1_fordr,3), &
                          N(p1_fordr)
REAL(8)                :: uu, dNshpfun(p1_fordr,nderiv) 
REAL(16)               :: dr(3), &
                          da(3), Jue, Jsu, mfactor, wg, q_i(4), &
                          q_bar(4), q_barc(4), q_delta(4), qda(4), &
		    	                qRot(4), wl_i(3), w_i(3), w_sect(3), & 
 			                    delta_rot_angle, delta_rot_dir(3), &
			                    fcentf(3), fgyro(3), fcg_temp(3), f_vec(6)
INTEGER                :: i, gp, icp, jcp, ispan, jj, vec6(6), inx

c1 = 0.5 * (elU(2) - elU(1))
c2 = 0.5 * (elU(2) + elU(1))

do gp = 1,ngpb
	eta = xigb(gp)
	wg  = wgb(gp)
	uu  = c1*eta + c2
	Jue = c1

	!ispan = FindSpan(nno, p1_fordr-1, uu, uKnot)
	!call DersBasisFuns(ispan, uu, p1_fordr-1, nderiv-1, & 
  !                                            uKnot, dNshpfun)

	inx = (gp-1)*p1_fordr
  do jj=1,p1_fordr
		dNshpfun(jj,:)= NdNe_array(inx+jj,:)
  enddo

 
	dr = 0.0D+00 !(/ 0. 0. 0. /)     
  do icp = 1,p1_fordr
  		 dr = dr + state_x(icp,1:3)*dNshpfun(icp,2)
  end do
  Jsu = SQRT( dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3) )
	
	! compute w in abs sys
	da = 0.0D+00 
	q_bar = state_x(1,4:7)
	do icp = 1,p1_fordr
		q_i = state_x(icp,4:7)
		call getConjugate(q_bar, q_barc)
		call quad_prod(q_barc, q_i, q_delta)
		call q_to_AngAxis(q_delta, delta_rot_angle, &
					delta_rot_dir)
		da  = da + delta_rot_dir * delta_rot_angle * &
					dNshpfun(icp,1)
	end do
 	call q_from_Rotv(da, qda)
	call quad_prod(q_bar, qda, qRot)
	
	! Compute w in abs csys
	w_sect = 0.0D+00
	do icp = 1, p1_fordr
		q_i  = state_x(icp, 4:7)
		wl_i = state_w(icp, 4:6)
		call Rotatef(q_i, wl_i, w_i)
		w_sect = w_sect + w_i * dNshpfun(icp,1)
	end do
	call RotateBack(qRot, w_sect, w_sect)
	
	! Compute quadtaric terms
	fcentf = 0.0D+00
	call Rotatef(qRot, fcentf, fcentf)
	if (Eflag .eqv. .false.) then
		fcg_temp = (/ rhof*(I_zz+I_yy)*w_sect(1), &
			rhof*I_yy*w_sect(2), rhof*I_zz*w_sect(3) /)
	else
		fcg_temp = (/ (rhoI_yy+rhoI_zz)*w_sect(1), &
			rhoI_yy*w_sect(2), rhoI_zz*w_sect(3) /)
	endif
	call vCross(w_sect, fcg_temp, fgyro)

	!print*,'inside compute_finer'
	!print*, "----------------finert", real(finert,4)
	
	!print*, Jue, Jsu, wg
	!print*, uae

	! Compute Inertia force
	mfactor = Jue * Jsu * wg
	do icp = 1, p1_fordr
	     vec6 = (/ (jj, jj=(icp-1)*6+1,icp*6)  /)
	     finert(vec6) = finert(vec6) + &
					mfactor * dNshpfun(icp,1) * (matmul(mInerCP, uae(icp,1:6)))
	        	     
	     ! add quadratic terms corresp. to gyroscopic force
	     f_vec = (/ fcentf, fgyro /)
			 finert(vec6) = finert(vec6) - mfactor * dNshpfun(icp,1) * f_vec
	end do

	!print*,"--finert--i,gp",gp
  !print'(F10.4)', finert
  !print*,"------mfactor", mfactor

end do
! print*,"--finert--i,gp",gp
! print'(F10.4)', finert
! print*,"------mfactor", mfactor
!print*, '', shape(uae)
!do icp = 1, ubound(uae,1)
!	print '(6F10.4)', uae(icp,:)
!enddo
!print*, ''
END SUBROUTINE compute_finer 



SUBROUTINE  compute_massE( mInerCP, elU, state_x, nderiv, ngpb, xigb, &
                         wgb, J_b, NdNe_array,  p1_fordr, nno, &
                         uKnot, mele )
IMPLICIT NONE

INTEGER, PARAMETER     :: rk=kind(1.0D+00)
INTEGER, INTENT(IN)    :: nderiv, ngpb, p1_fordr, nno
REAL(8), INTENT(IN)    :: uKnot(nno+p1_fordr)
REAL(16), INTENT(IN)   :: state_x(p1_fordr,7), &
												  xigb(ngpb), wgb(ngpb), mInerCP(6,6), &
                          NdNe_array(ngpb*p1_fordr, nderiv), J_b(ngpb) 
REAL, INTENT(IN)       :: elU(2)
REAL(16),INTENT(INOUT) :: mele(6*p1_fordr, 6*p1_fordr) 
REAL                   :: eta, c1, c2
REAL(8)                :: uu, dNshpfun(p1_fordr,nderiv) 
REAL(16)               :: Jue, Jsu, mfactor, wg
INTEGER                :: i, gp, icp, jcp, ispan, j, vec6(6), inx

c1 = 0.5 * ( elU(2) - elU(1) )
c2 = 0.5 * ( elU(2) + elU(1) )

do gp = 1,ngpb
	eta = xigb(gp)
	wg  = wgb(gp)
	uu  = c1*eta + c2
	Jue = c1
	Jsu = J_b(gp)
	
	! note: here basisfun can be used instead of DersBasisFuns
	!ispan = FindSpan(nno, p1_fordr-1, uu, uKnot)
	!call DersBasisFuns(ispan, uu, p1_fordr-1, nderiv-1, & 
  !                                            uKnot, dNshpfun)

	inx = (gp-1)*p1_fordr
  do icp=1,p1_fordr
    dNshpfun(icp,:)= NdNe_array(inx+icp,:)
  enddo

	do icp = 1,p1_fordr
		do jcp = 1,p1_fordr
			mfactor = dNshpfun(icp,1)*dNshpfun(jcp,1) * Jue * Jsu * wg
			vec6    = (/ (j, j=(icp-1)*6+1, 6*icp) /)
			mele(vec6, vec6) =  mele(vec6, vec6) + mfactor * mInerCP
		end do
	end do
end do

!print*, mfactor 
!do i = 1, 6*p1_fordr
!	print '(1x, 12E12.2)', (mele(i,jcp), jcp=1,6*p1_fordr)
!enddo
!STOP

END SUBROUTINE compute_massE




SUBROUTINE compute_fint( state_x, state_w, elU, ngps, xigs, wgs, J_s, &
               NdNe_array, p1_fordr, nno, uKnot, I_zz, I_yy, J, Ks_y, &
	             Ks_z, E_mod, G_mod, Arf, nderiv, fint, &
							 Eflag, E_modA, G_modA, EI_yy, EI_zz, GJxx, strain_gp0, istep)
IMPLICIT NONE

INTEGER, PARAMETER     :: rk=kind(1.0D+00)
INTEGER, INTENT(IN)    :: nderiv, ngps, p1_fordr, nno, istep
REAL(8), INTENT(IN)    :: uKnot(nno+p1_fordr), strain_gp0(ngps,3)
REAL(16), INTENT(IN)   :: state_x(p1_fordr,7), state_w(p1_fordr,6), &
												  xigs(ngps), wgs(ngps), J_s(ngps), &
                          NdNe_array(ngps*p1_fordr, nderiv)
REAL, INTENT(IN)       :: elU(2)
REAL(16),INTENT(IN)    :: E_mod, G_mod 
REAL(16), INTENT(IN)   :: I_zz, I_yy, J, Arf, Ks_y, Ks_z
LOGICAL, INTENT(IN)    :: Eflag
REAL(16),INTENT(IN)    :: E_modA, G_modA, EI_yy, EI_zz, GJxx 
REAL(16),INTENT(INOUT) :: fint(6*p1_fordr)
REAL(kind=8)           :: uu, dNshpfun(p1_fordr,nderiv) 
REAL(8)                :: eta, c1, c2 
REAL(16)               :: dr(3), &
                          da(3), drdt(3), dadt(3), Jue, Jsu, &
								  			  mfactor, wg, q_i(4), &
                          q_bar(4), q_barc(4), q_delta(4), qda(4), &
												  qRot(4), wl_i(3), w_i(3), w_sect(3), & 
													delta_rot_angle, delta_rot_dir(3), &
						 							strain_e0_gp(3), astrain_e_dt(3), &
						 							strain_k0_gp(3), astrain_k_dt(3), &
						 							astrain_e(3), R_trasp(3,3), Rrot(3,3), &
						 							astrain_k(3), astress_n(3), astress_m(3), &
						 							stress_n_abs(3), stress_m_abs(3), & 
						 							Force_i(3), Torque_m(3), Torque_n(3), &
                          Torque_i(3)
INTEGER                :: i, gp, icp, jcp, ispan, inx

c1 = 0.5 * (elU(2) - elU(1))
c2 = 0.5 * (elU(2) + elU(1))


do gp = 1,ngps
	eta = xigs(gp)
	wg  = wgs(gp)
	uu  = c1*eta + c2
	Jue = c1
	Jsu = J_s(gp)

	!ispan = FindSpan(nno, p1_fordr-1, uu, uKnot)
	!call DersBasisFuns(ispan, uu, p1_fordr-1, nderiv-1, & 
  !                       uKnot, dNshpfun)
	
	inx = (gp-1)*p1_fordr
  do icp=1,p1_fordr
    dNshpfun(icp,:)= NdNe_array(inx+icp,:)
  enddo


	! get dNshpfun from NdNe_array(later)
 	da = 0.0 
	q_bar = state_x(1, 4:7)
	do icp = 1,p1_fordr
		q_i = state_x(icp,4:7)
		call getConjugate(q_bar, q_barc)
		call quad_prod(q_barc, q_i, q_delta)
		call q_to_AngAxis(q_delta, delta_rot_angle, &
				delta_rot_dir)
		da  = da + delta_rot_dir * delta_rot_angle * &
				dNshpfun(icp,1)
	end do
 	call q_from_Rotv(da, qda)
	call quad_prod(state_x(1,4:7), qda, qRot)
	call rotMat_from_q(qRot, Rrot)

 	dr = 0.0
	do icp = 1,p1_fordr
		dr = dr + dNshpfun(icp,2) * state_x(icp,1:3)
	end do
	dr = dr / Jsu

	
  !if (istep .ge. 11) then	
	!	 !print*,"----------------dr, gp", real(dr,4), real(gp,4)
	!	 !print*,"----------------strain_gp0, gp", real(strain_gp0(gp,1:3),4), real(gp,4)
	!	 do icp=1,p1_fordr
	!			print*,"---------------- state_x(icp,4:7)", real(state_x(icp,4:7),4)
	!	 enddo
	!	 print*,"---------------- da", real(da,4)
	!endif

	drdt = 0.0
	do icp = 1, p1_fordr
		drdt = drdt + dNshpfun(icp,2) * state_x(icp,1:3)
	end do
	drdt = drdt / Jsu

	! Compute abs spline rotation
	da = 0.0
	do icp = 1,p1_fordr
		q_i = state_x(icp,4:7)
		call getConjugate(qRot, q_barc)
    call quad_prod(q_barc, q_i, q_delta)
    call q_to_AngAxis(q_delta, delta_rot_angle, &
                      delta_rot_dir)
    da  = da + delta_rot_dir * delta_rot_angle * &
                     dNshpfun(icp,2)
	end do
	da = da / Jsu

  !if (istep .ge. 11) then	
	!   !print*,"----------------da, gp", real(da,4), real(gp,4)
	!endif


	! Compute abs. rate of spline rotation
	dadt = 0.0
	do icp = 1,p1_fordr
		q_i  = state_x(icp, 4:7)
		wl_i = state_w(icp, 4:6)
		call Rotatef(q_i, wl_i, w_i)
		dadt = dadt + dNshpfun(icp,2) * w_i
	end do
	dadt = dadt / Jsu

	! Compute local translational strain
	R_trasp      = TRANSPOSE(Rrot)	
	strain_e0_gp = 0.0
	!astrain_e    = matmul(R_trasp, dr) - (/ 0.0, 0.995004, 0.0998334 /) - strain_e0_gp	
	!astrain_e    = matmul(R_trasp, dr) - (/ 0.0, 0.9950, 0.0995 /) - strain_e0_gp	
	!astrain_e    = matmul(R_trasp, dr) - (/ 0.0, 0.9928, 0.1197 /) - strain_e0_gp	
	!astrain_e    = matmul(R_trasp, dr) - (/ 0.0, 0.9511, 0.3090 /) - strain_e0_gp	
	!astrain_e    = matmul(R_trasp, dr) - (/ 0.0, 0.9995, -0.03140 /) - strain_e0_gp	
	!astrain_e    = matmul(R_trasp, dr) - (/ 0.0, 1.0, 0.0 /) - strain_e0_gp	! vec_ex
	!astrain_e    = matmul(R_trasp, dr) - (/ 0.0, 0.0, -1.0 /) - strain_e0_gp	! vec_ex
	!astrain_e    = matmul(R_trasp, dr) - (/ 0.0, 0.9877, 0.1564 /) - strain_e0_gp	! vec_ex
	astrain_e    = matmul(R_trasp, dr) - strain_gp0(gp,1:3) - strain_e0_gp	! vec_ex
	astrain_e_dt = matmul(R_trasp, drdt)
	!print*, strain_gp0(gp,1:3)

  !if (istep .ge. 11) then	
	!	print*, "------------R", real(R_trasp,4)
	!	!print*, "-------astrain_e, gp", real(astrain_e,4), real(gp,4)
	!	!!print*, "-------strain_e0_gp, gp", real(strain_e0_gp,4), real(gp,4)
	!	!print*, "-------matmul(R_trans,dp), gp", real(matmul(R_trasp, dr),4), real(gp,4)
	!	!print*,"--------Elag quant:   ", real(E_modA,4),  real(G_modA,4), GJxx,EI_yy
	!	!print*,"--------Elag wo quant:   ", real(E_mod*Arf,4),  real(G_mod*Arf,4),G_mod*J,E_mod*I_yy
	!endif


 	! Compute local torsional or curvature strain
 	strain_k0_gp = 0.0
	astrain_k    = da - strain_k0_gp	 
	astrain_k_dt = matmul(R_trasp, dadt)

	! Compute local cut force and torque
	if (Eflag .eqv. .false.) then
		astress_n    = 0.0
		astress_n(1) = E_mod * Arf * astrain_e(1)
		astress_n(2) = G_mod * Arf * Ks_y * astrain_e(2)
		astress_n(3) = G_mod * Arf * Ks_z * astrain_e(3)
		astress_m    = 0.0
		astress_m(1) = G_mod * J * astrain_k(1)
		astress_m(2) = E_mod * I_yy * astrain_k(2)
		astress_m(3) = E_mod * I_zz * astrain_k(3)
	else
		astress_n    = 0.0
		astress_n(1) = E_modA * astrain_e(1)
		astress_n(2) = G_modA * Ks_y * astrain_e(2)
		astress_n(3) = G_modA * Ks_z * astrain_e(3)
		astress_m    = 0.0
		astress_m(1) = GJxx  * astrain_k(1)
		astress_m(2) = EI_yy * astrain_k(2)
		astress_m(3) = EI_zz * astrain_k(3)
	end if


	!print*, "-------astres_n, gp", real(astress_n,4), real(gp,4)
	!print*, "-------astres_m, gp", real(astress_m,4), real(gp,4)

	!print*, "here astress_n--->"
        !write(*, '(E20.4)') G_mod * Arf* astrain_e(2)	
	!write(*,'(E20.2)') astrain_e
	!write(*,'(E20.2)') astress_n
	!write(*,'(E10.4)') E_mod, G_mod, Arf

	! Compute internal force and toruque in abs csys
	stress_n_abs = matmul(Rrot, astress_n)
	stress_m_abs = matmul(Rrot, astress_m)
	do icp = 1,p1_fordr
		Force_i  = stress_n_abs * dNshpfun(icp,2) * (wg*Jue)
		q_i      = state_x(icp, 4:7)
		Torque_m = stress_m_abs * dNshpfun(icp,2) * (wg*Jue)
		call vcross(dr, stress_n_abs, Torque_n)
		Torque_n = Torque_n * dNshpfun(icp,1) * (wg*Jue*Jsu)
		Torque_i = Torque_m - Torque_n
		call RotateBack(q_i, Torque_i, Torque_i)
		fint((icp-1)*6 + 1 : icp*6 - 3) = & 
	            fint((icp-1)*6 + 1 : icp*6 - 3) + Force_i
	        fint((icp-1)*6 + 4 : icp*6) = &
         	    fint((icp-1)*6 + 4 : icp*6) + Torque_i
	enddo

	!print*, 'astress_n :----'
	!write(*, *) fint
	!STOP
	!print*, 'strains'
	!do jcp = 1,3!p1_fordr
	!	write(*,'(F10.2)', advance='no') (Rrot(jcp, i), i=1,3)
	!	write(*,*) ''
	!end do
	!if (gp .EQ. 2) then
		!print*, fint
  		!print*, Arf 
		!print*, Torque_i
		!print*, 'gp=', gp
		!STOP
	!endif

end do ! gp loop

!print*, stress_n_abs
!print*, '=============', dNshpfun(icp,2), wg, Jue
!print*, Torque_i
!print*, "here astress_n--->"
!print*, "fint here" 
!write(*,'(E10.4)') astrain_e
!write(*,'(E10.4)')  astrain_e(1) * E_mod* Arf 
!print*, Ks_y 
!print*, 5.  / 6.
!write(*,'(E10.4)') G_mod * Arf * Ks_y * astrain_e(2)
!write(*,'(E10.4)') 1.0D-9 * G_mod* Arf * Ks_y! * astrain_e(2)
!write(*,'(E10.4)') G_mod* Arf * Ks_z * astrain_e(3)
!write(*,'(E10.4)') Torque_i 
 
!do jcp = 1,12
!	write(*,'(E10.4)') fint(jcp)
!enddo
!print*, 'all quantities'
!do i = 1,p1_fordr
!	write(*,'(F10.2)', advance='no') (dNshpfun(i,jcp), jcp=1,nderiv)
!	write(*,*)
!enddo
!print*, 'dNshpfun', dNshpfun(:,2)
!print*, 'state_x', state_x(:,1:3)
!print*, Torque_i

END SUBROUTINE compute_fint


END MODULE mod_cosseratFun
