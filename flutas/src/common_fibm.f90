module mod_common_fibm
use mod_param_fibm,		only: i1,j1,k1, pointmax,imax,jmax,kmax,nzero,itot,jtot,ktot,&
												nl, diametr, PI, nel, p1_fordr, ngps, ngpb, &
												nderiv, nno, ndof, nir, npmax, np, nl, nl2, &
												nl3, nl4, Eflag, np, c_type
!
!!real ,dimension(0:i1,0:j1,0:k1) :: unew, vnew, wnew, uo, vo, wo, uoo, voo, woo, pnew, po, dudt, dvdt, dwdt, &
!!                                   dudtold, dvdtold, dwdtold, &
real ,dimension(0:i1,0:j1,0:k1) :: forcex, forcey, forcez, forcexold, forceyold, forcezold, &
                                   uac, vac, wac, visc, ufft, vfft, wfft, vf
real ,dimension(1:pointmax,0:(imax/nzero)+1,0:(jmax/nzero)+1,0:(kmax/nzero)+1) :: cpx,cpy,cpz,cpnx,cpny,cpnz,cpa,&
 cpxt,cpyt,cpzt,cpdxdt,cpdydt,cpdzdt,cpfxl,cpfyl,cpfzl !v in suspension and collision
integer ,dimension(1:pointmax,0:(imax/nzero)+1,0:(jmax/nzero)+1,0:(kmax/nzero)+1) :: pid,pidp,pidl,pidrecv !v in collision
!!real(8) :: time,dt   !v not required
!!real:: dtold         !v not required
!!integer:: begin      !v not required
real(8) :: wi(itot+15), wj(jtot+15), wk(ktot+15)          !v solver and initsolver
real(8) :: wi2(2*itot+15), wj2(2*jtot+15), wk2(2*ktot+15) !v solver and initsolver
real, dimension(imax,jmax) :: xyrt                        !v solver and initsolver
real, dimension(kmax) :: a,b,c
real, dimension(kmax) :: ak,bk,ck                         !v solver and initsolver
real, dimension(jtot) :: aj,bj,cj
real :: forcextot,forceytot,forceztot                     !v forcing and interp_spread
!!real :: u_bulk,v_bulk,w_bulk,betafilament               !v only v_bulk is used in main
real :: betafilament                                    !v only v_bulk is used in main
real :: wallshearold,wallshearnew                          !v in rk3
real :: dpdx_sumrk,dudyt,dudyb,epsiloncol               !v comp_grad and forcing 
!!integer :: rkiter,istep,l,pnl,fcounter,lcounter         !v rkiter is used in many files 
integer :: l,pnl,fcounter,lcounter           !v rkiter is used in many files 
!real :: rkparalpha, length(npmax), !tensionplot(NL)         !v rkparalpha, tensionplot not used
real, dimension(np) :: theta_Spagetti,phi_Spagetti
real::s11l,s11g,s22l,s22g,s33l,s33g,s23l,s23g             !v in suspension
!****************************************tension variables for filament*******************************************
!
! particles
!
! for simplifying the communication between threads when there is a
! new master, the derived type 'particle' should be organized in this way:
!
!   (1st) real data that has to be communicated when there is a new master
!         is defined continuously
!
!   (2nd) integer data that has to be communicated when there is a new master
!         is defined contiguously
!
!   (3rd) real data that does not have to be communicated
!
!   (4th) integer data that does not have to be communicated


!
type particle
 real :: x,y,z,integralx, integraly,integralz
 real, dimension(nno) :: xll, yll, zll, qn1, qn2, qn3, qn4, & 
                        xfpo,yfpo,zfpo, qno1, qno2, qno3, qno4, &
                        xfpold,yfpold,zfpold, qnold1, qnold2, qnold3, qnold4, & 
                        dxdtl,dydtl,dzdtl, omg1, omg2, omg3,  & 
                        dxdto,dydto,dzdto, omgo1, omgo2, omgo3,  & 
                        ua,va,wa, omgd1, omgd2, omgd3,   &        
                        uao,vao,wao, omgdo1, omgdo2, omgdo3, &
                        dxds,dyds,dzds, &                        
                        fxll,fyll,fzll, mxl, myl, mzl,  &         
                        ull,vll,wll, omgl1, omgl2, omgl3,  &      
                        fcx,fcy,fcz, &              
                        xt1,yt1,zt1, &             
                        xt2,yt2,zt2, &
                        xt3,yt3,zt3, &
                        xt4,yt4,zt4, &
                        xtc,ytc,ztc, &             
                        fx1,fy1,fz1, &            
                        fx2,fy2,fz2, &
                        fx3,fy3,fz3, &
                        fx4,fy4,fz4!, &
                        !cp                        
 ! total ammount of reals to be communicated: 6+12*nl
 integer :: mslv
 integer, dimension(8) :: nb   
 ! quantities at parameteric points
 real, dimension (nl)  :: xfp, yfp, zfp, &
													dxdt, dydt, dzdt, &
													ul, vl, wl, &
													fxl, fyl, fzl
end type particle
type(particle), dimension(npmax) :: ap ! 'a particle' array
type(particle), dimension(npmax) :: sp ! send particle array (re-ordering of masters)      


type pd
  real, dimension(1:NL) :: xfp,yfp,zfp
end type pd
type(pd), dimension(1:npmax) :: p_data  !v in comp_tension
!
! The structure of 'particle_sumrk' follows the same criterion as 'particle'
! for its organization
!

type particle_interior
  real :: x,y,z
  real, dimension(nl+nl2+nl3+nl4) :: xfp,yfp,zfp,dvlagr
end type particle_interior
!
real :: dVlagr,dVeul
real :: radfp
real, dimension(1:NL) :: thetarc,phirc ! angular position of the Lfps
!
integer :: pmax, npmstr
!
!v-----------------------------------------------------------------------------------------------
! knot insertion values
REAL(8), ALLOCATABLE   :: uKnot(:), coefs(:,:), kntins(:), coefs_np(:,:,:), coefs_temp(:,:)
! Dummy / Temp Arrays
REAL(8), ALLOCATABLE   :: vec_dummy(:), uKnotp(:), coefsp(:,:), coefsp_temp(:,:), &
													coefsp_np(:,:,:), uKnoth(:), coefsh(:,:), coefsh_temp(:,:), &
													coefsh_np(:,:,:), uKnotad(:), coefsad(:,:), coefsad_temp(:,:)!, coefsad_np(:,:,:)

! Position etc related
REAL(16)               :: Arf, I_zz, I_yy, J_xx

! Connectivity related 
INTEGER                :: CONS(1:nel, 1:p1_fordr), cone(p1_fordr), &
                          ie_con(6*p1_fordr)
REAL                   :: unqKnot(1:nel+1), elRangeU(1:nel,2), &
                          elU(2), c1, c2
! Initial Setup array                             
REAL(16)               :: J_s(nel*ngps,np), J_b(nel*ngpb,np), &
													Jp_s(nel*ngps), Jp_b(nel*ngpb)
REAL(16)               :: NdN_array(nel*ngps*p1_fordr, nderiv), &
                          NdN_arrayb(nel*ngpb*p1_fordr, nderiv), &
                          mInerCP(6,6),mInerCP2(6,6), &
													NdN_temp(ngps*p1_fordr, nderiv), &
													NdNb_temp(ngpb*p1_fordr, nderiv)!, &
REAL(8) 							 ::	strain_k0(nel*ngps,3,np)

! Dummy / Temp Arrays
INTEGER                :: icp, jcp, iroco(3), vec3(3), vec6(6)
! Quadrature-rule related
REAL(16)               :: xigs(ngps), wgs(ngps), xigb(ngpb), wgb(ngpb)
! Position etc related 
REAL(16)               :: xn(nno,7), xn_old(nno,7), Xi(nno,7,np) !Xi(nno, 7), 
INTEGER, PARAMETER     :: all_dof(ndof) = (/ (ii, ii=1,ndof) /)
!-------INPUT-------
INTEGER                :: ir_dof(nir)
INTEGER, PARAMETER     :: es = 10000
REAL(16)               :: eN_mat(20), &
                          xn_storeG(nno, 7, es+1), Fexl(es), Mexl(es), &
                          Fexl0(es), Mexl0(es), Fext1(2), Fext2(es-2)

          
REAL(8)                :: dNshpfun(p1_fordr,nderiv)
INTEGER(4)             :: gp, ispan
REAL                   :: eta, Nshpfun(p1_fordr) 
REAL(8)                :: uu, dt_step, eN!, dt

end module mod_common_fibm
