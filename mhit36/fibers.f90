module fibers
! Single flexible fiber coupled to MHIT36 with the IGA/IBM method of
! "An efficient isogeometric/finite-difference immersed boundary method for the
! fluid-structure interactions of slender flexible structures", CMAME 418 (2024) 116495.
! Ported from FluTAS_IGA (Update_Pos.f90, ext_force.f90, interp_spread.f90,
! initparticles_fibm.f90, mod_setupInitial.f90) with the MPI master/slave logic removed.
!
! Split between host and device:
!  - GPU: interpolation of u,v,w at the Lagrangian points and spreading of the IBM force
!         (only nl points and 3*nl values move between host and device)
!  - CPU: Cosserat rod solver (Newmark + Newton-Raphson), double precision (quad with -DFIB_QUAD)
!
! MHIT36 staggering: p(i,j,k) at ((i-1)dx,(j-1)dx,(k-1)dx), u(i,j,k) at ((i-1.5)dx,(j-1)dx,(k-1)dx),
! v and w shifted by -dx/2 along y and z.
use fib_prec
use fib_param
use bspline
use mod_linspace,    only: quad_rule
use mod_cosseratFun
implicit none
private
public :: fib_init, fib_step, fib_output

! control-point state (same names as the ap(p)% components in FluTAS)
double precision, dimension(nno) :: xll, yll, zll, qn1, qn2, qn3, qn4, &
                                    xfpo, yfpo, zfpo, qno1, qno2, qno3, qno4, &
                                    dxdtl, dydtl, dzdtl, omg1, omg2, omg3, &
                                    dxdto, dydto, dzdto, omgo1, omgo2, omgo3, &
                                    ua, va, wa, omgd1, omgd2, omgd3, &
                                    uao, vao, wao, omgdo1, omgdo2, omgdo3, &
                                    fxll, fyll, fzll, mxl, myl, mzl
! Lagrangian points: positions, fiber velocity, fluid velocity, IBM force
double precision, allocatable, dimension(:) :: xfp, yfp, zfp, ul, vl, wl, fxl, fyl, fzl
double precision, dimension(nl) :: dxdt, dydt, dzdt

! IGA discretisation (constant in time)
real(8)  :: uKnot(nno+p1_fordr)
integer  :: CONS(nel,p1_fordr)
real     :: elRangeU(nel,2)
real(fk) :: xigs(ngps), wgs(ngps), xigb(ngpb), wgb(ngpb)
real(fk) :: J_s(nel*ngps), J_b(nel*ngpb)
real(fk) :: NdN_array(nel*ngps*p1_fordr,nderiv), NdN_arrayb(nel*ngpb*p1_fordr,nderiv)
real(fk) :: mInerCP(6,6), mInerCP2(6,6)
real(8)  :: strain_k0(nel*ngps,3)
real(fk) :: Arf, I_yy, I_zz, J_xx
! basis functions at the Lagrangian points and L2 projection matrix (points -> control points)
double precision :: Nxi(p1_fordr,nl), Amat(nno,nno)
integer          :: exi(nl)

double precision :: rhofl       ! fluid density
integer          :: nr_last     ! Newton iterations used at the last step

contains

!===========================================================================
subroutine fib_init(restart, tstart, rho)
use param, only: lx
implicit none
integer, intent(in)          :: restart, tstart
double precision, intent(in) :: rho
integer  :: i, j, p, ii, gp, icp, jcp, ispan, counter, inx, l
real(8)  :: coefsi(4,2), uKnoti(0:3), kntins(nel-1)
real(8), allocatable :: uKnotp(:), coefsp(:,:), coefs(:,:)
real(8)  :: dNshpfun(p1_fordr,nderiv), uu, dir(3), axis(3), ang
real(fk) :: Xi(nno,3), state_x0(p1_fordr,3), dr0(3), q0(4), R0(3,3), rotv(3)
real     :: elU(2), eta, c1, c2
integer  :: nhp

rhofl = rho

! cross section: area and second moments
if (c_type .eq. "c") then
    Arf  = acos(-1.0_fk) * (0.5_fk*diametr)**2
    I_zz = (acos(-1.0_fk)/4.0_fk) * (0.5_fk*diametr)**4
else
    Arf  = diametr**2
    I_zz = (diametr**4)/12.0_fk
endif
I_yy = I_zz
J_xx = I_zz + I_yy

! straight fiber: linear B-spline between the two end points, degree elevation, knot refinement
dir = dir0/sqrt(sum(dir0**2))
coefsi(1:3,1) = (/ xc0, yc0, zc0 /) - 0.5d0*L_fibr*dir
coefsi(1:3,2) = (/ xc0, yc0, zc0 /) + 0.5d0*L_fibr*dir
coefsi(4,:)   = 1.d0
uKnoti = (/ 0.d0, 0.d0, 1.d0, 1.d0 /)
nhp = 1 + deg_ele
allocate(uKnotp(0:2*(nhp+1)-1), coefsp(4,0:nhp))
call DegreeElevate(4, 1, 1, uKnoti, coefsi, deg_ele, nhp, uKnotp, coefsp)
do i = 1, nel-1
    kntins(i) = dble(i)/dble(nel)
enddo
allocate(coefs(4,0:nno-1))
call RefineKnotVector(4, nhp, nhp, uKnotp, coefsp, nel-2, kntins, uKnot, coefs)
do i = 1, nno
    Xi(i,1:3) = coefs(1:3,i-1)
enddo

! connectivity and element ranges in parametric space
do i = 1, nel
    do j = 1, p1_fordr
        CONS(i,j) = i + j - 1
    enddo
    elRangeU(i,1) = real(i-1)/real(nel)
    elRangeU(i,2) = real(i)/real(nel)
enddo
call quad_rule(ngps, xigs, wgs)
call quad_rule(ngpb, xigb, wgb)

! initial orientation of the material frame
q0 = (/ 1.0_fk, 0.0_fk, 0.0_fk, 0.0_fk /)
if (align_frame) then
    ! rotation taking the material axis 1 onto the fiber direction
    axis = (/ 0.d0, -dir(3), dir(2) /)
    ang  = acos(max(-1.d0, min(1.d0, dir(1))))
    if (sqrt(sum(axis**2)) .gt. 1.d-12) then
        rotv = axis/sqrt(sum(axis**2))*ang
    elseif (dir(1) .lt. 0.d0) then
        rotv = (/ 0.0_fk, acos(-1.0_fk), 0.0_fk /)
    else
        rotv = 0.0_fk
    endif
    call q_from_Rotv(rotv, q0)
endif
call rotMat_from_q(q0, R0)

! basis functions, Jacobians and initial strain at the Gauss points (setupInitial in FluTAS)
counter = 1
do i = 1, nel
    state_x0 = Xi(CONS(i,:),1:3)
    elU = elRangeU(i,1:2)
    c1  = 0.5 * (elU(2) - elU(1))
    c2  = 0.5 * (elU(2) + elU(1))
    do gp = 1, ngps
        eta   = xigs(gp)
        uu    = c1*eta + c2
        ispan = FindSpan(nno, p1_fordr-1, uu, uKnot)
        call DersBasisFuns(ispan, uu, p1_fordr-1, nderiv-1, uKnot, dNshpfun)
        inx = (i-1)*ngps*p1_fordr + (gp-1)*p1_fordr
        do icp = 1, p1_fordr
            NdN_array(inx+icp,:) = dNshpfun(icp,:)
        enddo
        dr0 = 0.0_fk
        do icp = 1, p1_fordr
            dr0 = dr0 + state_x0(icp,1:3)*dNshpfun(icp,2)
        enddo
        J_s(counter) = sqrt(sum(dr0**2))
        ! initial strain in the material frame (R0 = identity reproduces FluTAS)
        strain_k0(counter,1:3) = matmul(transpose(R0), dr0/J_s(counter))
        counter = counter + 1
    enddo
enddo
if (ngpb .ne. ngps) then
    counter = 1
    do i = 1, nel
        state_x0 = Xi(CONS(i,:),1:3)
        elU = elRangeU(i,1:2)
        c1  = 0.5 * (elU(2) - elU(1))
        c2  = 0.5 * (elU(2) + elU(1))
        do gp = 1, ngpb
            eta   = xigb(gp)
            uu    = c1*eta + c2
            ispan = FindSpan(nno, p1_fordr-1, uu, uKnot)
            call DersBasisFuns(ispan, uu, p1_fordr-1, nderiv-1, uKnot, dNshpfun)
            inx = (i-1)*ngpb*p1_fordr + (gp-1)*p1_fordr
            do icp = 1, p1_fordr
                NdN_arrayb(inx+icp,:) = dNshpfun(icp,:)
            enddo
            dr0 = 0.0_fk
            do icp = 1, p1_fordr
                dr0 = dr0 + state_x0(icp,1:3)*dNshpfun(icp,2)
            enddo
            J_b(counter) = sqrt(sum(dr0**2))
            counter = counter + 1
        enddo
    enddo
else
    J_b = J_s
    NdN_arrayb = NdN_array
endif

! mass matrices of the cross section (mInerCP2 includes buoyancy)
mInerCP  = 0.0_fk
mInerCP2 = 0.0_fk
if (.not. Eflag) then
    mInerCP(1,1) = rhof*Arf
    mInerCP(2,2) = rhof*Arf
    mInerCP(3,3) = rhof*Arf
    mInerCP(4,4) = rhof*J_xx
    mInerCP(5,5) = rhof*I_yy
    mInerCP(6,6) = rhof*I_zz
    mInerCP2(1,1) = (rhof - rho)*Arf
    mInerCP2(2,2) = (rhof - rho)*Arf
    mInerCP2(3,3) = (rhof - rho)*Arf
    mInerCP2(4,4) = rhof*J_xx
    mInerCP2(5,5) = rhof*I_yy
    mInerCP2(6,6) = rhof*I_zz
else
    mInerCP(1,1) = rhoA
    mInerCP(2,2) = rhoA
    mInerCP(3,3) = rhoA
    mInerCP(4,4) = rhoI_yy + rhoI_zz
    mInerCP(5,5) = rhoI_yy
    mInerCP(6,6) = rhoI_zz
    mInerCP2(1,1) = rhoA - rho*Arf
    mInerCP2(2,2) = rhoA - rho*Arf
    mInerCP2(3,3) = rhoA - rho*Arf
    mInerCP2(4,4) = rhof*(I_yy + I_zz)
    mInerCP2(5,5) = rhoI_yy
    mInerCP2(6,6) = rhoI_zz
endif

! Lagrangian points: nxie per element, at parametric coordinate (l-1/2)/nl
Amat = 0.d0
do l = 1, nl
    ii = (l-1)/nxie + 1
    exi(l) = ii
    uu = (dble(l) - 0.5d0)/dble(nl)
    ispan = FindSpan(nno, p1_fordr-1, uu, uKnot)
    call DersBasisFuns(ispan, uu, p1_fordr-1, nderiv-1, uKnot, dNshpfun)
    Nxi(:,l) = dNshpfun(:,1)
    do icp = 1, p1_fordr
        do jcp = 1, p1_fordr
            Amat(CONS(ii,icp),CONS(ii,jcp)) = Amat(CONS(ii,icp),CONS(ii,jcp)) + Nxi(icp,l)*Nxi(jcp,l)
        enddo
    enddo
enddo

! state at t=0 (fiber at rest) or from the restart file
if (restart .eq. 0) then
    xll = Xi(:,1)
    yll = Xi(:,2)
    zll = Xi(:,3)
    qn1 = q0(1)
    qn2 = q0(2)
    qn3 = q0(3)
    qn4 = q0(4)
    dxdtl = 0.d0; dydtl = 0.d0; dzdtl = 0.d0
    omg1  = 0.d0; omg2  = 0.d0; omg3  = 0.d0
    ua    = 0.d0; va    = 0.d0; wa    = 0.d0
    omgd1 = 0.d0; omgd2 = 0.d0; omgd3 = 0.d0
    if (pert_amp .ne. 0.d0) call init_mode(dir, q0)
else
    call fib_read_restart(tstart)
endif
xfpo = xll;  yfpo = yll;  zfpo = zll
qno1 = qn1;  qno2 = qn2;  qno3 = qn3;  qno4 = qn4
dxdto = dxdtl; dydto = dydtl; dzdto = dzdtl
omgo1 = omg1;  omgo2 = omg2;  omgo3 = omg3
uao = ua; vao = va; wao = wa
omgdo1 = omgd1; omgdo2 = omgd2; omgdo3 = omgd3
fxll = 0.d0; fyll = 0.d0; fzll = 0.d0
mxl  = 0.d0; myl  = 0.d0; mzl  = 0.d0
nr_last = 0

allocate(xfp(nl), yfp(nl), zfp(nl), ul(nl), vl(nl), wl(nl), fxl(nl), fyl(nl), fzl(nl))
xfp = 0.d0; yfp = 0.d0; zfp = 0.d0
ul  = 0.d0; vl  = 0.d0; wl  = 0.d0
fxl = 0.d0; fyl = 0.d0; fzl = 0.d0
!$acc enter data copyin(xfp,yfp,zfp,ul,vl,wl,fxl,fyl,fzl)

write(*,*) "Fiber: nel, nno, nl       ", nel, nno, nl
write(*,*) "Fiber: length, diameter   ", L_fibr, diametr
write(*,*) "Fiber: solver precision   ", precision(1.0_fk), "digits"
if (lx .lt. L_fibr) write(*,*) "Fiber: WARNING, fiber longer than the box"
if (.not. fib_coupling) write(*,*) "Fiber: WARNING, fluid-fiber coupling is OFF"
end subroutine fib_init



!===========================================================================
! Initial velocity with the shape of free-free bending mode pert_mode of an
! Euler-Bernoulli beam: v(s) = pert_amp*phi(s)/2 along n (|phi(0)| = |phi(L)| = 2),
! with the consistent rotation rate t x dv/ds. Control points are placed at
! their Greville abscissae.
subroutine init_mode(dir, q0)
implicit none
double precision, intent(in) :: dir(3)
real(fk), intent(in)         :: q0(4)
double precision, parameter  :: bLn(4) = (/ 4.730040744862704d0, 7.853204624095838d0, &
                                            10.99560783800167d0, 14.13716549125746d0 /)
double precision :: bL
double precision :: n(3), s, b, sig, phi, dphi, vel(3), wg(3)
real(fk)         :: wgl(3), wloc(3)
integer :: i, p
if (pert_mode .lt. 1 .or. pert_mode .gt. 4) stop "Fiber: pert_mode must be 1, 2, 3 or 4"
bL = bLn(pert_mode)
n = pert_dir - dot_product(pert_dir, dir)*dir
if (sqrt(sum(n**2)) .lt. 1.d-12) stop "Fiber: pert_dir must not be parallel to the fiber"
n = n/sqrt(sum(n**2))
b   = bL/L_fibr
sig = (cosh(bL) - cos(bL))/(sinh(bL) - sin(bL))
p   = p1_fordr - 1
do i = 1, nno
    s    = L_fibr*sum(uKnot(i+1:i+p))/dble(p)
    phi  = cosh(b*s) + cos(b*s) - sig*(sinh(b*s) + sin(b*s))
    dphi = b*(sinh(b*s) - sin(b*s) - sig*(cosh(b*s) + cos(b*s)))
    vel  = 0.5d0*pert_amp*phi*n
    ! rotation rate of the cross section, global frame -> material frame
    wg(1) = dir(2)*n(3) - dir(3)*n(2)
    wg(2) = dir(3)*n(1) - dir(1)*n(3)
    wg(3) = dir(1)*n(2) - dir(2)*n(1)
    wgl   = 0.5d0*pert_amp*dphi*wg
    call RotateBack(q0, wgl, wloc)
    dxdtl(i) = vel(1); dydtl(i) = vel(2); dzdtl(i) = vel(3)
    omg1(i)  = wloc(1); omg2(i) = wloc(2); omg3(i) = wloc(3)
enddo
write(*,*) "Fiber: initial velocity = bending mode", pert_mode, ", end amplitude", pert_amp
end subroutine init_mode



!===========================================================================
! One coupled fluid-fiber step, called on the predicted velocity u* (before the
! Poisson solve). Adds the IBM force to u,v,w and advances the fiber.
subroutine fib_step(istep, dt)
implicit none
integer, intent(in)          :: istep
double precision, intent(in) :: dt
integer :: l

! control points -> Lagrangian points (ap_ugp)
call cp_to_lag

! coupling off: the fiber evolves in vacuum (fxll stays zero) and the fluid is untouched
if (fib_coupling) then
    ! fluid velocity at the Lagrangian points (eulr2lagr), on the GPU
    !$acc update device(xfp,yfp,zfp)
    call fib_interp
    !$acc update self(ul,vl,wl)

    ! IBM force per unit volume at the Lagrangian points: beta*(u_fluid - u_fiber), beta=-1/dt
    do l = 1, nl
        fxl(l) = -(ul(l) - dxdt(l))/dt
        fyl(l) = -(vl(l) - dydt(l))/dt
        fzl(l) = -(wl(l) - dzdt(l))/dt
    enddo

    ! spread the force on the fluid and add it to u* (lagr2eulr + add_ibm_force), on the GPU
    !$acc update device(fxl,fyl,fzl)
    call fib_spread(dt)

    ! Lagrangian points -> control points (ap_fgp); the reaction acts on the fiber
    call lag_to_cp
endif

! advance the fiber (Newmark + Newton-Raphson)
call rod_update(istep, dt)
end subroutine fib_step



!===========================================================================
! positions and velocities of the Lagrangian points from the control points
subroutine cp_to_lag
implicit none
integer :: l, icp, c
do l = 1, nl
    xfp(l)  = 0.d0; yfp(l)  = 0.d0; zfp(l)  = 0.d0
    dxdt(l) = 0.d0; dydt(l) = 0.d0; dzdt(l) = 0.d0
    do icp = 1, p1_fordr
        c = CONS(exi(l),icp)
        xfp(l)  = xfp(l)  + Nxi(icp,l)*xll(c)
        yfp(l)  = yfp(l)  + Nxi(icp,l)*yll(c)
        zfp(l)  = zfp(l)  + Nxi(icp,l)*zll(c)
        dxdt(l) = dxdt(l) + Nxi(icp,l)*dxdtl(c)
        dydt(l) = dydt(l) + Nxi(icp,l)*dydtl(c)
        dzdt(l) = dzdt(l) + Nxi(icp,l)*dzdtl(c)
    enddo
enddo
end subroutine cp_to_lag



!===========================================================================
! L2 projection of the Lagrangian forces on the control points
subroutine lag_to_cp
implicit none
integer :: l, icp, c
double precision :: A(nno,nno), B(nno,3)
B = 0.d0
do l = 1, nl
    do icp = 1, p1_fordr
        c = CONS(exi(l),icp)
        B(c,1) = B(c,1) + fxl(l)*Nxi(icp,l)
        B(c,2) = B(c,2) + fyl(l)*Nxi(icp,l)
        B(c,3) = B(c,3) + fzl(l)*Nxi(icp,l)
    enddo
enddo
A = Amat
call solve_dp(nno, 3, A, B)
! force per unit volume acting on the fiber (FluTAS assumes rho=1)
fxll = rhofl*B(:,1)
fyll = rhofl*B(:,2)
fzll = rhofl*B(:,3)
end subroutine lag_to_cp



!===========================================================================
! Roma et al. (1999) 3-point discrete delta, in grid units
pure function kernel(r) result(phi)
!$acc routine seq
implicit none
double precision, intent(in) :: r
double precision :: phi, ar
ar = abs(r)
if (ar .gt. 1.5d0) then
    phi = 0.d0
elseif (ar .gt. 0.5d0) then
    phi = (5.d0 - 3.d0*ar - sqrt(max(0.d0, 1.d0 - 3.d0*(1.d0-ar)**2)))/6.d0
else
    phi = (1.d0 + sqrt(1.d0 - 3.d0*ar**2))/3.d0
endif
end function kernel



!===========================================================================
! stencil along one direction: s = position/dx, off = grid offset of the variable
! (1 for cell centres, 1.5 for faces). Returns the first (periodic) index and 3 weights.
subroutine stencil(s, off, i0, wgt)
!$acc routine seq
implicit none
double precision, intent(in)  :: s, off
integer, intent(out)          :: i0
double precision, intent(out) :: wgt(3)
integer :: ic, m
ic = nint(s + off)
do m = 1, 3
    wgt(m) = kernel(dble(ic + m - 2) - off - s)
enddo
i0 = ic - 2
end subroutine stencil



!===========================================================================
subroutine fib_interp
use param,    only: nx, dxi
use velocity, only: u, v, w
implicit none
integer :: l, ii, jj, kk, ic, jc, kc, is, js, ks, i, j, k
double precision :: sx, sy, sz, su, sv, sw
double precision :: wxc(3), wyc(3), wzc(3), wxs(3), wys(3), wzs(3), dxil

dxil = dxi
!$acc parallel loop gang vector present(u,v,w,xfp,yfp,zfp,ul,vl,wl) &
!$acc private(wxc,wyc,wzc,wxs,wys,wzs)
do l = 1, nl
    sx = xfp(l)*dxil
    sy = yfp(l)*dxil
    sz = zfp(l)*dxil
    call stencil(sx, 1.0d0, ic, wxc)
    call stencil(sy, 1.0d0, jc, wyc)
    call stencil(sz, 1.0d0, kc, wzc)
    call stencil(sx, 1.5d0, is, wxs)
    call stencil(sy, 1.5d0, js, wys)
    call stencil(sz, 1.5d0, ks, wzs)
    su = 0.d0
    sv = 0.d0
    sw = 0.d0
    do kk = 1, 3
        do jj = 1, 3
            do ii = 1, 3
                ! u: face in x, centre in y and z
                i = modulo(is+ii-1, nx) + 1
                j = modulo(jc+jj-1, nx) + 1
                k = modulo(kc+kk-1, nx) + 1
                su = su + u(i,j,k)*wxs(ii)*wyc(jj)*wzc(kk)
                ! v: face in y
                i = modulo(ic+ii-1, nx) + 1
                j = modulo(js+jj-1, nx) + 1
                sv = sv + v(i,j,k)*wxc(ii)*wys(jj)*wzc(kk)
                ! w: face in z
                j = modulo(jc+jj-1, nx) + 1
                k = modulo(ks+kk-1, nx) + 1
                sw = sw + w(i,j,k)*wxc(ii)*wyc(jj)*wzs(kk)
            enddo
        enddo
    enddo
    ul(l) = su
    vl(l) = sv
    wl(l) = sw
enddo
end subroutine fib_interp



!===========================================================================
! u* = u* + dt * f, with f = sum_l F_l * delta * Arf*ds/dV (as in FluTAS, no 1/rho)
subroutine fib_spread(dt)
use param,    only: nx, dx, dxi
use velocity, only: u, v, w
implicit none
double precision, intent(in) :: dt
integer :: l, ii, jj, kk, ic, jc, kc, is, js, ks, i, j, k
double precision :: sx, sy, sz, cf
double precision :: wxc(3), wyc(3), wzc(3), wxs(3), wys(3), wzs(3), dxil

cf = dt*dble(Arf)*ds/dx**3
dxil = dxi

!$acc parallel loop gang vector present(u,v,w,xfp,yfp,zfp,fxl,fyl,fzl) &
!$acc private(wxc,wyc,wzc,wxs,wys,wzs)
do l = 1, nl
    sx = xfp(l)*dxil
    sy = yfp(l)*dxil
    sz = zfp(l)*dxil
    call stencil(sx, 1.0d0, ic, wxc)
    call stencil(sy, 1.0d0, jc, wyc)
    call stencil(sz, 1.0d0, kc, wzc)
    call stencil(sx, 1.5d0, is, wxs)
    call stencil(sy, 1.5d0, js, wys)
    call stencil(sz, 1.5d0, ks, wzs)
    do kk = 1, 3
        do jj = 1, 3
            do ii = 1, 3
                i = modulo(is+ii-1, nx) + 1
                j = modulo(jc+jj-1, nx) + 1
                k = modulo(kc+kk-1, nx) + 1
                !$acc atomic update
                u(i,j,k) = u(i,j,k) + cf*fxl(l)*wxs(ii)*wyc(jj)*wzc(kk)
                i = modulo(ic+ii-1, nx) + 1
                j = modulo(js+jj-1, nx) + 1
                !$acc atomic update
                v(i,j,k) = v(i,j,k) + cf*fyl(l)*wxc(ii)*wys(jj)*wzc(kk)
                j = modulo(jc+jj-1, nx) + 1
                k = modulo(ks+kk-1, nx) + 1
                !$acc atomic update
                w(i,j,k) = w(i,j,k) + cf*fzl(l)*wxc(ii)*wyc(jj)*wzs(kk)
            enddo
        enddo
    enddo
enddo
end subroutine fib_spread



!===========================================================================
! Cosserat rod: Newmark time integration + Newton-Raphson (core of Update_Pos in FluTAS)
subroutine rod_update(istep, dt)
use param, only: lx
implicit none
integer, intent(in)          :: istep
double precision, intent(in) :: dt
real(fk) :: KG(ndof,ndof), FG(ndof), DispX(ndof), DispU(ndof)
real(fk) :: kele(rows_w,rows_w), mele(rows_w,rows_w)
real(fk) :: Q0(rows_w), Q1(rows_w), finert(rows_w), fgrv(rows_w), fdise(rows_w), fd_CP(rows_w)
real(fk) :: state_x(p1_fordr,7), state_w(p1_fordr,6), state_x_inc(p1_fordr,7)
real(fk) :: uae(p1_fordr,6), state_delta(rows_w)
real(fk) :: NdN_temp(ngps*p1_fordr,nderiv), NdNb_temp(ngpb*p1_fordr,nderiv)
real(fk) :: Jp_s(ngps), Jp_b(ngpb)
real(fk) :: q_old(4), rel_q(4), q_current(4), Drot_vec(3), eN
real(8)  :: strain_gp0(ngps,3)
real     :: elU(2)
integer  :: cone(p1_fordr), ie_con(rows_w), iroco(3)
integer  :: ii, jj, icp, n, nr_step, c
double precision :: c_v, c_a, c_vd, c_ad, length, xcen, ycen, zcen

! Newmark predictor coefficients
c_v  = 1.d0 - dble(Ngamma)/dble(Nbeta)
c_a  = dt*(1.d0 - dble(Ngamma)/(2.d0*dble(Nbeta)))
c_vd = 1.d0/(dble(Nbeta)*dt)
c_ad = 1.d0/(2.d0*dble(Nbeta)) - 1.d0

! start from the state at the previous time step
xll = xfpo;  yll = yfpo;  zll = zfpo
qn1 = qno1;  qn2 = qno2;  qn3 = qno3;  qn4 = qno4
dxdtl = dxdto*c_v + uao*c_a
dydtl = dydto*c_v + vao*c_a
dzdtl = dzdto*c_v + wao*c_a
omg1  = omgo1*c_v + omgdo1*c_a
omg2  = omgo2*c_v + omgdo2*c_a
omg3  = omgo3*c_v + omgdo3*c_a
ua    = -dxdto*c_vd - uao*c_ad
va    = -dydto*c_vd - vao*c_ad
wa    = -dzdto*c_vd - wao*c_ad
omgd1 = -omgo1*c_vd - omgdo1*c_ad
omgd2 = -omgo2*c_vd - omgdo2*c_ad
omgd3 = -omgo3*c_vd - omgdo3*c_ad

nr_step = 0
eN      = 1.0_fk
DispX   = 0.0_fk
do while (eN .gt. TOL)
    nr_step = nr_step + 1
    KG = 0.0_fk
    FG = 0.0_fk

    ! element loop
    do ii = 1, nel
        cone = CONS(ii,1:p1_fordr)
        do icp = 1, p1_fordr
            c = cone(icp)
            do n = 1, 6
                ie_con(6*(icp-1)+n) = 6*(c-1) + n
            enddo
            state_x(icp,:) = (/ xll(c), yll(c), zll(c), qn1(c), qn2(c), qn3(c), qn4(c) /)
            state_w(icp,:) = (/ dxdtl(c), dydtl(c), dzdtl(c), omg1(c), omg2(c), omg3(c) /)
            uae(icp,:)     = (/ ua(c), va(c), wa(c), omgd1(c), omgd2(c), omgd3(c) /)
            fd_CP(6*(icp-1)+1:6*icp) = (/ fxll(c), fyll(c), fzll(c), mxl(c), myl(c), mzl(c) /)
        enddo
        elU = elRangeU(ii,1:2)
        NdN_temp   = NdN_array ((ii-1)*ngps*p1_fordr+1:ii*ngps*p1_fordr,:)
        NdNb_temp  = NdN_arrayb((ii-1)*ngpb*p1_fordr+1:ii*ngpb*p1_fordr,:)
        strain_gp0 = strain_k0((ii-1)*ngps+1:ii*ngps,1:3)
        Jp_s = J_s((ii-1)*ngps+1:ii*ngps)
        Jp_b = J_b((ii-1)*ngpb+1:ii*ngpb)

        ! internal force
        Q0 = 0.0_fk
        call compute_fint(state_x, state_w, elU, ngps, xigs, wgs, Jp_s, NdN_temp, p1_fordr, nno, uKnot, &
                          I_zz, I_yy, J_xx, Ks_y, Ks_z, E_mod, G_mod, Arf, nderiv, Q0, Eflag, &
                          E_modA, G_modA, EI_yy, EI_zz, GJxx, strain_gp0, istep)

        ! element stiffness by finite differences
        kele = 0.0_fk
        state_delta = 0.0_fk
        do jj = 1, rows_w
            state_delta(jj) = Deltaf
            call stateIncrement(state_x, state_delta, p1_fordr, state_x_inc)
            Q1 = 0.0_fk
            call compute_fint(state_x_inc, state_w, elU, ngps, xigs, wgs, Jp_s, NdN_temp, p1_fordr, nno, uKnot, &
                              I_zz, I_yy, J_xx, Ks_y, Ks_z, E_mod, G_mod, Arf, nderiv, Q1, Eflag, &
                              E_modA, G_modA, EI_yy, EI_zz, GJxx, strain_gp0, istep)
            kele(:,jj) = (Q1 - Q0)/Deltaf
            state_delta(jj) = 0.0_fk
        enddo

        ! inertia
        mele = 0.0_fk
        call compute_massE(mInerCP, elU, state_x, nderiv, ngpb, xigb, wgb, Jp_b, NdNb_temp, &
                           p1_fordr, nno, uKnot, mele)
        finert = 0.0_fk
        call compute_finer(mInerCP, elU, state_x, state_w, nderiv, ngpb, xigb, wgb, rhof, &
                           I_yy, I_zz, J_xx, uae, NdNb_temp, p1_fordr, nno, uKnot, finert, &
                           Eflag, rhoA, rhoI_yy, rhoI_zz)
        kele = kele + (1.0_fk/(real(Nbeta,fk)*real(dt,fk)**2))*mele
        Q0   = Q0 + finert

        ! gravity and buoyancy
        if (inertia .eq. 1) then
            mele = 0.0_fk
            call compute_massE(mInerCP2, elU, state_x, nderiv, ngpb, xigb, wgb, Jp_b, NdNb_temp, &
                               p1_fordr, nno, uKnot, mele)
            fgrv = 0.0_fk
            do icp = 1, p1_fordr
                iroco = (/ 6*(icp-1)+1, 6*(icp-1)+2, 6*(icp-1)+3 /)
                fgrv(iroco) = matmul(mele(iroco,iroco), g_vec)
            enddo
            Q0 = Q0 - fgrv
        endif

        ! hydrodynamic (IBM) force
        fdise = 0.0_fk
        call compute_fdise(fd_CP, state_x(:,1:3), elU, nderiv, ngpb, xigb, wgb, NdNb_temp, &
                           p1_fordr, nno, uKnot, Arf, fdise)
        Q0 = Q0 + fdise

        ! assembly
        FG(ie_con) = FG(ie_con) + Q0
        KG(ie_con,ie_con) = KG(ie_con,ie_con) + kele
    enddo

    ! solve KG*DispU = -FG (free-free fiber: all dofs active)
    DispU = -FG
    call solve_fk(ndof, KG, DispU)
    DispX = DispX + DispU

    ! update positions, rotations, velocities and accelerations of the control points
    do ii = 1, nno
        n = 6*(ii-1)
        xll(ii) = xll(ii) + DispU(n+1)
        yll(ii) = yll(ii) + DispU(n+2)
        zll(ii) = zll(ii) + DispU(n+3)
        q_old    = (/ real(qn1(ii),fk), real(qn2(ii),fk), real(qn3(ii),fk), real(qn4(ii),fk) /)
        Drot_vec = DispU(n+4:n+6)
        call q_from_Rotv(Drot_vec, rel_q)
        call quad_prod(q_old, rel_q, q_current)
        qn1(ii) = q_current(1)
        qn2(ii) = q_current(2)
        qn3(ii) = q_current(3)
        qn4(ii) = q_current(4)

        dxdtl(ii) = dxdto(ii)*c_v + uao(ii)*c_a    + dble(Ngamma)*c_vd*DispX(n+1)
        dydtl(ii) = dydto(ii)*c_v + vao(ii)*c_a    + dble(Ngamma)*c_vd*DispX(n+2)
        dzdtl(ii) = dzdto(ii)*c_v + wao(ii)*c_a    + dble(Ngamma)*c_vd*DispX(n+3)
        omg1(ii)  = omgo1(ii)*c_v + omgdo1(ii)*c_a + dble(Ngamma)*c_vd*DispX(n+4)
        omg2(ii)  = omgo2(ii)*c_v + omgdo2(ii)*c_a + dble(Ngamma)*c_vd*DispX(n+5)
        omg3(ii)  = omgo3(ii)*c_v + omgdo3(ii)*c_a + dble(Ngamma)*c_vd*DispX(n+6)

        ua(ii)    = -dxdto(ii)*c_vd - uao(ii)*c_ad    + c_vd/dt*DispX(n+1)
        va(ii)    = -dydto(ii)*c_vd - vao(ii)*c_ad    + c_vd/dt*DispX(n+2)
        wa(ii)    = -dzdto(ii)*c_vd - wao(ii)*c_ad    + c_vd/dt*DispX(n+3)
        omgd1(ii) = -omgo1(ii)*c_vd - omgdo1(ii)*c_ad + c_vd/dt*DispX(n+4)
        omgd2(ii) = -omgo2(ii)*c_vd - omgdo2(ii)*c_ad + c_vd/dt*DispX(n+5)
        omgd3(ii) = -omgo3(ii)*c_vd - omgdo3(ii)*c_ad + c_vd/dt*DispX(n+6)
    enddo

    ! residual (evaluated before the update, as in FluTAS)
    eN = sqrt(sum(FG**2))
    if (nr_step .ge. nr_max) then
        write(*,*) "Fiber: Newton-Raphson not converged, residual", real(eN,8)
        stop
    endif
enddo
nr_last = nr_step
if (fib_out .gt. 0) then
    if (mod(istep,fib_out) .eq. 0) call fib_write_shape(istep)
endif

! store the new state
xfpo = xll;  yfpo = yll;  zfpo = zll
qno1 = qn1;  qno2 = qn2;  qno3 = qn3;  qno4 = qn4
dxdto = dxdtl; dydto = dydtl; dzdto = dzdtl
omgo1 = omg1;  omgo2 = omg2;  omgo3 = omg3
uao = ua; vao = va; wao = wa
omgdo1 = omgd1; omgdo2 = omgd2; omgdo3 = omgd3

! periodicity: shift the whole fiber when its centre leaves the box
xcen = xll((nno+1)/2)
ycen = yll((nno+1)/2)
zcen = zll((nno+1)/2)
if (xcen .ge. lx) then
    xll = xll - lx; xfpo = xfpo - lx
elseif (xcen .lt. 0.d0) then
    xll = xll + lx; xfpo = xfpo + lx
endif
if (ycen .ge. lx) then
    yll = yll - lx; yfpo = yfpo - lx
elseif (ycen .lt. 0.d0) then
    yll = yll + lx; yfpo = yfpo + lx
endif
if (zcen .ge. lx) then
    zll = zll - lx; zfpo = zfpo - lx
elseif (zcen .lt. 0.d0) then
    zll = zll + lx; zfpo = zfpo + lx
endif

! log: step, Newton iterations, length, centre position and velocity, position of the last end
if (mod(istep,fib_log) .eq. 0) then
    length = 0.d0
    do ii = 1, nno-1
        length = length + sqrt((xll(ii+1)-xll(ii))**2 + (yll(ii+1)-yll(ii))**2 + (zll(ii+1)-zll(ii))**2)
    enddo
    open(unit=56,file='./output/fiber_log.dat',position='append')
    write(56,'(I9,I4,10ES16.7)') istep, nr_last, length, xll((nno+1)/2), yll((nno+1)/2), zll((nno+1)/2), &
                                 dxdtl((nno+1)/2), dydtl((nno+1)/2), dzdtl((nno+1)/2), xll(nno), yll(nno), zll(nno)
    close(56)
endif
end subroutine rod_update



!===========================================================================
! Gaussian elimination with partial pivoting: A x = b (b overwritten by x)
subroutine solve_fk(n, A, b)
implicit none
integer, intent(in)     :: n
real(fk), intent(inout) :: A(n,n), b(n)
integer  :: i, k, p
real(fk) :: f, rowt(n), bt
do k = 1, n-1
    p = k - 1 + maxloc(abs(A(k:n,k)),1)
    if (p .ne. k) then
        rowt = A(k,:); A(k,:) = A(p,:); A(p,:) = rowt
        bt = b(k); b(k) = b(p); b(p) = bt
    endif
    do i = k+1, n
        f = A(i,k)/A(k,k)
        A(i,k:n) = A(i,k:n) - f*A(k,k:n)
        b(i) = b(i) - f*b(k)
    enddo
enddo
do i = n, 1, -1
    b(i) = (b(i) - sum(A(i,i+1:n)*b(i+1:n)))/A(i,i)
enddo
end subroutine solve_fk

subroutine solve_dp(n, m, A, B)
implicit none
integer, intent(in)             :: n, m
double precision, intent(inout) :: A(n,n), B(n,m)
integer :: i, k, p
double precision :: f, rowt(n), bt(m)
do k = 1, n-1
    p = k - 1 + maxloc(abs(A(k:n,k)),1)
    if (p .ne. k) then
        rowt = A(k,:); A(k,:) = A(p,:); A(p,:) = rowt
        bt = B(k,:); B(k,:) = B(p,:); B(p,:) = bt
    endif
    do i = k+1, n
        f = A(i,k)/A(k,k)
        A(i,k:n) = A(i,k:n) - f*A(k,k:n)
        B(i,:) = B(i,:) - f*B(k,:)
    enddo
enddo
do i = n, 1, -1
    B(i,:) = (B(i,:) - matmul(A(i,i+1:n),B(i+1:n,:)))/A(i,i)
enddo
end subroutine solve_dp



!===========================================================================
! Output: control points (text) + full state for restart (binary)
subroutine fib_output(t)
implicit none
integer, intent(in) :: t
character(len=40) :: namefile
call fib_write_shape(t)
write(namefile,'(a,i8.8,a)') './output/fibstate_',t,'.dat'
open(unit=57,file=namefile,form='unformatted',access='stream',status='replace')
write(57) nno, xll, yll, zll, qn1, qn2, qn3, qn4, dxdtl, dydtl, dzdtl, omg1, omg2, omg3, &
          ua, va, wa, omgd1, omgd2, omgd3
close(57)
end subroutine fib_output

! control points (text), also written every fib_out steps for movies
subroutine fib_write_shape(t)
implicit none
integer, intent(in) :: t
character(len=40) :: namefile
integer :: i
write(namefile,'(a,i8.8,a)') './output/fib_',t,'.dat'
open(unit=57,file=namefile,form='formatted',status='replace')
write(57,'(a)') '# x y z q1 q2 q3 q4 dxdt dydt dzdt fx fy fz (control points)'
do i = 1, nno
    write(57,'(13ES24.15)') xll(i), yll(i), zll(i), qn1(i), qn2(i), qn3(i), qn4(i), &
                            dxdtl(i), dydtl(i), dzdtl(i), fxll(i), fyll(i), fzll(i)
enddo
close(57)
end subroutine fib_write_shape

subroutine fib_read_restart(t)
implicit none
integer, intent(in) :: t
character(len=40) :: namefile
integer :: n
write(namefile,'(a,i8.8,a)') './output/fibstate_',t,'.dat'
open(unit=57,file=namefile,form='unformatted',access='stream',status='old')
read(57) n
if (n .ne. nno) stop "Fiber restart: number of control points does not match"
read(57) xll, yll, zll, qn1, qn2, qn3, qn4, dxdtl, dydtl, dzdtl, omg1, omg2, omg3, &
         ua, va, wa, omgd1, omgd2, omgd3
close(57)
write(*,*) "Fiber state read from ", namefile
end subroutine fib_read_restart

end module fibers
