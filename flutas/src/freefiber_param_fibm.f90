module mod_param_fibm
!use decomp_2d
!use mod_linspace
implicit none


integer, parameter :: ndims = 2
integer, dimension(ndims), parameter :: dims = (/1,2/)
integer, parameter :: itot=2, jtot=112, ktot=140
!integer, parameter :: itot=2,jtot=200,ktot=50
!integer, parameter :: itot=2,jtot=320,ktot=320
!integer, parameter :: itot=2,jtot=160,ktot=160
!---------
!integer, parameter :: itot=2,jtot=320,ktot=64
!--------
!integer, parameter :: itot=2,jtot=160,ktot=96
!integer, parameter :: itot=2,jtot=340,ktot=108
integer, parameter :: it1=itot+1, jt1 = jtot+1, kt1 = ktot+1
integer, parameter :: imax=itot/dims(1), jmax=jtot/dims(2),kmax=ktot
integer, parameter :: i1 = imax+1, j1 = jmax+1, k1 = kmax+1
real, parameter :: lx = 0.1, ly = 5.6, lz = 7.0
!real, parameter :: lx = 0.05,ly = 8.0,lz = 1.61
!real, parameter :: lx = 0.1,ly = 10.0,lz = 6.0
!real, parameter :: lx = 6.0,ly = 10.0,lz = 6.0
real, parameter :: dxi = itot/lx, dyi = jtot/ly, dzi = ktot/lz
real, parameter :: dx = 1./dxi, dy = 1./dyi, dz = 1./dzi
real, parameter :: pi = acos(-1.)
real, parameter :: kappa= 1000000.
!
! insert here info on the number of particles
!
integer, parameter :: np=1, npmax=3
!integer, parameter :: np=2, npmax=5
!integer, parameter :: np=5, npmax=9
integer, parameter :: npi=1,npj=1,npk=1
integer,parameter::nrow=3,ncolumn=2
!v--------------------------------------------------------------------------------------
!Input: Length, diam, rhof, nel, deg_ele, Nbeta-Ngamma, E_mod, G_mod, Ks_y, g_vec, BCs
REAL, PARAMETER        :: L_fibr= 1.4, diametr=0.02, aspectratio=diametr/L_fibr
!REAL, PARAMETER       :: L_fibr=0.7, diametr=0.0212, aspectratio=diametr/L_fibr
!REAL, PARAMETER       :: L_fibr= 1.0,  diametr=0.1, aspectratio=diametr/L_fibr
REAL, PARAMETER        :: rhof  = 5.0 !*pi*.25*(aspectratio**2.) 
INTEGER, PARAMETER     :: nel   = 28, deg_ele=0, add_ele=0 !v num of elements should be even
!REAL(16), PARAMETER    :: E_mod= 7.8228D+04, G_mod = 3.0088D+04, &
!REAL(16), PARAMETER    :: E_mod= 5.7367e+04, G_mod = 2.6076e+04, &
!REAL(16), PARAMETER    :: E_mod= 5.2152e+04, G_mod = 2.0058D+04, & !2.173e+02, &
!----
!REAL(16), PARAMETER    :: E_mod=  1.4D+06, G_mod = 0.4667D+06, &
!----
!REAL(16), PARAMETER    :: E_mod= 4.572D+05, G_mod = 1.7585D+05, & !
!REAL(16), PARAMETER    :: E_mod= 4.572D+03, G_mod = 1.7585D+03, & !<0.7 more NRs
!REAL(16), PARAMETER    :: E_mod= 4.572D+05, G_mod = 1.7585D+05, & !~0.95 (rho5,bval45)
!REAL(16), PARAMETER    :: E_mod= 2.8D+06, G_mod = 1.0769D+06, & !~0.7 small defc
!REAL(16), PARAMETER    :: E_mod= 4.572D+06, G_mod = 1.7585D+06, & !~0.65 small defc
!REAL(16), PARAMETER    :: E_mod= 7.76e+06, G_mod = 2.9846D+06, & !~0.6 def
!REAL(16), PARAMETER    :: E_mod= 1.572D+07, G_mod = 0.6046D+07, & !~hv conv issue
!REAL(16), PARAMETER    :: E_mod= 4.572D+08, G_mod = 1.7585D+08, & !<conv issue 4m strt
!----
REAL(16), PARAMETER    :: E_mod= 4.572D+05, G_mod = 1.7585D+05, & !
                          nu   = (E_mod/(2. * G_mod)) - 1.0, &
                          Ks_y = 5.0/6.0, & !6.0*(1.0 + nu)/(7.0 +6.0*nu), &
                          Ks_z = Ks_y
!REAL(16), PARAMETER    :: g_vec(3) = (/ 0.0, 0.5, 0.0 /)
REAL(16), PARAMETER    :: g_vec(3) = (/ 0.0, 0.0, 0.0 /)
REAL                   :: Nbeta=0.25, Ngamma=0.5
REAL(8)                :: uKnoti(0:3) = (/ 0., 0., 1., 1. /)
INTEGER                :: p_ordr=1, noCP=2
real, parameter        :: ds=L_fibr/(1.*nel)
! ControlPoints and DOFs
INTEGER, PARAMETER     :: nno = nel + (1+deg_ele) + add_ele*nel, ndof=6*nno, &
                          p1_fordr = 1+1+deg_ele+add_ele, rows_w=6*p1_fordr
!integer, parameter    :: nl= nel+(1+deg_ele), nderiv = 2 ! nl=nno
integer, parameter     :: nl=nno, nderiv = 2 ! nl=nno
! End point forces/moments and Distributed force related
!INTEGER, PARAMETER     :: bc0_dof(6) = (/ 1, 2, 3, 4, 5, 6 /), &
INTEGER, PARAMETER     :: bc0_dof(3) = (/ 1, 2, 3 /), &
                          bcL_dof(6) = (/ ndof-5, ndof-4, &
                          ndof-3, ndof-2, ndof-1, ndof /), &
                          Fext_vec(3)  = (/0, 0, 0 /), &
                          Mext_vec(3)  = (/0, 0, 0/), &
                          Fext0_vec(3) = (/0, 0, 0/), &
                          Mext0_vec(3) = (/0, 0, 0/), &
                          fmD_vec(6)   = (/0, 0, 0, 0, 0, 0/)
! Adjust nir as per BCs
character(1),parameter :: fboundary="f"  !c=clamped,h=hinged, and f=free 
!INTEGER, PARAMETER     :: nir=ndof-ubound(bc0_dof,1) !-ubound(bcL_dof,1), & 
INTEGER, PARAMETER     :: nir=ndof!-ubound(bcL_dof,1)  
! NR Related
REAL(16), PARAMETER    :: TOL=1.0D-07, Deltaf = 1.0D-10
INTEGER                :: l_step=0, nr_step=0
! Quadrature-rule related
!INTEGER, PARAMETER    :: noGPs=CEILING(p1_fordr/2.0), & !Exact 
INTEGER, PARAMETER     :: noGPs=p1_fordr-1, & !reduced
                          ngps=noGPs, ngpb=noGPs
! Bending and axial rigidity                          
LOGICAL, PARAMETER     :: Eflag = .true.
!LOGICAL, PARAMETER    :: Eflag = .false.
REAL(16), PARAMETER    :: E_modA=rhof*2.0D+04, G_modA=rhof*0.76923D+04, &
                          EI_yy=rhof*5.0D-01, EI_zz=EI_yy, GJxx=rhof*3.846D-01
REAL(16), PARAMETER    :: rhoA=rhof*3.1416D-04, rhoI_yy=rhof*7.854D-9, rhoI_zz=rhoI_yy
!v--------------------------------------------------------------------------------------
real,parameter:: gammafilament=0.0000, alphafilament=-0.0D+06, fr=0.0,rhofilament=pi*.25*(aspectratio**2.)!*rhof
!better define rhofilament inside common_fibm so that area confusion is gone
!v alphafilament used in forcing; gammafilamaent and fr are used in Update_Pos.f90
real,parameter:: fdiam=aspectratio
integer, parameter :: nfd = 2             !v ??used in 'suspension.f90' many times
!character(1),parameter::TypeFlow ="d"     !v which type of flow
character(1),parameter::solver="f"        !v which solver
integer,parameter::filsub=1               !v ??used only in 'Update_Pos.f90'
integer,parameter::flforce=1              !v used in 'forcing.f90'
!integer,parameter::sharpdelta=0          !v used in 'interp_spred' (maybe related to IBM detla)

!*************************************boundary condition of first point of filament***********************************
!
integer, parameter :: send_real = 6+21*nl !,send_int = 1+nqmax    !v used in 6 files/routines, everytime for MPI purposes
! amount of data to be send from mstr to slve: see common file
!
! type of initial velocity field (see init.f90)
! iniu = 'cou' --> plane Couette flow
!      = 'poi' --> plane Poiseuille flow
!      = 'zer' --> zero velocity everywhere
!      = 'log' --> logarithmic profile + random noise
!!character(len=3), parameter :: iniu = 'poi' !      !v type of initial velocity field
! type of BC z = lz is set to free-slip if
! isfreeslip is true
!!logical, parameter :: isfreeslip = .false.
!!real, parameter :: diam=1. ! lengths are scaled with the particle diameter
!!real, parameter :: Reb=200.0
!!real, parameter ::Rep=Reb/lz
!!real, parameter ::tw=(180./2800.)**2.
! desired bulk Reynolds number based on superficial bulk velocity and channel height
!!real, parameter :: visc0=1./Rep
! characteristic velocity scale and the sphere diameter as characteristic length scale
!!real, parameter :: bulk_v_sup = Rep*visc0 !equal to 1. for chosen value of visc
!!real, parameter :: gacc = 0.!-9.81*diam/bulk_v_sup**2.!diam/vscale**2.
!!real, parameter :: gaccx = gacc*cos(pi/2.), gaccy = gacc*cos(pi/2.), &
!!                        gaccz = gacc*sin(pi/2.)
!!real, parameter :: ratiorho = 8.
real, parameter :: radius = .5, offset =(( sqrt(3.*(1.5**2)) )/dxi + 0.01/dxi)*.2
!!real, parameter :: volp = (4./3.)*pi*radius**3.       !!!!!v for shpere not filament
!!real, parameter :: mominert = (2./5.)*volp*radius**2.
!!real, parameter :: lref = 1., uref = sqrt(lref), tref = lref/uref
!!real,parameter:: shear=1
!********************************************************************filament forcing***********************************************
integer,parameter::inertia=0   ! 0 filament is neutrally bouyant and 1 when there is density difference between filament & fluid
!integer,parameter::inertia=0   ! 0 filament is neutrally bouyant and 1 when there is density difference between filament & fluid
integer,parameter::ifcurved=0   ! 1 when curved shape will be given to initial shape of the filaments
integer,parameter::ifattached=0 ! 1 when filaments are arranged near wall
!***********************lubrication parametrs*******************************
!!character(len=3), parameter :: lubr = 'aid'
integer,parameter :: nzero= 3, iffric=0
integer, parameter :: pointmax = int((np*nl/real(itot*jtot*ktot))*(nzero**3.)*1.5)+30
real,parameter :: cutoff=1.25*aspectratio,dsurf=0.,mufric=.05
!!real,parameter :: cl= 5.,betacol=15.,dist0=0.!col is non dimensional Beta*De in the Morse potential model
!!real, parameter :: retrac = 0.3/dxi
!
!integer, parameter :: nl = nint((pi/3.)*(12.*(((radius-retrac)*dxi)**2)+1.))

!     nr lfp's of second shell
      integer NL2
      parameter(NL2=nint((pi/3.)*(12.*(((radius-1.5/dxi)*dxi)**2.) + 1. )))
!     nr lfp's of third shell
      integer NL3
      parameter(NL3=nint((pi/3.)*(12.*(((radius-2.5/dxi)*dxi)**2.) + 1. )))
!     nr lfp's of fourth shell
      integer NL4
      parameter(NL4=nint((pi/3.)*(12.*(((radius-3.5/dxi)*dxi)**2.) + 1. )))
!     nr lfp's of 1st-4th shell
      integer NLtot
      parameter(NLtot=NL+NL2+NL3+NL4)
real, parameter :: solidity = (1.*np)*(4./3.)*pi*(radius**3)/(lx*ly*lz)
!
character(len=5), parameter :: datadir = 'data/'
!
!!real, parameter, dimension(3,2) :: rkcoeff = reshape((/ 32./60., 25./60., 45./60., 0., -17./60., -25./60. /), shape(rkcoeff))
!!real, parameter, dimension(3) :: rkcoeffab = rkcoeff(:,1)+rkcoeff(:,2)
!
!Output parameters
!integer,parameter :: ioutchk = 10/1, iout1d = 1/1,iout2d = 20 ,ioutfld = 200 , ioutfil=10000000
!!integer,parameter :: ioutchk=10/1, iout1d=100/1, iout2d=100, ioutfld=500, ioutfil=10000000

!
! set collision parameters
!
!
end module mod_param_fibm
