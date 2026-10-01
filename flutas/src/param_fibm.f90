module mod_param_fibm
!use decomp_2d
!use mod_linspace
implicit none


integer, parameter :: ndims = 2
integer, dimension(ndims), parameter :: dims = (/2,2/)   !<---------
integer, parameter :: itot=128,jtot=128,ktot=128           !<---------
integer, parameter :: it1=itot+1, jt1 = jtot+1, kt1 = ktot+1
integer, parameter :: imax=itot/dims(1), jmax=jtot/dims(2),kmax=ktot
integer, parameter :: i1 = imax+1, j1 = jmax+1, k1 = kmax+1
real, parameter :: lx = 6.28,ly = 6.28,lz = 6.28           !<---------
real, parameter :: dxi = itot/lx, dyi = jtot/ly, dzi = ktot/lz
real, parameter :: dx = 1./dxi, dy = 1./dyi, dz = 1./dzi
real, parameter :: pi = acos(-1.)
real, parameter :: kappa= 1000000.
!
! insert here info on the number of particles
!
integer, parameter :: np=1, npmax=3
!integer, parameter :: np=5, npmax=9
integer, parameter :: npi=1,npj=1,npk=1
integer,parameter  ::nrow=3,ncolumn=2
character(1),parameter :: c_type= "c" 
REAL, PARAMETER        :: L_fibr= 1., diametr=0.025, aspectratio=diametr/L_fibr
REAL, PARAMETER        :: rhof  = 10. 
INTEGER, PARAMETER     :: nel   = 12, deg_ele=0, add_ele=0 !num of elements should be even
REAL(16), PARAMETER    :: E_mod = 218.0D+03, G_mod = 83.846D+03, & !
                          nu    = (E_mod/(2. * G_mod)) - 1.0, &
                          Ks_y  = 5.0/6.0, & !6.0*(1.0 + nu)/(7.0 +6.0*nu), &
                          Ks_z  = Ks_y
REAL(16), PARAMETER    :: g_vec(3) = (/ 0.0, -9.81, 0.0 /)
!REAL(16), PARAMETER   :: g_vec(3) = (/ 0.0, 0.0, 0.0 /)
REAL                   :: Nbeta=0.25, Ngamma=0.5
REAL(8)                :: uKnoti(0:3) = (/ 0., 0., 1., 1. /)
INTEGER                :: p_ordr=1, noCP=2
! ControlPoints and DOFs
INTEGER, PARAMETER     :: nno=nel+ (1+deg_ele) + add_ele*nel, ndof=6*nno, &
                          p1_fordr=1+1+deg_ele+add_ele, rows_w=6*p1_fordr
! parametric points
integer, parameter     :: nxie=10, nl=nel*nxie, nderiv=2   ! nl = points actually used (was nel*nxie+1)
real                   :: nxi_tvec(nl), nxi_vec(nl)
real, parameter        :: ds = L_fibr/(1.*nl)               ! arc length of each Lagrangian point
! End point forces/moments and Distributed force related
INTEGER, PARAMETER    :: bc0_dof(6) = (/ 1, 2, 3, 4, 5, 6 /), &
!INTEGER, PARAMETER   :: bc0_dof(3) = (/ 1, 2, 3 /), &
                          bcL_dof(6) = (/ ndof-5, ndof-4, &
                          ndof-3, ndof-2, ndof-1, ndof /), &
                          Fext_vec(3)  = (/0, 0, 0 /), &
                          Mext_vec(3)  = (/0, 0, 0/), &
                          Fext0_vec(3) = (/0, 0, 0/), &
                          Mext0_vec(3) = (/0, 0, 0/), &
                          fmD_vec(6)   = (/0, 0, 0, 0, 0, 0/)

character(1),parameter :: fboundary="f"  !c=clamped,h=hinged, and f=free 
INTEGER, PARAMETER     :: nir=ndof
! NR Related
REAL(16), PARAMETER    :: TOL=1.0D-07, Deltaf = 1.0D-09
INTEGER                :: l_step=0, nr_step=0
!INTEGER, PARAMETER    :: noGPs=CEILING(p1_fordr/2.0), & !Exact 
INTEGER, PARAMETER     :: noGPs=p1_fordr-1, &     !reduced
                          ngps=noGPs, ngpb=noGPs
LOGICAL, PARAMETER     :: Eflag = .false.
REAL(16), PARAMETER    :: E_modA=25.6, G_modA=9.8459, &
                          EI_yy=0.001, EI_zz=EI_yy, GJxx=7.6921D-04
REAL(16), PARAMETER    :: rhoA=rhof*4.9807D-04, rhoI_yy=rhof*1.9175D-08, rhoI_zz=rhoI_yy

!v--------------------------------------------------------------------------------------
real,parameter::gammafilament=0.0000,alphafilament=-0.0D+06,fr=0.0,rhofilament=pi*.25*(aspectratio**2.)
real,parameter::fdiam=aspectratio
integer, parameter :: nfd = 2             
character(1),parameter::solver="f"        
integer,parameter::filsub=1               
integer,parameter::flforce=1              
!integer,parameter::sharpdelta=0          

!*************************************boundary condition of first point of filament***********************************
integer, parameter :: send_real = 6+57*nno !,send_int = 1+nqmax    !v used in 6 files/routines, everytime for MPI purposes
real, parameter :: radius = .5, offset =(( sqrt(3.*(1.5**2)) )/dxi + 0.01/dxi)*.2
!********************************************************************filament forcing***********************************************
integer,parameter::inertia=1   ! 0 filament is neutrally bouyant and 1 when there is density difference between filament & fluid
integer,parameter::ifcurved=0   ! 1 when curved shape will be given to initial shape of the filaments
integer,parameter::ifattached=0 ! 1 when filaments are arranged near wall
!***********************lubrication parametrs*******************************
!!character(len=3), parameter :: lubr = 'aid'
integer,parameter :: nzero= 3, iffric=0
integer, parameter :: pointmax = int((np*nl/real(itot*jtot*ktot))*(nzero**3.)*1.5)+30
real,parameter :: cutoff=1.25*aspectratio,dsurf=0.,mufric=.05
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
end module mod_param_fibm
