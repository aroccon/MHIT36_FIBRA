! Precision used by the fiber (Cosserat rod) solver, which runs on the CPU.
! Default: double precision (the finite-difference stiffness step Deltaf is set accordingly).
! With -DFIB_QUAD the solver uses REAL(16) as in FluTAS, if the compiler supports it.
module fib_prec
#ifdef FIB_QUAD
    integer, parameter :: fk_quad = selected_real_kind(30)
    integer, parameter :: fk = merge(fk_quad, kind(1.d0), fk_quad > 0)
#else
    integer, parameter :: fk = kind(1.d0)
#endif
end module fib_prec



! Fiber parameters (single fiber, ported from FluTAS param_fibm.f90)
! Lengths are in the same units as the MHIT36 box (lx = 2*pi)
module fib_param
    use fib_prec
    implicit none
    ! geometry and density of the fiber
    double precision, parameter :: L_fibr = 1.d0          ! fiber length
    double precision, parameter :: diametr = 0.025d0      ! fiber diameter
    real, parameter             :: rhof = 10.             ! fiber density (only used when Eflag=.false.)
    character(1), parameter     :: c_type = "c"           ! cross section: c=circular, s=square
    ! initial position of the fiber centre and direction of its axis
    double precision, parameter :: xc0 = 3.14159265358979d0
    double precision, parameter :: yc0 = 3.14159265358979d0
    double precision, parameter :: zc0 = 3.14159265358979d0
    double precision, parameter :: dir0(3) = (/ 0.d0, 0.d0, 1.d0 /)
    ! .false.: material frame = global frame at t=0, as in FluTAS (the rod axis is then
    !          the material z axis, so axial/bending/torsion get the FluTAS stiffnesses)
    ! .true. : material axis 1 aligned with the fiber tangent at t=0
    logical, parameter          :: align_frame = .true.
    ! IGA discretisation
    integer, parameter          :: nel = 12, deg_ele = 0  ! elements, degree elevation
    integer, parameter          :: p1_fordr = 2 + deg_ele ! control points per element
    integer, parameter          :: nno = nel + 1 + deg_ele, ndof = 6*nno
    integer, parameter          :: nir = ndof             ! free-free fiber: all dofs active
    integer, parameter          :: rows_w = 6*p1_fordr
    integer, parameter          :: nderiv = 2
    integer, parameter          :: ngps = p1_fordr-1, ngpb = ngps ! reduced integration
    ! Lagrangian (IBM) points: nxie per element, at the centre of nl equal segments
    integer, parameter          :: nxie = 10, nl = nel*nxie
    double precision, parameter :: ds = L_fibr/dble(nl)
    ! material properties
    real(fk), parameter         :: E_mod = 218.0D+03, G_mod = 83.846D+03
    real(fk), parameter         :: Ks_y = 5.0/6.0, Ks_z = Ks_y
    logical, parameter          :: Eflag = .false.        ! .true.: use the rigidities below
    real(fk), parameter         :: E_modA = 25.6, G_modA = 9.8459, &
                                   EI_yy = 0.001, EI_zz = EI_yy, GJxx = 7.6921D-04
    real(fk), parameter         :: rhoA = rhof*4.9807D-04, rhoI_yy = rhof*1.9175D-08, rhoI_zz = rhoI_yy
    ! gravity (buoyancy term active when inertia=1)
    integer, parameter          :: inertia = 0                     ! 1: gravity and buoyancy on the fiber
    real(fk), parameter         :: g_vec(3) = (/ 0.0, 0.0, 0.0 /)
    ! Newmark parameters
    real, parameter             :: Nbeta = 0.25, Ngamma = 0.5
    ! Newton-Raphson: tolerance, max iterations and step of the finite-difference stiffness
    real(fk), parameter         :: TOL = 1.0D-07
    integer, parameter          :: nr_max = 50
    real(fk), parameter         :: Deltaf = merge(1.0D-09, 1.0D-07, precision(1._fk) > 20)
    ! output frequency of the fiber log (time steps)
    integer, parameter          :: fib_log = 100
    ! output frequency of the fiber shape files output/fib_XXXXXXXX.dat (time steps, 0 = only at dump)
    integer, parameter          :: fib_out = 100
    ! .false.: fluid-fiber coupling off (no interpolation/spreading, zero hydrodynamic load):
    !          the fiber evolves in vacuum, to validate the structural solver alone.
    !          For a test without external load also set g_vec = 0 (or inertia = 0).
    logical, parameter          :: fib_coupling = .true.
    ! initial transverse velocity with the shape of a free-free bending mode (pert_mode = 1..4)
    ! (amplitude at the fiber ends, along pert_dir projected normal to the fiber); 0 = fiber at rest
    double precision, parameter :: pert_amp = 0.d0
    double precision, parameter :: pert_dir(3) = (/ 0.d0, 1.d0, 0.d0 /)
    integer, parameter          :: pert_mode = 1
end module fib_param
