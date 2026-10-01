module mod_Update_Pos
use mpi
use mod_param_fibm,		only: i1, j1, k1, nl, np, ndims, ndof, nir, p1_fordr, &
												filsub, l_step, nr_step, Ngamma, Nbeta, PI, &
												TOL, nel, ngps, ngpb, nno, Ks_y, Ks_z, &
												E_mod, G_mod, nderiv, rows_w, Deltaf, &
												inertia, rhof, g_vec, nir, datadir, L_fibr, dx, &
												dy, dz, ds, dims, lx, ly, lz, offset, send_real, &
												Eflag, rhoA, rhoI_yy, rhoI_zz, E_modA, EI_yy, &
												EI_zz, G_modA, GJxx
use mod_common_fibm, only: npmax, es, pmax, ap, sp, eN, CONS, cone, &
                        ie_con, icp, jcp, elU, c1, c2, J_s, J_b, &
                        NdN_array, uKnot, I_zz, I_yy, J_xx, elRangeU, &
                        mInerCP, NdN_arrayb, vec3, vec6, ir_dof, lcounter, &
                        Arf, iroco, dt_step, npmstr, xigs, xigb, wgs, wgb, &
                        mInerCP2, NdN_temp, NdNb_temp, Jp_s, Jp_b, strain_k0
												!forcex, forcey, forcez
use mod_ext_force,   only: ext_force
use mod_common_mpi,  only: MPI_STATUS_SIZE, MPI_INTEGER, neighbor, comm_cart, &
												status, ierr, MPI_SUM, myid, coords
use mod_param,         only: rho_sp
use profiler
!use mod_interp_spread, only:    
!use mod_comp_tension    
!use mod_collisions      
!v----
!use mod_commonIGA
use bspline 
use mod_linspace
use mod_initialSetup
use mod_cosseratFun
use mod_linalg

implicit none
private
public Update_Pos
contains
!
subroutine Update_Pos(unew, vnew, wnew, dt, istep, time, forcex, forcey, forcez)
implicit none
!
real(8), intent(in), dimension(0:i1, 0:j1, 0:k1)  :: unew, vnew, wnew
real(8), intent(in)                               :: dt,time !,f_t12
integer, intent(in)                               :: istep
real, intent(inout), dimension(0:i1, 0:j1, 0:k1)  :: forcex, forcey, forcez  
!
integer                                           :: do1, do2, do3
integer:: p,g, pp
real,dimension(NL):: aeex,awwx,sm,dxs,dys,dzs  
real :: maxerrorx,maxerrort,Cont_x,Cont_y,distx,disty,length
real::sf,ur
integer::lstart,ntp,nbp,nt,nbb

type neighbour2                 
  real :: x,y,z
  real, dimension(nno) :: xll,yll,zll!, qN1, qN2, qN3, qN4
end type neighbour2
type(neighbour2), dimension(1:npmax,0:8) :: anb2
                                
integer :: tag,l,k,i            
integer :: nrrequests           
integer :: arrayrequests(1:30)
integer :: arraystatuses(MPI_STATUS_SIZE,1:30) 

real :: ax,ay   

real :: leftbound,rightbound,frontbound,backbound 
integer :: nbrec2,nbrec,nbsend 
integer, dimension(ndims) :: proccoords 
integer :: procrank                    
integer :: counter !delete after checking things

integer :: count_mstr,count_slve,count_mstr_all,count_slve_loc 
logical :: found_mstr

integer, dimension(0:8) :: pmax_nb
integer, dimension(1:npmax,0:8) :: mslv_nb, newmaster_nb   
integer :: idp,idp_nb,nb                                    
real :: boundleftnb,boundrightnb,boundfrontnb,boundbacknb   
integer,parameter :: nprocs = dims(1)*dims(2)               
! Local and Global arrays
REAL(16)              :: KG(ndof,ndof), FG(ndof), &
                         FGr(nir), DUr(nir), &
                         kele(6*p1_fordr, 6*p1_fordr), &
                         mele(6*p1_fordr, 6*p1_fordr), &
                         finert(6*p1_fordr), &
                         fgrv(6*p1_fordr), fdise(6*p1_fordr), &
                         fdval(es)                            
REAL(16)              :: Q0(6*p1_fordr), Q1(6*p1_fordr)
integer, parameter    :: dp = kind(1.d0)!, KGr(nir,nir)
REAL*8                :: KGr(nir, nir) , WORK(nir)
REAL(8)               :: strain_kp0(nel*ngps,3), strain_gp0(ngps,3)
REAL(8)               :: vecp_ex(np,3) != &
												 !RESHAPE( (/ 0.0, 0.0, 1.0,  0.0, 0.0, -1.0 /), (/np,3/))
REAL(16)              :: state_x(p1_fordr,7), state_w(p1_fordr,6),&
                         state_x_inc(p1_fordr, 7), uae(p1_fordr,6), &
                         state_x0(p1_fordr,3), dr0(3), &
                         state_delta(6*p1_fordr), DispX(ndof), &
                         DispU(ndof), Disp(ndof), Drot_vec(3), &
                         rel_q(4), q_old(4), q_current(4), &                          
                         rel_q_temp(4), un(nno,6), un_old(nno,6),   & 
                         ua(nno,6), ua_old(nno,6), fd_CP(6*p1_fordr)!, Drot_vec(3)
INTEGER               :: ii, jj, inx, info, IPIV(nir) 
external                DGETRF, DGETRI,  SGETRF, SGETRI    

pmax_nb(0)        = pmax  !v pmax defined in mod_common, initparticles -> pmax = i!count_mstr+count_slve
!$omp workshare
mslv_nb(1:pmax,0) = ap(1:pmax)%mslv
!$omp end workshare
do nb=1,8
  nbsend = 4+nb
  nbrec  = nb
  if (nbsend .gt. 8) nbsend = nbsend - 8
  call MPI_SENDRECV(pmax_nb(0),     1, MPI_INTEGER, neighbor(nbsend), 1, &
                    pmax_nb(nbrec), 1, MPI_INTEGER, neighbor(nbrec),  1, &
                    comm_cart, status, ierr)
  call MPI_SENDRECV(mslv_nb(1,0),     pmax,           MPI_INTEGER, neighbor(nbsend), 2, &
                    mslv_nb(1,nbrec), pmax_nb(nbrec), MPI_INTEGER, neighbor(nbrec),  2, &
                    comm_cart, status, ierr)
enddo



dt_step = dt

if (istep==1) then
  do pp=1,pmax                      
    if (ap(pp)%mslv>0) then          

      ap(pp)%xfpo(:)=ap(pp)%xll(:)  
      ap(pp)%yfpo(:)=ap(pp)%yll(:)
      ap(pp)%zfpo(:)=ap(pp)%zll(:)
      ap(pp)%qno1(:)=ap(pp)%qn1(:)
      ap(pp)%qno2(:)=ap(pp)%qn2(:)
      ap(pp)%qno3(:)=ap(pp)%qn3(:)
      ap(pp)%qno4(:)=ap(pp)%qn4(:)

      ap(pp)%xfpold(:)=ap(pp)%xll(:)  
      ap(pp)%yfpold(:)=ap(pp)%yll(:)
      ap(pp)%zfpold(:)=ap(pp)%zll(:)
      ap(pp)%qnold1(:)=ap(pp)%qn1(:)
      ap(pp)%qnold2(:)=ap(pp)%qn2(:)
      ap(pp)%qnold3(:)=ap(pp)%qn3(:)
      ap(pp)%qnold4(:)=ap(pp)%qn4(:)

      ap(pp)%dxdto(:) = ap(pp)%dxdtl(:)  
      ap(pp)%dydto(:) = ap(pp)%dydtl(:)
      ap(pp)%dzdto(:) = ap(pp)%dzdtl(:)
      ap(pp)%omgo1(:) = ap(pp)%omg1(:)
      ap(pp)%omgo2(:) = ap(pp)%omg2(:)
      ap(pp)%omgo3(:) = ap(pp)%omg3(:)
        
      ap(pp)%uao(:)    = ap(pp)%ua(:) 
      ap(pp)%vao(:)    = ap(pp)%va(:)
      ap(pp)%wao(:)    = ap(pp)%wa(:)
      ap(pp)%omgdo1(:) = ap(pp)%omgd1(:)
      ap(pp)%omgdo2(:) = ap(pp)%omgd2(:)
      ap(pp)%omgdo3(:) = ap(pp)%omgd3(:)      

    endif
  enddo
endif



!===========================================================v MAIN LOOP v==========================================================!
do lcounter=1,filsub      

	call profiler_start("IBM_O", tag = .true., tag_color = COLOR_RED)
	call ext_force(unew,vnew,wnew,dt, istep, time, forcex, forcey, forcez)
	call profiler_stop("IBM_O")


	call profiler_start("FIBER_T", tag = .true., tag_color = COLOR_RED)
  do pp=1,pmax
    if (ap(pp)%mslv.gt.0) then

 		  !--- Initialize l_step, nr_step and eN for each fiber 
			l_step  = istep !l_step + 1 
 		  nr_step = 0
 		  eN      = 1.0D+00
      DispX   = 0.0D+00
			do icp =1,nel*ngps
				do jcp=1,3
					strain_kp0(icp,jcp) = strain_k0(icp,jcp,pp)
				enddo
			enddo

      ap(pp)%xll(:)=ap(pp)%xfpo(:)     
      ap(pp)%yll(:)=ap(pp)%yfpo(:)
      ap(pp)%zll(:)=ap(pp)%zfpo(:)
      ap(pp)%qn1(:)=ap(pp)%qno1(:)
      ap(pp)%qn2(:)=ap(pp)%qno2(:)
      ap(pp)%qn3(:)=ap(pp)%qno3(:)
      ap(pp)%qn4(:)=ap(pp)%qno4(:)

      ap(pp)%dxdtl(:) = ap(pp)%dxdto(:) * (1.0 - (Ngamma/Nbeta)) &  
        + ap(pp)%uao(:) * dt_step * (1.0 - (Ngamma/(2.0*Nbeta)))
      ap(pp)%dydtl(:) = ap(pp)%dydto(:) * (1.0 - (Ngamma/Nbeta)) &
        + ap(pp)%vao(:) * dt_step * (1.0 - (Ngamma/(2.0*Nbeta)))
      ap(pp)%dzdtl(:) = ap(pp)%dzdto(:) * (1.0 - (Ngamma/Nbeta)) &
        + ap(pp)%wao(:) * dt_step * (1.0 - (Ngamma/(2.0*Nbeta)))
      ap(pp)%omg1(:) = ap(pp)%omgo1(:) * (1.0 - (Ngamma/Nbeta)) &
        + ap(pp)%omgdo1(:) * dt_step * (1.0 - (Ngamma/(2.0*Nbeta)))
      ap(pp)%omg2(:) = ap(pp)%omgo2(:) * (1.0 - (Ngamma/Nbeta)) &
        + ap(pp)%omgdo2(:) * dt_step * (1.0 - (Ngamma/(2.0*Nbeta)))
      ap(pp)%omg3(:) = ap(pp)%omgo3(:) * (1.0 - (Ngamma/Nbeta)) &
        + ap(pp)%omgdo3(:) * dt_step * (1.0 - (Ngamma/(2.0*Nbeta)))
      
      ap(pp)%ua(:) = -ap(pp)%dxdto(:) * (1.0/(Nbeta*dt_step)) & 
 	      - ap(pp)%uao(:) * ((1.0/(2.0*Nbeta))-1.0)
      ap(pp)%va(:) = -ap(pp)%dydto(:) * (1.0/(Nbeta*dt_step)) &
 	      - ap(pp)%vao(:) * ((1.0/(2.0*Nbeta))-1.0)
      ap(pp)%wa(:) = -ap(pp)%dzdto(:) * (1.0/(Nbeta*dt_step)) &
 	      - ap(pp)%wao(:) * ((1.0/(2.0*Nbeta))-1.0)
      ap(pp)%omgd1(:) = -ap(pp)%omgo1(:) * (1.0/(Nbeta*dt_step)) &
 	      - ap(pp)%omgdo1(:) * ((1.0/(2.0*Nbeta))-1.0)
      ap(pp)%omgd2(:) = -ap(pp)%omgo2(:) * (1.0/(Nbeta*dt_step)) &
 	      - ap(pp)%omgdo2(:) * ((1.0/(2.0*Nbeta))-1.0)  
      ap(pp)%omgd3(:) = -ap(pp)%omgo3(:) * (1.0/(Nbeta*dt_step)) &
 	      - ap(pp)%omgdo3(:) * ((1.0/(2.0*Nbeta))-1.0)
       

      !---- NR loop
      do while (eN .GT. TOL)  

        nr_step = nr_step + 1

        ! Initialize: FG, KG, MG
        KG = 0.0D+00
        FG = 0.0D+00

        ! Element loop
        do ii = 1,nel
          
          ! connectivity arrays          
          cone   = CONS(ii,1:p1_fordr)
          ie_con = (/ 6*cone, 6*cone-5, 6*cone-4, 6*cone-3, & 
						          6*cone-2, 6*cone-1 /)
          call sort_pick(6*p1_fordr, ie_con)
          
          ! Initialize: disp+rot and vel variables
          do icp = 1,p1_fordr
            state_x(icp,1) = ap(pp)%xll(cone(icp))
            state_x(icp,2) = ap(pp)%yll(cone(icp))
            state_x(icp,3) = ap(pp)%zll(cone(icp))
            state_x(icp,4) = ap(pp)%qn1(cone(icp))
            state_x(icp,5) = ap(pp)%qn2(cone(icp))
            state_x(icp,6) = ap(pp)%qn3(cone(icp))
            state_x(icp,7) = ap(pp)%qn4(cone(icp))            
            
            state_w(icp,1) = ap(pp)%dxdtl(cone(icp))
            state_w(icp,2) = ap(pp)%dydtl(cone(icp))
            state_w(icp,3) = ap(pp)%dzdtl(cone(icp))
            state_w(icp,4) = ap(pp)%omg1(cone(icp))
            state_w(icp,5) = ap(pp)%omg2(cone(icp))
            state_w(icp,6) = ap(pp)%omg3(cone(icp))

            uae(icp,1) = ap(pp)%ua(cone(icp))
            uae(icp,2) = ap(pp)%va(cone(icp))
            uae(icp,3) = ap(pp)%wa(cone(icp))
            uae(icp,4) = ap(pp)%omgd1(cone(icp))
            uae(icp,5) = ap(pp)%omgd2(cone(icp))
            uae(icp,6) = ap(pp)%omgd3(cone(icp))

						fd_CP(6*(icp-1)+1) = ap(pp)%fxll(cone(icp))
						fd_CP(6*(icp-1)+2) = ap(pp)%fyll(cone(icp))
						fd_CP(6*(icp-1)+3) = ap(pp)%fzll(cone(icp))
						fd_CP(6*(icp-1)+4) = ap(pp)%mxl(cone(icp))
						fd_CP(6*(icp-1)+5) = ap(pp)%myl(cone(icp))
						fd_CP(6*(icp-1)+6) = ap(pp)%mzl(cone(icp))
          end do

          ! Compute GPs in pametric space
          elU = elRangeU(ii,1:2)
          c1  = 0.5 * (elU(2) - elU(1))
          c2  = 0.5 * (elU(2) + elU(1))
										
					! Compute basis function array for ngps,ngpb
					inx = (ii-1)*ngps*p1_fordr 
					do jj=1,ngps*p1_fordr
						NdN_temp(jj,:)=NdN_array(inx+jj,:)
					enddo
					!
					inx = (ii-1)*ngpb*p1_fordr 
          do jj=1,ngpb*p1_fordr
            NdNb_temp(jj,:)=NdN_arrayb(inx+jj,:)
          enddo
					
					! stain in material frame
					strain_gp0 = strain_kp0((ii-1)*ngps+1:ii*ngps,1:3)

          !---- Compute internal force, Q0
          Q0 = 0.0D+00


					do icp=1,ngps
            Jp_s(icp) = J_s((nel-1)*ngps+icp, pp)
          enddo

					call compute_fint( state_x, state_w, elU, ngps, &
          xigs, wgs, Jp_s, NdN_temp, p1_fordr, nno, uKnot, &
          I_zz, I_yy, J_xx, Ks_y, Ks_z, E_mod, G_mod, Arf, &
          nderiv, Q0, Eflag, E_modA, G_modA, EI_yy, EI_zz, &
					GJxx, strain_gp0, istep)
					

          kele        = 0.0D+00
          state_x_inc = 0.0D+00
          state_delta = 0.0D+00          
					!
					do icp=1,ngps
            Jp_s(icp) = J_s((nel-1)*ngps+icp, pp)
          enddo
          do jj = 1, rows_w
            state_delta(jj) = state_delta(jj) + Deltaf
            call stateIncrement(state_x, state_delta, &
                                  p1_fordr, state_x_inc)	 				    	

            Q1 = 0.0D+00
            call compute_fint( state_x_inc, state_w, elRangeU(ii,:), ngps, &
                  xigs, wgs, Jp_s, NdN_temp, p1_fordr, nno, uKnot, &
                  I_zz, I_yy, J_xx, Ks_y, Ks_z, E_mod, G_mod, Arf, &
                  nderiv, Q1, Eflag, E_modA, G_modA, EI_yy, EI_zz, &
									GJxx, strain_gp0, istep)

            ! Compute stifness matrix
            kele(:,jj)   = (Q1 - Q0) / Deltaf
            state_delta = state_delta(jj) - Deltaf
          enddo ! rows_w loop for kele
					!
          
					!--- Compute inertia terms: mele and finert                
          mele = 0.0D+00
					do icp=1,ngpb
            Jp_b(icp) = J_b((nel-1)*ngpb+icp, pp)
          enddo
          call compute_massE( mInerCP, elRangeU(ii,1:2), state_x, &
          nderiv, ngpb, xigb, wgb, Jp_b, NdNb_temp, p1_fordr, nno, &
					uKnot, mele )
					!

					!
          finert = 0.0D+00
          call compute_finer( mInerCP, elRangeU(ii,1:2), &
          state_x, state_w, nderiv, ngpb, xigb, wgb, rhof, &
          I_yy, I_zz, J_xx, uae, NdNb_temp, p1_fordr, nno, &
					uKnot, finert, Eflag, rhoA, rhoI_yy, rhoI_zz )
					!

          ! Incorporate inerta terms
          kele = kele + (1.0/(Nbeta*dt_step**2)) * mele
          Q0   = Q0   + finert

          !---- Compute force due to body force
          if (inertia .EQ. 1) then
            mele = 0.0D+00
					  do icp=1,ngpb
              Jp_b(icp) = J_b((nel-1)*ngpb+icp, pp)
            enddo
            call compute_massE( mInerCP2, elRangeU(ii,1:2), state_x, &
              nderiv, ngpb, xigb, wgb, Jp_b, &
					  	NdNb_temp, p1_fordr, nno, uKnot, mele )

					  fgrv = 0.0D+00
            do icp = 1,p1_fordr 
                do jcp = 1,3 !p1_fordr
                    iroco(jcp) = 6*(icp-1)+jcp
                end do  !jcp loop
                fgrv(iroco) = matmul(mele(iroco, iroco), g_vec)
            end do ! icp loop            
            Q0 = Q0 - fgrv

          endif ! gflag

          fdise = 0.0D+00
          call compute_fdise( fd_CP, state_x, elRangeU(ii,1:2), &
					nderiv, ngpb, xigb, wgb, NdNb_temp, p1_fordr, nno, uKnot, &
					Arf, fdise )
          Q0 = Q0 + fdise

          !--- Global Assembly
          FG(ie_con)         = FG(ie_con) + Q0
          KG(ie_con, ie_con) = KG(ie_con, ie_con) + kele

        end do ! nel loop  

        ! Apply boundary conditions
        KGr = KG(ir_dof, ir_dof)
        FGr = FG(ir_dof)

        !--- Solve
        !
        call DGETRF(nir,nir,KGr,nir,IPIV,info)        
        call DGETRI(nir,KGr,nir,IPIV,WORK,nir,info)
        DUr           = -1 * matmul(KGr,FGr)
				DispU         = 0.0D+00 
        DispU(ir_dof) = DUr
        DispX         = DispX + DispU
        !


        do ii = 1, nno

          vec3           = (/ (jj, jj=(ii-1)*6+1, ii*6-3) /)          
          ap(pp)%xll(ii) = ap(pp)%xll(ii) + DispU(vec3(1))  
          ap(pp)%yll(ii) = ap(pp)%yll(ii) + DispU(vec3(2))
          ap(pp)%zll(ii) = ap(pp)%zll(ii) + DispU(vec3(3))          
          q_old          = (/ ap(pp)%qN1(ii), ap(pp)%qN2(ii), &
                             ap(pp)%qN3(ii), ap(pp)%qN4(ii)/)
          vec3           = (/ (jj, jj = (ii-1)*6+4, ii*6) /)
          Drot_vec       = DispU(vec3)
          call q_from_Rotv(Drot_vec, rel_q)
          call quad_prod(q_old, rel_q, q_current)          
          ap(pp)%qN1(ii) = q_current(1)
          ap(pp)%qN2(ii) = q_current(2)
          ap(pp)%qN3(ii) = q_current(3)
          ap(pp)%qN4(ii) = q_current(4)

          vec6 = (/ (jj, jj = (ii-1)*6+1, ii*6) /)
          ap(pp)%dxdtl(ii) = ap(pp)%dxdto(ii) * (1.0 - (Ngamma/Nbeta)) & 
             + ap(pp)%uao(ii) * dt_step * (1.0 - (Ngamma/(2.0*Nbeta)))
          ap(pp)%dydtl(ii) = ap(pp)%dydto(ii) * (1.0 - (Ngamma/Nbeta)) &
            + ap(pp)%vao(ii) * dt_step * (1.0 - (Ngamma/(2.0*Nbeta)))  
          ap(pp)%dzdtl(ii) = ap(pp)%dzdto(ii) * (1.0 - (Ngamma/Nbeta)) &
            + ap(pp)%wao(ii) * dt_step * (1.0 - (Ngamma/(2.0*Nbeta)))
          ap(pp)%omg1(ii) = ap(pp)%omgo1(ii) * (1.0 - (Ngamma/Nbeta)) &
            + ap(pp)%omgdo1(ii) * dt_step * (1.0 - (Ngamma/(2.0*Nbeta)))
          ap(pp)%omg2(ii) = ap(pp)%omgo2(ii) * (1.0 - (Ngamma/Nbeta)) &
            + ap(pp)%omgdo2(ii) * dt_step * (1.0 - (Ngamma/(2.0*Nbeta)))
          ap(pp)%omg3(ii) = ap(pp)%omgo3(ii) * (1.0 - (Ngamma/Nbeta)) &
            + ap(pp)%omgdo3(ii) * dt_step * (1.0 - (Ngamma/(2.0*Nbeta)))
          
          ap(pp)%dxdtl(ii) = ap(pp)%dxdtl(ii) + (Ngamma/(Nbeta*dt_step)) &
            * DispX(vec6(1))
          ap(pp)%dydtl(ii) = ap(pp)%dydtl(ii) + (Ngamma/(Nbeta*dt_step)) &
            * DispX(vec6(2))
          ap(pp)%dzdtl(ii) = ap(pp)%dzdtl(ii) + (Ngamma/(Nbeta*dt_step)) &
            * DispX(vec6(3))  
          ap(pp)%omg1(ii) = ap(pp)%omg1(ii) + (Ngamma/(Nbeta*dt_step)) &
            * DispX(vec6(4))  
          ap(pp)%omg2(ii) = ap(pp)%omg2(ii) + (Ngamma/(Nbeta*dt_step)) &
            * DispX(vec6(5))
          ap(pp)%omg3(ii) = ap(pp)%omg3(ii) + (Ngamma/(Nbeta*dt_step)) &
            * DispX(vec6(6))
                    

          ap(pp)%ua(ii) = -ap(pp)%dxdto(ii) * (1.0/(Nbeta*dt_step)) &
            - ap(pp)%uao(ii) * ((1.0/(2.0*Nbeta))-1.0)
          ap(pp)%va(ii) = -ap(pp)%dydto(ii) * (1.0/(Nbeta*dt_step)) &
            - ap(pp)%vao(ii) * ((1.0/(2.0*Nbeta))-1.0)  
          ap(pp)%wa(ii) = -ap(pp)%dzdto(ii) * (1.0/(Nbeta*dt_step)) &
            - ap(pp)%wao(ii) * ((1.0/(2.0*Nbeta))-1.0)
          ap(pp)%omgd1(ii) = -ap(pp)%omgo1(ii) * (1.0/(Nbeta*dt_step)) &
            - ap(pp)%omgdo1(ii) * ((1.0/(2.0*Nbeta))-1.0)            
          ap(pp)%omgd2(ii) = -ap(pp)%omgo2(ii) * (1.0/(Nbeta*dt_step)) &
            - ap(pp)%omgdo2(ii) * ((1.0/(2.0*Nbeta))-1.0)
          ap(pp)%omgd3(ii) = -ap(pp)%omgo3(ii) * (1.0/(Nbeta*dt_step)) &
            - ap(pp)%omgdo3(ii) * ((1.0/(2.0*Nbeta))-1.0)
            

          ap(pp)%ua(ii) = ap(pp)%ua(ii) + (1.0/(Nbeta*(dt_step**2))) &
            * DispX(vec6(1))
          ap(pp)%va(ii) = ap(pp)%va(ii) + (1.0/(Nbeta*(dt_step**2))) &
            * DispX(vec6(2))
          ap(pp)%wa(ii) = ap(pp)%wa(ii) + (1.0/(Nbeta*(dt_step**2))) &
            * DispX(vec6(3))
          ap(pp)%omgd1(ii) = ap(pp)%omgd1(ii) + (1.0/(Nbeta*(dt_step**2))) &
            * DispX(vec6(4))
          ap(pp)%omgd2(ii) = ap(pp)%omgd2(ii) + (1.0/(Nbeta*(dt_step**2))) &
            * DispX(vec6(5))
          ap(pp)%omgd3(ii) = ap(pp)%omgd3(ii) + (1.0/(Nbeta*(dt_step**2))) &
            * DispX(vec6(6))
        end do  ! end of ii = 1, nno

        ! Compute eN
        eN = norm2(FGr)
        !if (mod(l_step, 500) .eq. 0 .or. nr_step .ge. 50) then
        !        write(*,*) "FResidual at NR_itr ", nr_step, 'is:  ', REAL(log10(eN),4)
        !endif
        
				if (nr_step .ge. 50) then
					print*,"__________nr_step > 50________"
					call mpi_abort(comm_cart,ierr,ierr) !stop
				endif

      end do ! NR loop

      !if ( mod(l_step, 1000) .eq. 0 .and. pp .eq. np) then 
      !  write(*,*) "Load step", l_step, "used", nr_step, "Newto-steps"
      !  !print*, "-----------------------------------------"
      !endif

      ! Update: Disp, vel, acc
      ap(pp)%xfpo(:)=ap(pp)%xll(:)  
      ap(pp)%yfpo(:)=ap(pp)%yll(:)
      ap(pp)%zfpo(:)=ap(pp)%zll(:)
      ap(pp)%qno1(:)=ap(pp)%qn1(:)
      ap(pp)%qno2(:)=ap(pp)%qn2(:)
      ap(pp)%qno3(:)=ap(pp)%qn3(:)
      ap(pp)%qno4(:)=ap(pp)%qn4(:)

      ap(pp)%dxdto(:)=ap(pp)%dxdtl(:)  
      ap(pp)%dydto(:)=ap(pp)%dydtl(:)
      ap(pp)%dzdto(:)=ap(pp)%dzdtl(:)
      ap(pp)%omgo1(:)=ap(pp)%omg1(:)
      ap(pp)%omgo2(:)=ap(pp)%omg2(:)
      ap(pp)%omgo3(:)=ap(pp)%omg3(:)

      ap(pp)%uao(:)   =ap(pp)%ua(:)  
      ap(pp)%vao(:)   =ap(pp)%va(:)  
      ap(pp)%wao(:)   =ap(pp)%wa(:)
      ap(pp)%omgdo1(:)=ap(pp)%omgd1(:)
      ap(pp)%omgdo2(:)=ap(pp)%omgd2(:)
      ap(pp)%omgdo3(:)=ap(pp)%omgd3(:)

      
			!*****************************************************updating velocities**************************************************


			ur=(real(lcounter)/((filsub/2)+1))
			if (ur>1.) ur=1.
			ur=1.-ur
			!ur=1.-(real(lcounter)/filsub)

    endif   ! conditional statement if ap(pp)%mslv>0 
  enddo     ! loop over p=1,pmax

  call profiler_stop("FIBER_T")
enddo       ! loop over lcounter = 1, filesub


!v below, ua is acceleratin of nodes
do pp=1,pmax

  if (ap(pp)%mslv.gt.0) then


    !********************************************************updating old values***********************************************

    ap(pp)%xfpold(:) = ap(pp)%xfpo(:)     !v values at t_(n-2)
    ap(pp)%yfpold(:) = ap(pp)%yfpo(:)
    ap(pp)%zfpold(:) = ap(pp)%zfpo(:) 
    ap(pp)%qNold1(:) = ap(pp)%qNo1(:)     !v values at t_(n-2)
    ap(pp)%qNold2(:) = ap(pp)%qNo2(:)     !v values at t_(n-2)
    ap(pp)%qNold3(:) = ap(pp)%qNo3(:)     !v values at t_(n-2)
    ap(pp)%qNold4(:) = ap(pp)%qNo4(:)     !v values at t_(n-2)
    

    !********************************************computing length of filament**************************************************
    !v after 
    !length=0.
    !do l=1,NL-1
    !  length=length+sqrt((((ap(pp)%xfp(l+1)-ap(pp)%xfp(l))**2.) + ((ap(pp)%yfp(l+1)-ap(pp)%yfp(l))**2.)  +  &
    !    ((ap(pp)%zfp(l+1)-ap(pp)%zfp(l))**2.)))
    !enddo
    !print*,"length for particle",ap(pp)%mslv,"=",length
    !print*,"relative length deviation=",abs(length-1.)/1.
    !if (floor(length+.8)/=1. .and. floor(length)/=0.)then
       !print*,&
       !"divergence occured in filament equations...Length=",length,'number',ap(pp)%mslv,'position:', ap(pp)%x,ap(pp)%y,ap(pp)%z
    !  do l=1,nl
    !  enddo
    !  !call mpi_finalize(ierr)
    !  !stop
    !endif

		!********************************************computing length of filament**************************************************
    length=0.
    do l=1,nno-1
      length=length+sqrt((((ap(pp)%xll(l+1)-ap(pp)%xll(l))**2.)+((ap(pp)%yll(l+1)-ap(pp)%yll(l))**2.)+ &
      ((ap(pp)%zll(l+1)-ap(pp)%zll(l))**2.)))
    enddo
    !print*,"length for particle",ap(pp)%mslv,"=",length
    !print*,"relative length deviation=",abs(length-L_fibr)/1.
		!print*,"-----------------------------------------------------------------------------"


    !*******************************************************tensionplot********************************************************

    if (mod(istep,100)==0.and.ap(pp)%mslv==1) then
      open(102,file=datadir//'pos.txt',position='append')
      write(102,'(3E15.7)')  time,abs(atan((ap(pp)%zll(nno)-ap(pp)%z)/(ap(pp)%yll(nno)-ap(pp)%y))/pi), &
        sqrt((ap(pp)%dydtl(nno)**2)+(ap(pp)%dzdtl(nno)**2))*2.
      close(102)
    endif

    ap(pp)%x = ap(pp)%xll(int((nno+1)/2.))
    ap(pp)%y = ap(pp)%yll(int((nno+1)/2.)) 
    ap(pp)%z = ap(pp)%zll(int((nno+1)/2.))

    if (ap(pp)%x .gt. lx) then
			ap(pp)%x = ap(pp)%x - lx
			do l=1,nno
        ap(pp)%xll(l) = ap(pp)%xll(l) - lx
        ap(pp)%xfpo(l) = ap(pp)%xfpo(l) - lx
        ap(pp)%xfpold(l) = ap(pp)%xfpold(l) - lx
      enddo
    endif

    if (ap(pp)%x .lt. 0.) then
      ap(pp)%x = ap(pp)%x + lx
      do l=1,nno
        ap(pp)%xll(l) = ap(pp)%xll(l) + lx
        ap(pp)%xfpo(l) = ap(pp)%xfpo(l) + lx
        ap(pp)%xfpold(l) = ap(pp)%xfpold(l) + lx
      enddo
    endif

    if (ap(pp)%y .gt. ly) then            
      ap(pp)%y = ap(pp)%y - ly
      do l=1,nno
        ap(pp)%yll(l)    = ap(pp)%yll(l) - ly
        ap(pp)%yfpo(l)   = ap(pp)%yfpo(l) - ly
        ap(pp)%yfpold(l) = ap(pp)%yfpold(l) - ly
      enddo
    endif

    if (ap(pp)%y .lt. 0.) then
      ap(pp)%y = ap(pp)%y + ly
      do l=1,nno
        ap(pp)%yll(l)    = ap(pp)%yll(l) + ly
        ap(pp)%yfpo(l)   = ap(pp)%yfpo(l) + ly
        ap(pp)%yfpold(l) = ap(pp)%yfpold(l) + ly
      enddo
    endif

		if (ap(pp)%z .gt. lz) then            
      ap(pp)%z = ap(pp)%z - lz
      do l=1,nno
        ap(pp)%zll(l)    = ap(pp)%zll(l) - lz
        ap(pp)%zfpo(l)   = ap(pp)%zfpo(l) - lz
        ap(pp)%zfpold(l) = ap(pp)%zfpold(l) - lz
      enddo
    endif

    if (ap(pp)%z .lt. 0.) then
      ap(pp)%z = ap(pp)%z + lz
      do l=1,nno
        ap(pp)%zll(l)    = ap(pp)%zll(l) + lz
        ap(pp)%zfpo(l)   = ap(pp)%zfpo(l) + lz
        ap(pp)%zfpold(l) = ap(pp)%zfpold(l) + lz
      enddo
    endif

    !**************************************************************************************************************************

  endif
enddo



nt =0
ntp=0
nbb=0
nbp=0

call mpi_allreduce(ntp, nt,  1, mpi_integer, mpi_sum, comm_cart, ierr)
call mpi_allreduce(nbp, nbb, 1, mpi_integer, mpi_sum, comm_cart, ierr)



! particle positions were updated. Now check if there are new masters
!
!
!
do pp=1,pmax                      !v for all master-slave parts
  newmaster_nb(pp,0) = 0          !v Initialized newmasters with null
  if (ap(pp)%mslv.gt.0) then      !v for all master part
    ! myid was master of particle ap(pp)%mslv at previous time step n
    ax = 0.5
    ay = 0.5
    if (ap(pp)%x.eq.lx) ax = 0.51
    if (ap(pp)%x.eq.0)  ax = 0.49
    if (ap(pp)%y.eq.ly)  ay = 0.51
    if (ap(pp)%y.eq.0)   ay = 0.49
    proccoords(1) = nint( dims(1)*ap(pp)%x/lx - ax )
    proccoords(2) = nint( dims(2)*ap(pp)%y/ly - ay )
    call MPI_CART_RANK(comm_cart, proccoords, procrank, ierr)  !v Input: comm_cart, procoords. Outupt: procrank, ierr
                !v MP_CART_RANK:  Determines process rank in communicator given cartesian location
                !v comm_art: communicator with Cartesian structure,
                !v proccorrds: Integer array (of size ndims) specifying Cartesian coords of a process
                !v procrank: rank of specified process
    !v if (procrank .ne. myid  .and.  ap(pp)%x>=0.  .and.  ap(pp)%x<=lx) then
    if (procrank .ne. myid) then
      ! particle ap(pp)%mslv has a new master at time step n+1
      do nb=1,8
        if (procrank .eq. neighbor(nb)) then
          newmaster_nb(pp,0) = nb
        endif
      enddo
    endif
  endif
enddo






!!!!!!!!!!!!!!!!!!!!!!!!!!MPI_integer should be changed!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! exchange data
!

do nb=1,8
  nbsend = nb
  nbrec  = nb+4
  if (nbrec .gt. 8) nbrec = nbrec-8
  call MPI_SENDRECV(newmaster_nb(1,0),    pmax,           MPI_INTEGER, neighbor(nbsend), 1, &
                   newmaster_nb(1,nbrec), pmax_nb(nbrec), MPI_INTEGER, neighbor(nbrec),  1, &
                   comm_cart, status,ierr)
  !v MPI_SendRecv(): send and receive simultaneouly, to shift operation across a chain of processes, mpi
                                                              ! interger in this case (as mentioned above)
enddo
!



do pp=1,pmax
  nrrequests = 0
  if (newmaster_nb(pp,0).gt.0) then
    nbsend = newmaster_nb(pp,0)
    !v tag = ap(pp)%mslv*10+nbsend
    tag = ap(pp)%mslv!*10+nbsend
    nrrequests = nrrequests + 1
    call MPI_ISEND(ap(pp)%x,send_real,MPI_REAL8,neighbor(nbsend),tag,comm_cart,arrayrequests((nrrequests-1)*2+1),ierr)
    !call MPI_ISEND(ap(pp)%x,send_real,MPI_REAL8,neighbor(nbsend),tag,comm_cart,arrayrequests(nrrequests),ierr)
    ap(pp)%mslv = -ap(pp)%mslv ! master is a slave now
  endif

  if(mslv_nb(pp,0).lt.0) then
    idp = -ap(pp)%mslv
    do nbrec = 1,8
      nbrec2 = nbrec + 4
      if(nbrec2 .gt. 8) nbrec2 = nbrec2-8
      do i=1,pmax_nb(nbrec)
        idp_nb = mslv_nb(i,nbrec)
        if(newmaster_nb(i,nbrec) .eq. nbrec2.and.idp.eq.idp_nb) then
          ap(pp)%mslv = -ap(pp)%mslv ! slave became a master
          nrrequests = nrrequests + 1
          !v tag = ap(pp)%mslv*10+nbrec2
          tag = ap(pp)%mslv!*10+nbrec2
          call MPI_IRECV(ap(pp)%x,send_real,MPI_REAL8,neighbor(nbrec),tag,comm_cart,arrayrequests((nrrequests-1)*2+1),ierr)
          !call MPI_IRECV(ap(pp)%x,send_real,MPI_REAL8,neighbor(nbrec),tag,comm_cart,arrayrequests(nrrequests),ierr)
        endif
      enddo
    enddo
  endif
  nrrequests = nrrequests
  call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,ierr)
enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Masters are known now: ap(pp)%mslv > 0.
! Next step: determine slaves and master/slave neighbors
!
forall(p=1:npmax,k=0:8) mslv_nb(p,k) = 0
mslv_nb(1:pmax,0) = ap(1:pmax)%mslv  ! new masters, but new slaves not all determined yet!
! process might remain master, but slave might lose particle
anb2(1:pmax,0)%x  = ap(1:pmax)%x  ! from new master
anb2(1:pmax,0)%y  = ap(1:pmax)%y  ! from new master
anb2(1:pmax,0)%z  = ap(1:pmax)%z  ! from new master

do l=1,nno
  anb2(1:pmax,0)%xll(l) = ap(1:pmax)%xll(l) ! from new master
  anb2(1:pmax,0)%yll(l) = ap(1:pmax)%yll(l) ! from new master
  anb2(1:pmax,0)%zll(l) = ap(1:pmax)%zll(l) ! from new master
enddo

do nb=1,8
  nbsend = 4+nb
  if (nbsend .gt. 8) nbsend = nbsend - 8
  nbrec  = nb
  call MPI_SENDRECV(mslv_nb(1,0),pmax,MPI_INTEGER,neighbor(nbsend),1, &
    mslv_nb(1,nbrec),pmax_nb(nbrec),MPI_INTEGER,neighbor(nbrec),1, &
    comm_cart,status,ierr)
  call MPI_SENDRECV(anb2(1,0)%x,pmax*3*(1+nno),MPI_REAL8,neighbor(nbsend),2, &
    anb2(1,nbrec)%x,pmax_nb(nbrec)*3*(1+nno),MPI_REAL8,neighbor(nbrec),2, &
    comm_cart,status,ierr)
  !send x,y,z -> 3*pmax contiguous info
  !(see definition of type neighbor2 in the begining of the subroutine)
enddo
!
! recompute particle positions because of periodic b.c.'s in the x and y-direction
!
nb=1
do p=1,pmax_nb(nb)
  if (mslv_nb(p,nb) .gt. 0) then
    boundleftnb = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
    if (anb2(p,nb)%x .lt. boundleftnb) then
			anb2(p,nb)%x = anb2(p,nb)%x + lx
      do l=1,nno
        anb2(p,nb)%xll(l) = anb2(p,nb)%xll(l) + lx
      enddo
    endif
  endif
enddo
nb=2
do p=1,pmax_nb(nb)
    if (mslv_nb(p,nb) .gt. 0) then
      boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
      boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
      if (anb2(p,nb)%x .lt. boundleftnb) then
        anb2(p,nb)%x = anb2(p,nb)%x + lx
        do l=1,nno
          anb2(p,nb)%xll(l) = anb2(p,nb)%xll(l) + lx
        enddo
      endif
      if (anb2(p,nb)%y .gt. boundbacknb) then
        anb2(p,nb)%y = anb2(p,nb)%y - ly
        do l=1,nno
          anb2(p,nb)%yll(l) = anb2(p,nb)%yll(l) - ly
        enddo
      endif
    endif
enddo
nb=3
do p=1,pmax_nb(nb)
  if (mslv_nb(p,nb) .gt. 0) then
    boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
    if (anb2(p,nb)%y .gt. boundbacknb) then
      anb2(p,nb)%y = anb2(p,nb)%y - ly
      do l=1,nno
        anb2(p,nb)%yll(l) = anb2(p,nb)%yll(l) - ly
      enddo
    endif
  endif
enddo
nb=4
do p=1,pmax_nb(nb)
  if (mslv_nb(p,nb) .gt. 0) then
    boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
    boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
    if (anb2(p,nb)%x .gt. boundrightnb) then
       anb2(p,nb)%x = anb2(p,nb)%x - lx
       do l=1,nno
         anb2(p,nb)%xll(l) = anb2(p,nb)%xll(l) - lx
       enddo
    endif
    if (anb2(p,nb)%y .gt. boundbacknb) then
      anb2(p,nb)%y = anb2(p,nb)%y - ly
      do l=1,nno
        anb2(p,nb)%yll(l) = anb2(p,nb)%yll(l) - ly
      enddo
    endif
  endif
enddo
nb=5
do p=1,pmax_nb(nb)
  if (mslv_nb(p,nb) .gt. 0) then
    boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
    if (anb2(p,nb)%x .gt. boundrightnb) then
       anb2(p,nb)%x = anb2(p,nb)%x - lx
       do l=1,nno
         anb2(p,nb)%xll(l) = anb2(p,nb)%xll(l) - lx
       enddo
    endif
  endif
enddo
nb=6
do p=1,pmax_nb(nb)
  if (mslv_nb(p,nb) .gt. 0) then
    boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
    boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
    if (anb2(p,nb)%x .gt. boundrightnb) then
       anb2(p,nb)%x = anb2(p,nb)%x - lx
       do l=1,nno
         anb2(p,nb)%xll(l) = anb2(p,nb)%xll(l) - lx
       enddo
    endif

    if (anb2(p,nb)%y .lt. boundfrontnb) then
      anb2(p,nb)%y = anb2(p,nb)%y + ly
      do l=1,nno
        anb2(p,nb)%yll(l) = anb2(p,nb)%yll(l) + ly
      enddo
    endif
  endif
enddo
nb=7
do p=1,pmax_nb(nb)
  if (mslv_nb(p,nb) .gt. 0) then
    boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
    if (anb2(p,nb)%y .lt. boundfrontnb) then
      anb2(p,nb)%y = anb2(p,nb)%y + ly
      do l=1,nno
        anb2(p,nb)%yll(l) = anb2(p,nb)%yll(l) + ly
      enddo
    endif
  endif
enddo
nb=8
do p=1,pmax_nb(nb)
  if (mslv_nb(p,nb) .gt. 0) then
    boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
    boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
    if (anb2(p,nb)%x .lt. boundleftnb) then
      anb2(p,nb)%x = anb2(p,nb)%x + lx
      do l=1,nno
        anb2(p,nb)%xll(l) = anb2(p,nb)%xll(l) + lx
      enddo
    endif

    if (anb2(p,nb)%y .lt. boundfrontnb) then
      anb2(p,nb)%y = anb2(p,nb)%y + ly
      do l=1,nno
        anb2(p,nb)%yll(l) = anb2(p,nb)%yll(l) + ly
      enddo
    endif
  endif
enddo




!
! important info that should be kept:
!
sp(1:pmax)%mslv      = ap(1:pmax)%mslv
sp(1:pmax)%x         = ap(1:pmax)%x
sp(1:pmax)%y         = ap(1:pmax)%y
sp(1:pmax)%z         = ap(1:pmax)%z
sp(1:pmax)%integralx = ap(1:pmax)%integralx
sp(1:pmax)%integraly = ap(1:pmax)%integraly
sp(1:pmax)%integralz = ap(1:pmax)%integralz


forall(l=1:nno)
  sp(1:pmax)%xll(l)  = ap(1:pmax)%xll(l)
  sp(1:pmax)%yll(l)  = ap(1:pmax)%yll(l)
  sp(1:pmax)%zll(l)  = ap(1:pmax)%zll(l) !v add qNs part here
  sp(1:pmax)%qN1(l)  = ap(1:pmax)%qN1(l) 
  sp(1:pmax)%qN2(l)  = ap(1:pmax)%qN2(l) 
  sp(1:pmax)%qN3(l)  = ap(1:pmax)%qN3(l) 
  sp(1:pmax)%qN4(l)  = ap(1:pmax)%qN4(l) 
  
  sp(1:pmax)%xfpo(l) = ap(1:pmax)%xfpo(l)
  sp(1:pmax)%yfpo(l) = ap(1:pmax)%yfpo(l)
  sp(1:pmax)%zfpo(l) = ap(1:pmax)%zfpo(l)
  sp(1:pmax)%qNo1(l) = ap(1:pmax)%qNo1(l) 
  sp(1:pmax)%qNo2(l) = ap(1:pmax)%qNo2(l) 
  sp(1:pmax)%qNo3(l) = ap(1:pmax)%qNo3(l) 
  sp(1:pmax)%qNo4(l) = ap(1:pmax)%qNo4(l) 

  sp(1:pmax)%xfpold(l) = ap(1:pmax)%xfpold(l)
  sp(1:pmax)%yfpold(l) = ap(1:pmax)%yfpold(l)
  sp(1:pmax)%zfpold(l) = ap(1:pmax)%zfpold(l)
  sp(1:pmax)%qNold1(l) = ap(1:pmax)%qNold1(l) 
  sp(1:pmax)%qNold2(l) = ap(1:pmax)%qNold2(l) 
  sp(1:pmax)%qNold3(l) = ap(1:pmax)%qNold3(l) 
  sp(1:pmax)%qNold4(l) = ap(1:pmax)%qNold4(l) 
  
  sp(1:pmax)%dxdtl(l)= ap(1:pmax)%dxdtl(l)
  sp(1:pmax)%dydtl(l)= ap(1:pmax)%dydtl(l)
  sp(1:pmax)%dzdtl(l)= ap(1:pmax)%dzdtl(l)
  sp(1:pmax)%omg1(l) = ap(1:pmax)%omg1(l)
  sp(1:pmax)%omg2(l) = ap(1:pmax)%omg2(l)
  sp(1:pmax)%omg3(l) = ap(1:pmax)%omg3(l)  
    
  sp(1:pmax)%ua(l)    = ap(1:pmax)%ua(l)
  sp(1:pmax)%va(l)    = ap(1:pmax)%va(l)
  sp(1:pmax)%wa(l)    = ap(1:pmax)%wa(l)
  sp(1:pmax)%omgd1(l) = ap(1:pmax)%omgd1(l)
  sp(1:pmax)%omgd2(l) = ap(1:pmax)%omgd2(l)
  sp(1:pmax)%omgd3(l) = ap(1:pmax)%omgd3(l)

  sp(1:pmax)%fxll(l) = ap(1:pmax)%fxll(l)
  sp(1:pmax)%fyll(l) = ap(1:pmax)%fyll(l)
  sp(1:pmax)%fzll(l) = ap(1:pmax)%fzll(l)
  sp(1:pmax)%mxl(l) = ap(1:pmax)%mxl(l)
  sp(1:pmax)%myl(l) = ap(1:pmax)%myl(l)
  sp(1:pmax)%mzl(l) = ap(1:pmax)%mzl(l)

  sp(1:pmax)%ull(l)    = ap(1:pmax)%ull(l)
  sp(1:pmax)%vll(l)    = ap(1:pmax)%vll(l)
  sp(1:pmax)%wll(l)    = ap(1:pmax)%wll(l)
  sp(1:pmax)%omgl1(l) = ap(1:pmax)%omgl1(l)
  sp(1:pmax)%omgl2(l) = ap(1:pmax)%omgl2(l)
  sp(1:pmax)%omgl3(l) = ap(1:pmax)%omgl3(l)  
  
  sp(1:pmax)%dxdto(l) = ap(1:pmax)%dxdto(l)
  sp(1:pmax)%dydto(l) = ap(1:pmax)%dydto(l)
  sp(1:pmax)%dzdto(l) = ap(1:pmax)%dzdto(l)
  sp(1:pmax)%omgo1(l) = ap(1:pmax)%omgo1(l)
  sp(1:pmax)%omgo2(l) = ap(1:pmax)%omgo2(l)
  sp(1:pmax)%omgo3(l) = ap(1:pmax)%omgo3(l)  

  sp(1:pmax)%uao(l)    = ap(1:pmax)%uao(l)
  sp(1:pmax)%vao(l)    = ap(1:pmax)%vao(l)
  sp(1:pmax)%wao(l)    = ap(1:pmax)%wao(l)
  sp(1:pmax)%omgdo1(l) = ap(1:pmax)%omgdo1(l)
  sp(1:pmax)%omgdo2(l) = ap(1:pmax)%omgdo2(l)
  sp(1:pmax)%omgdo3(l) = ap(1:pmax)%omgdo3(l)
end forall


!
! clear structure ap for re-ordering:
!
ap(1:pmax)%mslv = 0
ap(1:pmax)%x = 0.
ap(1:pmax)%y = 0.
ap(1:pmax)%z = 0.
ap(1:pmax)%integralx = 0.
ap(1:pmax)%integraly = 0.
ap(1:pmax)%integralz = 0.

forall (pp=1:pmax)
  ap(pp)%xll(:)  = 0.
  ap(pp)%yll(:)  = 0.
  ap(pp)%zll(:)  = 0.
  ap(pp)%qN1(:)  = 0.
  ap(pp)%qN2(:)  = 0.
  ap(pp)%qN3(:)  = 0.
  ap(pp)%qN4(:)  = 0.
  
  ap(pp)%xfpo(:) = 0.
  ap(pp)%yfpo(:) = 0.
  ap(pp)%zfpo(:) = 0.
  ap(pp)%qNo1(:)  = 0.
  ap(pp)%qNo2(:)  = 0.
  ap(pp)%qNo3(:)  = 0.
  ap(pp)%qNo4(:)  = 0.

  ap(pp)%xfpold(:) = 0.
  ap(pp)%yfpold(:) = 0.
  ap(pp)%zfpold(:) = 0.
  ap(pp)%qNold1(:) = 0.
  ap(pp)%qNold2(:) = 0.
  ap(pp)%qNold3(:) = 0.
  ap(pp)%qNold4(:) = 0.
  
  ap(pp)%dxdtl(:) = 0.
  ap(pp)%dydtl(:) = 0.
  ap(pp)%dzdtl(:) = 0.
  ap(pp)%omg1(:) = 0.
  ap(pp)%omg2(:) = 0.
  ap(pp)%omg3(:) = 0.
   
  ap(pp)%ua(:)    = 0.
  ap(pp)%va(:)    = 0.
  ap(pp)%wa(:)    = 0.
  ap(pp)%omgd1(:) = 0.
  ap(pp)%omgd2(:) = 0.
  ap(pp)%omgd3(:) = 0.
  
  ap(pp)%fxll(:)  = 0.
  ap(pp)%fyll(:)  = 0.
  ap(pp)%fzll(:)  = 0.
  ap(pp)%mxl(:)  = 0.
  ap(pp)%myl(:)  = 0.
  ap(pp)%mzl(:)  = 0.

  ap(pp)%ull(:)    = 0.
  ap(pp)%vll(:)    = 0.
  ap(pp)%wll(:)    = 0.
  ap(pp)%omgl1(:) = 0.
  ap(pp)%omgl2(:) = 0.
  ap(pp)%omgl3(:) = 0.

  ap(pp)%dxdto(:) = 0.
  ap(pp)%dydto(:) = 0.
  ap(pp)%dzdto(:) = 0.
  ap(pp)%omgo1(:) = 0.
  ap(pp)%omgo2(:) = 0.
  ap(pp)%omgo3(:) = 0.
    
  ap(pp)%uao(:)    = 0.
  ap(pp)%vao(:)    = 0.
  ap(pp)%wao(:)    = 0.
  ap(pp)%omgdo1(:) = 0.
  ap(pp)%omgdo2(:) = 0.
  ap(pp)%omgdo3(:) = 0.

  ap(pp)%nb(:)     = 0
end forall
!
leftbound  = (coords(1)  )*lx/(1.*dims(1)) ! left  boundary of process myid
rightbound = (coords(1)+1)*lx/(1.*dims(1)) ! right boundary of process myid
frontbound = (coords(2)  )*ly/(1.*dims(2)) ! front boundary of process myid
backbound  = (coords(2)+1)*ly/(1.*dims(2)) ! back  boundary of process myid




!v what is the purpose of binsearch? to find master proc?
i = 0
count_mstr = 0
count_slve = 0
do idp = 1,np
  found_mstr = .false.

  call binsearch(idp,mslv_nb(1:pmax_nb(0),0),pmax_nb(0),found_mstr,p)

  !  do p=1,pmax
  !    if(sp(p)%mslv.eq.idp) then
  !      found_mstr = .true.
  if(found_mstr) then        !v below, why information is transferred from slave to master (sp to ap)?
    if (p==0) print*,idp,found_mstr,pmax_nb(0)
    i = i + 1
    ap(i)%mslv      = sp(p)%mslv

    ap(i)%x         = sp(p)%x
    ap(i)%y         = sp(p)%y
    ap(i)%z         = sp(p)%z

    ap(i)%integralx = sp(p)%integralx
    ap(i)%integraly = sp(p)%integraly
    ap(i)%integralz = sp(p)%integralz

    do l=1,nno
      ap(i)%xll(l)  = sp(p)%xll(l)
      ap(i)%yll(l)  = sp(p)%yll(l)
      ap(i)%zll(l)  = sp(p)%zll(l)
      ap(i)%qN1(l)  = sp(p)%qN1(l)
      ap(i)%qN2(l)  = sp(p)%qN2(l)
      ap(i)%qN3(l)  = sp(p)%qN3(l)
      ap(i)%qN4(l)  = sp(p)%qN4(l)
      
      ap(i)%xfpo(l) = sp(p)%xfpo(l)
      ap(i)%yfpo(l) = sp(p)%yfpo(l)
      ap(i)%zfpo(l) = sp(p)%zfpo(l)
      ap(i)%qNo1(l) = sp(p)%qNo1(l)
      ap(i)%qNo2(l) = sp(p)%qNo2(l)
      ap(i)%qNo3(l) = sp(p)%qNo3(l)
      ap(i)%qNo4(l) = sp(p)%qNo4(l)
      
      ap(i)%xfpold(l) = sp(p)%xfpold(l)
      ap(i)%yfpold(l) = sp(p)%yfpold(l)
      ap(i)%zfpold(l) = sp(p)%zfpold(l)
      ap(i)%qNold1(l) = sp(p)%qNold1(l)
      ap(i)%qNold2(l) = sp(p)%qNold2(l)
      ap(i)%qNold3(l) = sp(p)%qNold3(l)
      ap(i)%qNold4(l) = sp(p)%qNold4(l)
      
      ap(i)%dxdtl(l) = sp(p)%dxdtl(l)
      ap(i)%dydtl(l) = sp(p)%dydtl(l)
      ap(i)%dzdtl(l) = sp(p)%dzdtl(l)
      ap(i)%omg1(l) = sp(p)%omg1(l)
      ap(i)%omg2(l) = sp(p)%omg2(l)
      ap(i)%omg3(l) = sp(p)%omg3(l)
            
      ap(i)%ua(l)    = sp(p)%ua(l)
      ap(i)%va(l)    = sp(p)%va(l)
      ap(i)%wa(l)    = sp(p)%wa(l)
      ap(i)%omgd1(l) = sp(p)%omgd1(l)
      ap(i)%omgd2(l) = sp(p)%omgd2(l)
      ap(i)%omgd3(l) = sp(p)%omgd3(l)

      ap(i)%fxll(l)  = sp(p)%fxll(l)
      ap(i)%fyll(l)  = sp(p)%fyll(l)
      ap(i)%fzll(l)  = sp(p)%fzll(l)
      ap(i)%mxl(l)  = sp(p)%mxl(l)
      ap(i)%myl(l)  = sp(p)%myl(l)
      ap(i)%mzl(l)  = sp(p)%mzl(l)

      ap(i)%ull(l)   = sp(p)%ull(l)
      ap(i)%vll(l)   = sp(p)%vll(l)
      ap(i)%wll(l)   = sp(p)%wll(l)
      ap(i)%omgl1(l) = sp(p)%omgl1(l)
      ap(i)%omgl2(l) = sp(p)%omgl2(l)
      ap(i)%omgl3(l) = sp(p)%omgl3(l)

      ap(i)%dxdto(l) = sp(p)%dxdto(l)
      ap(i)%dydto(l) = sp(p)%dydto(l)
      ap(i)%dzdto(l) = sp(p)%dzdto(l)
      ap(i)%omgo1(l) = sp(p)%omgo1(l)
      ap(i)%omgo2(l) = sp(p)%omgo2(l)
      ap(i)%omgo3(l) = sp(p)%omgo3(l)
            
      ap(i)%uao(l)    = sp(p)%uao(l)
      ap(i)%vao(l)    = sp(p)%vao(l)
      ap(i)%wao(l)    = sp(p)%wao(l)
      ap(i)%omgdo1(l) = sp(p)%omgdo1(l)
      ap(i)%omgdo2(l) = sp(p)%omgdo2(l)
      ap(i)%omgdo3(l) = sp(p)%omgdo3(l)
    enddo

    count_mstr = count_mstr + 1

		do l=1,nno

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
          if ( (ap(i)%xll(l)+offset) .ge. rightbound ) then
            ap(i)%nb(8) = 1 !neighbor 8 is slave of particle ap(i)%mslv
          endif
        endif

      enddo

    else
      count_slve_loc = 0

      nb=5
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then

          do l=1,nno

            if ( (anb2(k,nb)%xll(l)+offset) .gt. leftbound ) then
              if ( ((anb2(k,nb)%yll(l)+offset) .ge. frontbound) .and. ((anb2(k,nb)%yll(l)-offset) .le. backbound ) ) then

                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1


                goto 206

              endif
            endif

          enddo
        endif

206 continue
      endif


      nb=6
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then

          do l=1,nno

            if ( (anb2(k,nb)%yll(l)-offset) .lt. backbound ) then
              if  ( (anb2(k,nb)%xll(l)+offset) .ge. leftbound ) then

                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1


                goto 207

              endif
            endif

          enddo
        endif

207 continue
      endif

  nb=7
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then

          do l=1,nno

            if ( (anb2(k,nb)%yll(l)-offset) .lt. backbound ) then
              if ( ((anb2(k,nb)%xll(l)+offset) .ge. leftbound) .and. ((anb2(k,nb)%xll(l)-offset) .le. rightbound) ) then

                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1


                goto 208

              endif
              endif

          enddo
        endif

208 continue
      endif


      nb=8
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then

          do l=1,nno

            if ( (anb2(k,nb)%yll(l)-offset) .lt. backbound ) then
              if  ( (anb2(k,nb)%xll(l)-offset) .le. rightbound ) then

                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1


                goto 201

              endif
            endif

          enddo
        endif

201 continue
      endif

   nb=1
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then

          do l=1,nno

            if ( (anb2(k,nb)%xll(l)-offset) .lt. rightbound ) then
              if ( ((anb2(k,nb)%yll(l)+offset) .ge. frontbound) .and. ((anb2(k,nb)%yll(l)-offset) .le. backbound) ) then

                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1


                goto 202

              endif
            endif

          enddo
        endif

202 continue
      endif


      nb=2
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then

          do l=1,nno

            if ( (anb2(k,nb)%xll(l)-offset) .lt. rightbound ) then
              if  ( (anb2(k,nb)%yll(l)+offset) .ge. frontbound ) then

                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1


                goto 203

              endif
            endif

          enddo
        endif

203 continue
      endif

  nb=3
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then

          do l=1,nno

            if ( (anb2(k,nb)%yll(l)+offset) .gt. frontbound ) then
              if ( ((anb2(k,nb)%xll(l)+offset) .ge. leftbound) .and. ((anb2(k,nb)%xll(l)-offset) .le. rightbound) ) then

                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1



                goto 204

              endif
            endif

          enddo
        endif

204 continue
      endif


      nb=4
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then

          do l=1,nno

            if ( (anb2(k,nb)%yll(l)+offset) .gt. frontbound ) then
              if  ( (anb2(k,nb)%xll(l)+offset) .ge. leftbound ) then

                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1



                goto 205

              endif
            endif

          enddo
        endif

205 continue
      endif

      if(count_slve_loc.ne.0) count_slve = count_slve + 1

    endif
enddo
!
! the new value of pmax is:
!
pmax   = count_mstr + count_slve
npmstr = count_mstr
!
!check if number of masters yield np
!
call MPI_ALLREDUCE(count_mstr, count_mstr_all, 1, MPI_INTEGER, MPI_SUM, comm_cart, ierr)


do pp=1,pmax
  if (ap(pp)%mslv .le. 0) then
    do l=1,nno
      ap(pp)%xll(l)    = 0.
      ap(pp)%yll(l)    = 0.
      ap(pp)%zll(l)    = 0.
      ap(pp)%qN1(:)    = 0.
      ap(pp)%qN2(:)    = 0.
      ap(pp)%qN3(:)    = 0.
      ap(pp)%qN4(:)    = 0.
      
      ap(pp)%xfpo(l)   = 0.
      ap(pp)%yfpo(l)   = 0.
      ap(pp)%zfpo(l)   = 0.
      ap(pp)%qNo1(:)   = 0.
      ap(pp)%qNo2(:)   = 0.
      ap(pp)%qNo3(:)   = 0.
      ap(pp)%qNo4(:)   = 0.
      
      ap(pp)%xfpold(l) = 0.
      ap(pp)%yfpold(l) = 0.
      ap(pp)%zfpold(l) = 0.
      ap(pp)%qNold1(:) = 0.
      ap(pp)%qNold2(:) = 0.
      ap(pp)%qNold3(:) = 0.
      ap(pp)%qNold4(:) = 0.
      
      ap(pp)%dxdtl(l)   = 0.
      ap(pp)%dydtl(l)   = 0.
      ap(pp)%dzdtl(l)   = 0.
      ap(pp)%omg1(:)   = 0.
      ap(pp)%omg2(:)   = 0.
      ap(pp)%omg3(:)   = 0.
      
      ap(pp)%ua(l)     = 0.
      ap(pp)%va(l)     = 0.
      ap(pp)%wa(l)     = 0.
      ap(pp)%omgd1(:)  = 0.
      ap(pp)%omgd2(:)  = 0.
      ap(pp)%omgd3(:)  = 0.

      ap(pp)%fxll(l)    = 0.
      ap(pp)%fyll(l)    = 0.
      ap(pp)%fzll(l)    = 0.
      ap(pp)%mxl(:)    = 0.
      ap(pp)%myl(:)    = 0.
      ap(pp)%mzl(:)    = 0.

      ap(pp)%ull(l)     = 0.
      ap(pp)%vll(l)     = 0.
      ap(pp)%wll(l)     = 0.      
      ap(pp)%omgl1(:)  = 0.
      ap(pp)%omgl2(:)  = 0.
      ap(pp)%omgl3(:)  = 0.
            
      ap(pp)%dxdto(:) = 0.
      ap(pp)%dydto(:) = 0.
      ap(pp)%dzdto(:) = 0.
      ap(pp)%omgo1(:) = 0.
      ap(pp)%omgo2(:) = 0.
      ap(pp)%omgo3(:) = 0.
        
      ap(pp)%uao(:)    = 0.
      ap(pp)%vao(:)    = 0.
      ap(pp)%wao(:)    = 0.
      ap(pp)%omgdo1(:) = 0.
      ap(pp)%omgdo2(:) = 0.
      ap(pp)%omgdo3(:) = 0.
    enddo
  endif
enddo

!write(6,'(A7,I5,A8,I5,A18,I5,A11,A8,I5)') 'Thread ', myid, ' masters ', count_mstr, &
!' and is slave for ', count_slve, ' particles. ', ' pmax = ', pmax

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
!  write(6,*) 'The particles were successfully initialized in thread ', myid, ' !'
endif


return
end subroutine Update_Pos



subroutine binsearch(ival,array,idim,found,index)
integer, intent(in), dimension(1:) :: array
integer, intent(in) :: ival,idim
logical, intent(out) :: found
integer, intent(out) :: index
integer :: start,finish,range,mid

start = 1
finish = idim
range = finish-start
mid = (start+finish)/2

do while( abs(array(mid)) .ne. ival .and. range .gt.  0)

  if (ival .gt. abs(array(mid))) then
    start = mid + 1
  else
    finish = mid - 1
  endif
  range = finish - start
  mid = (start + finish)/2
enddo

if(array(mid).ne.ival) then
  found = .false.
  index = 0
else
  found = .true.
  index = mid
endif

return
end subroutine binsearch


end module mod_Update_Pos
