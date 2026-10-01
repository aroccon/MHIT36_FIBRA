module mod_suspension
use mod_common_fibm
use mod_common_mpi
use mod_interp_spread
implicit none
private
public suspension
contains
!
subroutine suspension
!**************************************************************************************************************************
!the subroutine saves normal stress differences and mean value of end to end and maximum distance for filaments
!**************************************************************************************************************************
!
implicit none
integer i,j,k,p,xcoor,zcoor
real::st1,st2,st3,st4,stt,sigma11,sigma22,sigma33,sigma23
real::st5,st6,st7,st8,st9,st10,st11,st12,st13,st14,st15,st16
real::st17,st18,st19,st20,stf1,stf2,stf3,stf4,factor
real::sm,sm1,sm2,sm3,sm4,eeg,eel,eee,org,orll,dd
real::maxf,maxl,maxg,el,eg,dis,lwbound1,lwbound2,minl,ming,dzf
real:: dxf,ddd,ncpav,ncpavl
real,dimension(nl)::dxds,dyds,dzds
real,dimension(kmax)::vmean,vmeanl,uu,vv,ww,vw
real,dimension(kmax)::uul,vvl,wwl,vwl
real,dimension(itot/nfd,ktot/nfd)::fd,fdl,ed,edl,or,orl
real,dimension(np)::ncpl,ncp
real ,dimension(0:i1,0:j1,0:k1) :: vd
character(7) filenumber
real :: leftbound,rightbound,frontbound,backbound
!

leftbound   = (coords(1)  )*lx/(1.*dims(1)) !left  boundary of process myid
rightbound  = (coords(1)+1)*lx/(1.*dims(1)) !right boundary of process myid
frontbound  = (coords(2)  )*ly/(1.*dims(2)) !front boundary of process myid
backbound   = (coords(2)+1)*ly/(1.*dims(2)) !back  boundary of process myid



stt=0.
el=0.;eg=0.;st1=0.;st2=0.;st3=0.;st4=0.;stf1=0.;stf2=0.;stf3=0.;stf4=0.
st5=0.;st6=0.;st7=0.;st8=0.;st9=0.;st10=0.;st11=0.;st12=0.
!******************************************************************computing normal stress differences**********************************************************************
!**************************************************computing viscous forces whole domain ********************************************
!do k=1,kmax
!do j=1,jmax
!do i=1,imax
!s11(i,j,k)=-pnew(i,j,k)+(2.*visc0*(unew(i,j,k)-unew(i-1,j,k))*dxi)
!s22(i,j,k)=-pnew(i,j,k)+(2.*visc0*(vnew(i,j,k)-vnew(i,j-1,k))*dyi)
!s33(i,j,k)=-pnew(i,j,k)+(2.*visc0*(wnew(i,j,k)-wnew(i,j,k-1))*dzi)
!s23(i,j,k)=(visc0*((((vnew(i,j,k+1)+vnew(i,j-1,k+1))-(vnew(i,j,k-1)+vnew(i,j-1,k-1)))*dzi*.25)+(((wnew(i,j+1,k)+wnew(i,j+1,k-1))- & 
!(wnew(i,j-1,k)+wnew(i,j-1,k-1)))*dyi*.25)))
!st1=st1+s11(i,j,k)*dx*dy*dz
!st2=st2+s22(i,j,k)*dx*dy*dz
!st3=st3+s33(i,j,k)*dx*dy*dz
!st4=st4+s23(i,j,k)*dx*dy*dz
!enddo
!enddo
!enddo
!call mpi_allreduce(st1,stt,1,mpi_real8,mpi_sum,comm_cart,error)
!st1=stt
!call mpi_allreduce(st2,stt,1,mpi_real8,mpi_sum,comm_cart,error)
!st2=stt
!call mpi_allreduce(st3,stt,1,mpi_real8,mpi_sum,comm_cart,error)
!st3=stt
!call mpi_allreduce(st4,stt,1,mpi_real8,mpi_sum,comm_cart,error)
!st4=stt
!**************************************************computing viscous forces inside particles ********************************************
!call eulr2lagrst

!do p=1,pmax
!if (ap(p)%mslv>0) then
!do l=1,nl
!if (l==1.or.l==nl) then
!st5=st5+(ap(p)%s11l(l)*pi*.25*(aspectratio**2.)*ds*.5)
!st6=st6+(ap(p)%s22l(l)*pi*.25*(aspectratio**2.)*ds*.5)
!st7=st7+(ap(p)%s33l(l)*pi*.25*(aspectratio**2.)*ds*.5)
!st8=st8+(ap(p)%s23l(l)*pi*.25*(aspectratio**2.)*ds*.5)
!else
!st5=st5+(ap(p)%s11l(l)*pi*.25*(aspectratio**2.)*ds)
!st6=st6+(ap(p)%s22l(l)*pi*.25*(aspectratio**2.)*ds)
!st7=st7+(ap(p)%s33l(l)*pi*.25*(aspectratio**2.)*ds)
!st8=st8+(ap(p)%s23l(l)*pi*.25*(aspectratio**2.)*ds)
!endif
!enddo
!endif
!enddo
!call mpi_allreduce(st5,stt,1,mpi_real8,mpi_sum,comm_cart,error)
!st5=stt
!call mpi_allreduce(st6,stt,1,mpi_real8,mpi_sum,comm_cart,error)
!st6=stt
!call mpi_allreduce(st7,stt,1,mpi_real8,mpi_sum,comm_cart,error)
!st7=stt
!call mpi_allreduce(st8,stt,1,mpi_real8,mpi_sum,comm_cart,error)
!st8=stt
!*********************************************************total fluid stress************************************************************
!stf1=st1/(lx*ly*lz)
!stf2=st2/(lx*ly*lz)
!stf3=st3/(lx*ly*lz)
!stf4=st4/(lx*ly*lz)

!**************************************************computing contribution of FSI forces **********************************************
!!v yet to modify this part (to include moment part)
do p=1,pmax
if (ap(p)%mslv>0) then
do l=1,nl
if (l==1.or.l==nl) then
st1=st1-(ap(p)%fxl(l)*ap(p)%xfp(l)*ds*pi*.125*(aspectratio**2.))
st2=st2-(ap(p)%fyl(l)*ap(p)%yfp(l)*ds*pi*.125*(aspectratio**2.))
st3=st3-(ap(p)%fzl(l)*ap(p)%zfp(l)*ds*pi*.125*(aspectratio**2.))
st4=st4-(ap(p)%fyl(l)*ap(p)%zfp(l)*ds*pi*.125*(aspectratio**2.))
else
st1=st1-(ap(p)%fxl(l)*ap(p)%xfp(l)*ds*pi*.25*(aspectratio**2.))
st2=st2-(ap(p)%fyl(l)*ap(p)%yfp(l)*ds*pi*.25*(aspectratio**2.))
st3=st3-(ap(p)%fzl(l)*ap(p)%zfp(l)*ds*pi*.25*(aspectratio**2.))
st4=st4-(ap(p)%fyl(l)*ap(p)%zfp(l)*ds*pi*.25*(aspectratio**2.))
endif

enddo
endif
enddo
call mpi_allreduce(st1,stt,1,mpi_real8,mpi_sum,comm_cart,ierr)
st1=stt/(lx*ly*lz)
call mpi_allreduce(st2,stt,1,mpi_real8,mpi_sum,comm_cart,ierr)
st2=stt/(lx*ly*lz)
call mpi_allreduce(st3,stt,1,mpi_real8,mpi_sum,comm_cart,ierr)
st3=stt/(lx*ly*lz)
call mpi_allreduce(st4,stt,1,mpi_real8,mpi_sum,comm_cart,ierr)
st4=stt/(lx*ly*lz)
!**************************************************computing contribution of lubrication forces **********************************************
do p=1,pmax
if (ap(p)%mslv>0) then
do l=1,nl
if (l==1.or.l==nl) then
st5=st5+(ap(p)%fcx(l)*ap(p)%xfp(l)*ds*pi*.125*(aspectratio**2.))
st6=st6+(ap(p)%fcy(l)*ap(p)%yfp(l)*ds*pi*.125*(aspectratio**2.))
st7=st7+(ap(p)%fcz(l)*ap(p)%zfp(l)*ds*pi*.125*(aspectratio**2.))
st8=st8+(ap(p)%fcy(l)*ap(p)%zfp(l)*ds*pi*.125*(aspectratio**2.))
else
st5=st5+(ap(p)%fcx(l)*ap(p)%xfp(l)*ds*pi*.25*(aspectratio**2.))
st6=st6+(ap(p)%fcy(l)*ap(p)%yfp(l)*ds*pi*.25*(aspectratio**2.))
st7=st7+(ap(p)%fcz(l)*ap(p)%zfp(l)*ds*pi*.25*(aspectratio**2.))
st8=st8+(ap(p)%fcy(l)*ap(p)%zfp(l)*ds*pi*.25*(aspectratio**2.))
endif

enddo
endif
enddo
call mpi_allreduce(st5,stt,1,mpi_real8,mpi_sum,comm_cart,ierr)
st5=stt/(lx*ly*lz)
call mpi_allreduce(st6,stt,1,mpi_real8,mpi_sum,comm_cart,ierr)
st6=stt/(lx*ly*lz)
call mpi_allreduce(st7,stt,1,mpi_real8,mpi_sum,comm_cart,ierr)
st7=stt/(lx*ly*lz)
call mpi_allreduce(st8,stt,1,mpi_real8,mpi_sum,comm_cart,ierr)
st8=stt/(lx*ly*lz)
!********************************************computing contribution of filament acceleration *****************************************
!st17=0.;st18=0.;st19=0.;st20=0.
!factor=.25*pi*(aspectratio**2.)
!do p=1,pmax
!if (ap(p)%mslv>0) then
!do l=1,nl
!if (l==1.or.l==nl) then
!st17=st17+(ap(p)%sxx(l)*ap(p)%xfp(l)*ds*factor*.5)
!st18=st18+(ap(p)%syy(l)*ap(p)%yfp(l)*ds*factor*.5)
!st19=st19+(ap(p)%szz(l)*ap(p)%zfp(l)*ds*factor*.5)
!st20=st20+(ap(p)%syy(l)*ap(p)%zfp(l)*ds*factor*.5)
!else
!st17=st17+(ap(p)%sxx(l)*ap(p)%xfp(l)*ds*factor)
!st18=st18+(ap(p)%syy(l)*ap(p)%yfp(l)*ds*factor)
!st19=st19+(ap(p)%szz(l)*ap(p)%zfp(l)*ds*factor)
!st20=st20+(ap(p)%syy(l)*ap(p)%zfp(l)*ds*factor)
!endif

!enddo
!endif
!enddo
!call mpi_allreduce(st17,stt,1,mpi_real8,mpi_sum,comm_cart,error)
!st17=stt/(lx*ly*lz)
!call mpi_allreduce(st18,stt,1,mpi_real8,mpi_sum,comm_cart,error)
!st18=stt/(lx*ly*lz)
!call mpi_allreduce(st19,stt,1,mpi_real8,mpi_sum,comm_cart,error)
!st19=stt/(lx*ly*lz)
!call mpi_allreduce(st20,stt,1,mpi_real8,mpi_sum,comm_cart,error)
!st20=stt/(lx*ly*lz)

!**************************************************writing the stresses**********************************************
sigma11=st1!-st17
sigma22=st2!-st18
sigma33=st3!-st19
sigma23=st4!-st20


if (myid==0) then


 open(71,file=datadir//'stress.txt',position='append')
 write(71,'(6E16.8)') time,sigma11,sigma22,sigma33,sigma23,&
 ((1./3.)*(sigma11+sigma22+sigma33))
 close(71)


! open(73,file=datadir//'normalstress_dif.txt',position='append')
! write(73,'(7E16.8)') time,1.+(sigma23/visc0),sigma22-sigma33,sigma33-sigma11, &
! 1.+((sigma23+st8)/visc0),sigma22-sigma33+st6-st7,sigma33-sigma11+st7-st5
! close(73)
endif

!*******************************************************computing mean value of maximum and end to end distance**********************************************************************
maxl=0.;maxg=0.;eel=0.;eeg=0.
minl=1.;ming=1.

!***************************************************maximum distance******************************************
do p=1,pmax
if (ap(p)%mslv>0) then
maxf=0.
do i=1,nl
do j=i,nl
if (i/=j) then
maxf=max(maxf,sqrt(((ap(p)%xfp(i)-ap(p)%xfp(j))**2.)+((ap(p)%yfp(i)-ap(p)%yfp(j))**2.)+&
((ap(p)%zfp(i)-ap(p)%zfp(j))**2.)))
endif
enddo
enddo
maxl=maxl+maxf

sm=0.
do l=2,NL-1
sm=sm+((((ap(p)%xfp(l+1)+ap(p)%xfp(l-1)-(2.*ap(p)%xfp(l)))/(ds**2.))**2.)+ &
(((ap(p)%yfp(l+1)+ap(p)%yfp(l-1)-(2.*ap(p)%yfp(l)))/(ds**2.))**2.)+&
(((ap(p)%zfp(l+1)+ap(p)%zfp(l-1)-(2.*ap(p)%zfp(l)))/(ds**2.))**2.))*.5*ds
enddo

el=el+sm

endif
enddo


call mpi_allreduce(maxl,maxg,1,mpi_real8,mpi_sum,comm_cart,ierr)
call mpi_allreduce(el,eg,1,mpi_real8,mpi_sum,comm_cart,ierr)

maxg=maxg/np
!**************************************************filament distribution****************************************
fd=0.
fdl=0.
dzf=(lz*nfd)/ktot
dxf=(lx*nfd)/itot
do p=1,pmax
if (ap(p)%mslv>0) then



xcoor=int(ap(p)%x/dxf)+1
zcoor=int(ap(p)%z/dzf)+1 

if (xcoor<1) xcoor=1
if (zcoor<1) zcoor=1
if (xcoor>(itot/nfd)) xcoor=(itot/nfd)
if (zcoor>(ktot/nfd)) zcoor=(ktot/nfd)


 fdl(xcoor,zcoor)=fdl(xcoor,zcoor)+1.


endif
enddo

do i=1,itot/nfd
do k=1,ktot/nfd
call mpi_allreduce(fdl(i,k),fd(i,k),1,mpi_real8,mpi_sum,comm_cart,ierr)
enddo
enddo

!**************************************************end to end distance****************************************
ed=0.
edl=0.
do p=1,pmax
if (ap(p)%mslv>0) then


eee=sqrt(((ap(p)%xfp(nl)-ap(p)%xfp(1))**2.)+((ap(p)%yfp(nl)-ap(p)%yfp(1))**2.)+&
((ap(p)%zfp(nl)-ap(p)%zfp(1))**2.))
eel=eel+eee
minl=min(minl,eee)

xcoor=int(ap(p)%x/dxf)+1
zcoor=int(ap(p)%z/dzf)+1


if (xcoor<1) xcoor=1
if (zcoor<1) zcoor=1
if (xcoor>(itot/nfd)) xcoor=(itot/nfd)
if (zcoor>(ktot/nfd)) zcoor=(ktot/nfd)


edl(xcoor,zcoor)=edl(xcoor,zcoor)+eee


endif
enddo
call mpi_allreduce(eel,eeg,1,mpi_real8,mpi_sum,comm_cart,ierr)
call mpi_allreduce(minl,ming,1,mpi_real8,mpi_min,comm_cart,ierr)
eeg=eeg/np
do i=1,itot/nfd
do k=1,ktot/nfd
call mpi_allreduce(edl(i,k),ed(i,k),1,mpi_real8,mpi_sum,comm_cart,ierr)
enddo
enddo

do i=1,itot/nfd
do k=1,ktot/nfd
if (fd(i,k)/=0) then
ed(i,k)=ed(i,k)/fd(i,k)
else
ed(i,k)=0.
endif
enddo
enddo
!****************************************************orientation*********************************************
or=0.
orl=0.
org=0.
orll=0.
do p=1,pmax
if (ap(p)%mslv>0) then

xcoor=int(ap(p)%x/dxf)+1
zcoor=int(ap(p)%z/dzf)+1


if (xcoor<1) xcoor=1
if (zcoor<1) zcoor=1
if (xcoor>(itot/nfd)) xcoor=(itot/nfd)
if (zcoor>(ktot/nfd)) zcoor=(ktot/nfd)

dd=sqrt(((ap(p)%xfp(nl)-ap(p)%xfp(1))**2.)+((ap(p)%yfp(nl)-ap(p)%yfp(1))**2.)+&
((ap(p)%zfp(nl)-ap(p)%zfp(1))**2.))

ddd=sqrt(((ap(p)%xfp(nl)-ap(p)%xfp(1))**2.)+&
((ap(p)%zfp(nl)-ap(p)%zfp(1))**2.))

eee=asin(min(ddd,1.)/dd)/(pi*.5)

orl(xcoor,zcoor)=orl(xcoor,zcoor)+eee

orll=orll+eee

endif
enddo
call mpi_allreduce(orll,org,1,mpi_real8,mpi_sum,comm_cart,ierr)
org=org/np
do i=1,itot/nfd
do k=1,ktot/nfd
call mpi_allreduce(orl(i,k),or(i,k),1,mpi_real8,mpi_sum,comm_cart,ierr)
enddo
enddo

do i=1,itot/nfd
do k=1,ktot/nfd
if (fd(i,k)/=0) then
or(i,k)=or(i,k)/fd(i,k)
else
or(i,k)=0.
endif
enddo
enddo
!****************************************************number of contact points*********************************************
ncp=0.;ncpl=0.;ncpavl=0.;ncpav=0.
do p=1,pmax
if (ap(p)%mslv>0) then
do l=1,nl
if (ap(p)%cp(l)>0.) then
ncpavl=ncpavl+ap(p)%cp(l)
ncpl(ap(p)%mslv)=ncpl(ap(p)%mslv)+ap(p)%cp(l)
endif
enddo
endif
enddo
call mpi_allreduce(ncpavl,ncpav,1,mpi_real8,mpi_sum,comm_cart,ierr)
ncpav=ncpav/np


do p=1,np
call mpi_allreduce(ncpl(p),ncp(p),1,mpi_real8,mpi_sum,comm_cart,ierr)
enddo
!*************************************************************************************************************************

if (myid==0) then


!!open(73,file=datadir//'distance.txt',position='append')
!!write(73,'(4E16.8)') time,eeg,maxg,ming
!!close(73)
!!
!!
!!open(74,file=datadir//'energy.txt',position='append')
!!write(74,'(2E16.8)') time,eg/np
!! close(74)
!!
!!
!!
!!open(75,file=datadir//'orientation.txt',position='append')
!!write(75,'(2E16.8)') time,org
!!close(75)
!!
!!open(76,file=datadir//'contact_points.txt',position='append')
!!write(76,'(2E16.8)') time,ncpav
!!close(76)

!!if (mod(istep,2000)==0.or.istep==1) then
!!
!!write(filenumber,'(i7.7)') istep
!!open(77,file=datadir//'fdis'//filenumber//'.txt')
!!do i=1,itot/nfd
!!do k=1,ktot/nfd
!!write(77,'(5E16.8)') (i-.5)*dxf,(k-.5)*dzf,(fd(i,k)*(itot/nfd)*(ktot/nfd))/np,ed(i,k),or(i,k)
!!enddo
!!enddo
!! close(77)
!!
!!endif




!!if (mod(istep,10)==0) then
!!
!!write(filenumber,'(i7.7)') istep
!!open(78,file=datadir//'cp'//filenumber//'.txt')
!!do p=1,np
!!write(78,'(E16.8)') ncp(p)
!!enddo
!!close(78)
!!endif





endif

return
end subroutine suspension
! 
end module mod_suspension
