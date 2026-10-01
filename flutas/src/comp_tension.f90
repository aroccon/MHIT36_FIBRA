module mod_comp_tension
use mpi
use mod_common_fibm
use mod_common_mpi
use mod_param_fibm
implicit none
private
public comp_tension
contains
!
subroutine comp_tension(p)
implicit none
integer,intent(in)::p
integer::l
real,dimension(NL-1)::  d2xdsdt,d2ydsdt,d2zdsdt,dxdsold,dydsold,dzdsold
real,dimension(NL-1)::  d2dt2xs,dxdsnew,dydsnew,dzdsnew
real,dimension(NL-1)::  dfxds,dfyds,dfzds,bt,st,duads,dvads,dwads,dacc
real,dimension(NL-1)::  dfcxds,dfcyds,dfczds
real,dimension(NL-1)::  dxds,dyds,dzds,dxdsp,dydsp,dzdsp
real,dimension(NL-1)::  dsc
real,dimension(NL-1,NL-1)::ft









!**************************************************p_data is predicted data****************************************************
if (inertia==0) then
  !    p_data(p)%xfp(:)=ap(p)%xfp(:)+(ap(p)%ul(:)*dt)
  !    p_data(p)%yfp(:)=ap(p)%yfp(:)+(ap(p)%vl(:)*dt)
  !    p_data(p)%zfp(:)=ap(p)%zfp(:)+(ap(p)%wl(:)*dt)
  p_data(p)%xfp(:)=(2.*ap(p)%xfpo(:))-ap(p)%xfpold(:)
  p_data(p)%yfp(:)=(2.*ap(p)%yfpo(:))-ap(p)%yfpold(:)
  p_data(p)%zfp(:)=(2.*ap(p)%zfpo(:))-ap(p)%zfpold(:)
elseif (inertia==1) then 
  p_data(p)%xfp(:)=(2.*ap(p)%xfpo(:))-ap(p)%xfpold(:)
  p_data(p)%yfp(:)=(2.*ap(p)%yfpo(:))-ap(p)%yfpold(:)
  p_data(p)%zfp(:)=(2.*ap(p)%zfpo(:))-ap(p)%zfpold(:)
endif

 
do l=1,NL
  if (l==1) then
    !******************************************************free*************************************************************        
    if (fboundary=="f") then
      d3xds3(l)=(p_data(p)%xfp(l+2)+p_data(p)%xfp(l)-(2.0*p_data(p)%xfp(l+1)))/(ds**3.0)
      d3yds3(l)=(p_data(p)%yfp(l+2)+p_data(p)%yfp(l)-(2.0*p_data(p)%yfp(l+1)))/(ds**3.0)
      d3zds3(l)=(p_data(p)%zfp(l+2)+p_data(p)%zfp(l)-(2.0*p_data(p)%zfp(l+1)))/(ds**3.0)
    !*****************************************************clamped***********************************************************
    elseif (fboundary=="c") then
      d3xds3(l)=(((p_data(p)%xfp(l)+p_data(p)%xfp(l+2)-(2.0*p_data(p)%xfp(l+1)))/(ds**2.0)) &
      -((p_data(p)%xfp(l+2)-p_data(p)%xfp(l))/(2.*(ds**2.))))/ds
      d3yds3(l)=(((p_data(p)%yfp(l)+p_data(p)%yfp(l+2)-(2.0*p_data(p)%yfp(l+1)))/(ds**2.0)) &
      -((p_data(p)%yfp(l+2)-p_data(p)%yfp(l))/(2.*(ds**2.))))/ds
      d3zds3(l)=(((p_data(p)%zfp(l)+p_data(p)%zfp(l+2)-(2.0*p_data(p)%zfp(l+1)))/(ds**2.0)) &
      -((((p_data(p)%zfp(l+2)-p_data(p)%zfp(l))/(2.*ds))-1.)/ds))/ds
    !******************************************************hinged***********************************************************
    elseif (fboundary=="h") then
      d3xds3(l)=(p_data(p)%xfp(l)+p_data(p)%xfp(l+2)-(2.0*p_data(p)%xfp(l+1)))/(ds**3.0)
      d3yds3(l)=(p_data(p)%yfp(l)+p_data(p)%yfp(l+2)-(2.0*p_data(p)%yfp(l+1)))/(ds**3.0)
      d3zds3(l)=(p_data(p)%zfp(l)+p_data(p)%zfp(l+2)-(2.0*p_data(p)%zfp(l+1)))/(ds**3.0)
    endif
  elseif (l==2) then
    !******************************************************hinged***********************************************************
    if (fboundary=="h".or.fboundary=="f") then
      d3xds3(l)=(p_data(p)%xfp(l)+p_data(p)%xfp(l+2)-(2.0*p_data(p)%xfp(l+1)))/(2.*(ds**3.0))
      d3yds3(l)=(p_data(p)%yfp(l)+p_data(p)%yfp(l+2)-(2.0*p_data(p)%yfp(l+1)))/(2.*(ds**3.0))
      d3zds3(l)=(p_data(p)%zfp(l)+p_data(p)%zfp(l+2)-(2.0*p_data(p)%zfp(l+1)))/(2.*(ds**3.0))
    !******************************************************clamped**********************************************************
    elseif (fboundary=="c") then
      d3xds3(l)=(((p_data(p)%xfp(l)+p_data(p)%xfp(l+2)-(2.0*p_data(p)%xfp(l+1)))/(ds**2.0)) &
      -((p_data(p)%xfp(l+1)-p_data(p)%xfp(l-1))/(2.*(ds**2.))))/(2.*ds)
      d3yds3(l)=(((p_data(p)%yfp(l)+p_data(p)%yfp(l+2)-(2.0*p_data(p)%yfp(l+1)))/(ds**2.0)) &
      -((p_data(p)%yfp(l+1)-p_data(p)%yfp(l-1))/(2.*(ds**2.))))/(2.*ds)
      d3zds3(l)=(((p_data(p)%zfp(l)+p_data(p)%zfp(l+2)-(2.0*p_data(p)%zfp(l+1)))/(ds**2.0)) &
      -((((p_data(p)%zfp(l+1)-p_data(p)%zfp(l-1))/(2.*ds))-1.)/ds))/(2.*ds)
    endif
  elseif (l==NL-1) then
    d3xds3(l)=-(p_data(p)%xfp(l+1)+p_data(p)%xfp(l-1)-(2.0*p_data(p)%xfp(l)))/(ds**3.0)
    d3yds3(l)=-(p_data(p)%yfp(l+1)+p_data(p)%yfp(l-1)-(2.0*p_data(p)%yfp(l)))/(ds**3.0)
    d3zds3(l)=-(p_data(p)%zfp(l+1)+p_data(p)%zfp(l-1)-(2.0*p_data(p)%zfp(l)))/(ds**3.0)
  elseif (l==NL) then
    d3xds3(l)= 0.
    d3yds3(l)= 0.
    d3zds3(l)= 0.
  else
    d3xds3(l)=((p_data(p)%xfp(l)+p_data(p)%xfp(l+2)-(2.0*p_data(p)%xfp(l+1)))-(p_data(p)%xfp(l)+p_data(p)%xfp(l-2) &
              -(2.0*p_data(p)%xfp(l-1))))/(2.0*(ds**3.0))
    d3yds3(l)=((p_data(p)%yfp(l)+p_data(p)%yfp(l+2)-(2.0*p_data(p)%yfp(l+1)))-(p_data(p)%yfp(l)+p_data(p)%yfp(l-2) &
              -(2.0*p_data(p)%yfp(l-1))))/(2.0*(ds**3.0))
    d3zds3(l)=((p_data(p)%zfp(l)+p_data(p)%zfp(l+2)-(2.0*p_data(p)%zfp(l+1)))-(p_data(p)%zfp(l)+p_data(p)%zfp(l-2) &
              -(2.0*p_data(p)%zfp(l-1))))/(2.0*(ds**3.0))
  endif
enddo

do l=1,NL
  if (l==1) then
    d4xds4(l)=(d3xds3(l+1)-d3xds3(l))/ds
    d4yds4(l)=(d3yds3(l+1)-d3yds3(l))/ds
    d4zds4(l)=(d3zds3(l+1)-d3zds3(l))/ds
  elseif (l==NL) then
    d4xds4(l)= -d3xds3(l-1)/ds
    d4yds4(l)= -d3yds3(l-1)/ds
    d4zds4(l)= -d3zds3(l-1)/ds
  else
    d4xds4(l)=(d3xds3(l+1)-d3xds3(l-1))/(2.0*ds)
    d4yds4(l)=(d3yds3(l+1)-d3yds3(l-1))/(2.0*ds)
    d4zds4(l)=(d3zds3(l+1)-d3zds3(l-1))/(2.0*ds)
  endif
enddo


do l=1,NL
  if (l==1) then
    d5xds5(l)=(d4xds4(l+1)-d4xds4(l))/ds
    d5yds5(l)=(d4yds4(l+1)-d4yds4(l))/ds
    d5zds5(l)=(d4zds4(l+1)-d4zds4(l))/ds
  elseif (l==NL) then
    d5xds5(l)= (d4xds4(l)-d4xds4(l-1))/ds
    d5yds5(l)= (d4yds4(l)-d4yds4(l-1))/ds
    d5zds5(l)= (d4zds4(l)-d4zds4(l-1))/ds
  else
    d5xds5(l)=(d4xds4(l+1)-d4xds4(l-1))/(2.0*ds)
    d5yds5(l)=(d4yds4(l+1)-d4yds4(l-1))/(2.0*ds)
    d5zds5(l)=(d4zds4(l+1)-d4zds4(l-1))/(2.0*ds)
  endif
enddo


     

!    ap(p)%xcur(:)=d4xds4(:)
!    ap(p)%ycur(:)=d4yds4(:)
!    ap(p)%zcur(:)=d4zds4(:)


do l=1,NL-1
  dxds(l) = ( ap(p)%xfpo(l+1) - ap(p)%xfpo(l) ) / ds      
  dyds(l) = ( ap(p)%yfpo(l+1) - ap(p)%yfpo(l) ) / ds
  dzds(l) = ( ap(p)%zfpo(l+1) - ap(p)%zfpo(l) ) / ds
  dxdsold(l) = ( ap(p)%xfpold(l+1) - ap(p)%xfpold(l) ) / ds      
  dydsold(l) = ( ap(p)%yfpold(l+1) - ap(p)%yfpold(l) ) / ds
  dzdsold(l) = ( ap(p)%zfpold(l+1) - ap(p)%zfpold(l) ) / ds
  dxdsp(l) = ( p_data(p)%xfp(l+1) - p_data(p)%xfp(l) ) / ds      
  dydsp(l) = ( p_data(p)%yfp(l+1) - p_data(p)%yfp(l) ) / ds
  dzdsp(l) = ( p_data(p)%zfp(l+1) - p_data(p)%zfp(l) ) / ds
  d2xdsdt(l)= (ap(p)%dxdt(l+1)-ap(p)%dxdt(l))/ds
  d2ydsdt(l)= (ap(p)%dydt(l+1)-ap(p)%dydt(l))/ds
  d2zdsdt(l)= (ap(p)%dzdt(l+1)-ap(p)%dzdt(l))/ds
  duads(l)=   (ap(p)%ua(l+1)-ap(p)%ua(l))/ds
  dvads(l)=   (ap(p)%va(l+1)-ap(p)%va(l))/ds
  dwads(l)=   (ap(p)%wa(l+1)-ap(p)%wa(l))/ds
  dacc(l)=    (duads(l)*dxdsp(l))+(dvads(l)*dydsp(l))+(dwads(l)*dzdsp(l))
  d2dt2xs(l)=(1.-(2.*((dxds(l)**2.)+(dyds(l)**2.)+(dzds(l)**2.)))+((dxdsold(l)**2.)+ &
              (dydsold(l)**2.)+(dzdsold(l)**2.))) /(2.*(dt**2.))
  dfxds(l)= (ap(p)%fxl(l+1) - ap(p)%fxl(l))/ds
  dfyds(l)= (ap(p)%fyl(l+1) - ap(p)%fyl(l))/ds
  dfzds(l)= (ap(p)%fzl(l+1) - ap(p)%fzl(l))/ds
  dfcxds(l)= (ap(p)%fcx(l+1) - ap(p)%fcx(l))/ds
  dfcyds(l)= (ap(p)%fcy(l+1) - ap(p)%fcy(l))/ds
  dfczds(l)= (ap(p)%fcz(l+1) - ap(p)%fcz(l))/ds
enddo
  
  
st= 0.
aex= 0.
awx= 0.
apx= 0.
ft= 0.
bt= 0.

!if (inertia==1) then


do l=1,NL-1
  if (l==1) then
    if (fboundary=="h".or.fboundary=="c") then 
      apx(l)=-(((p_data(p)%xfp(l+1)-p_data(p)%xfp(l))/(ds**3.))*dxdsp(l)) &
              -(((p_data(p)%yfp(l+1)-p_data(p)%yfp(l))/(ds**3.))*dydsp(l)) &
              -(((p_data(p)%zfp(l+1)-p_data(p)%zfp(l))/(ds**3.))*dzdsp(l))
      aex(l)= (((p_data(p)%xfp(l+2)-p_data(p)%xfp(l+1))/(ds**3.))*dxdsp(l)) &
              +(((p_data(p)%yfp(l+2)-p_data(p)%yfp(l+1))/(ds**3.))*dydsp(l)) &
              +(((p_data(p)%zfp(l+2)-p_data(p)%zfp(l+1))/(ds**3.))*dzdsp(l))

    elseif (fboundary=="f") then

      apx(l)=3.*(-(((p_data(p)%xfp(l+1)-p_data(p)%xfp(l))/(ds**3.))*dxdsp(l)) &
                  -(((p_data(p)%yfp(l+1)-p_data(p)%yfp(l))/(ds**3.))*dydsp(l)) &
                  -(((p_data(p)%zfp(l+1)-p_data(p)%zfp(l))/(ds**3.))*dzdsp(l)))
      aex(l)= (((p_data(p)%xfp(l+2)-p_data(p)%xfp(l+1))/(ds**3.))*dxdsp(l)) &
              +(((p_data(p)%yfp(l+2)-p_data(p)%yfp(l+1))/(ds**3.))*dydsp(l)) &
              +(((p_data(p)%zfp(l+2)-p_data(p)%zfp(l+1))/(ds**3.))*dzdsp(l))
    endif
  elseif (l==NL-1) then

    apx(l)=3.*(-(((p_data(p)%xfp(l+1)-p_data(p)%xfp(l))/(ds**3.))*dxdsp(l)) &
                -(((p_data(p)%yfp(l+1)-p_data(p)%yfp(l))/(ds**3.))*dydsp(l)) &
                -(((p_data(p)%zfp(l+1)-p_data(p)%zfp(l))/(ds**3.))*dzdsp(l)))
    awx(l)= (((p_data(p)%xfp(l)-p_data(p)%xfp(l-1))/(ds**3.))*dxdsp(l)) &
            +(((p_data(p)%yfp(l)-p_data(p)%yfp(l-1))/(ds**3.))*dydsp(l)) &
            +(((p_data(p)%zfp(l)-p_data(p)%zfp(l-1))/(ds**3.))*dzdsp(l))
  else
    
    apx(l)=-2.*((((p_data(p)%xfp(l+1)-p_data(p)%xfp(l))*dxdsp(l))/(ds**3.)) &
                +(((p_data(p)%yfp(l+1)-p_data(p)%yfp(l))*dydsp(l))/(ds**3.)) &
                +(((p_data(p)%zfp(l+1)-p_data(p)%zfp(l))*dzdsp(l))/(ds**3.)))    
    aex(l)= (((p_data(p)%xfp(l+2)-p_data(p)%xfp(l+1))*dxdsp(l))/(ds**3.)) &
            +(((p_data(p)%yfp(l+2)-p_data(p)%yfp(l+1))*dydsp(l))/(ds**3.)) &
            +(((p_data(p)%zfp(l+2)-p_data(p)%zfp(l+1))*dzdsp(l))/(ds**3.))
    awx(l)= (((p_data(p)%xfp(l)-p_data(p)%xfp(l-1))*dxdsp(l))/(ds**3.)) &
            +(((p_data(p)%yfp(l)-p_data(p)%yfp(l-1))*dydsp(l))/(ds**3.)) &
            +(((p_data(p)%zfp(l)-p_data(p)%zfp(l-1))*dzdsp(l))/(ds**3.)) 
  endif
enddo
  



do l=1,NL-1
  if (l==1) then
    if (fboundary=="h") then

      st(l)= d2dt2xs(l)-(dacc(l)*abs(inertia-1))
      st(l)= st(l)-((d2xdsdt(l)**2.)+(d2ydsdt(l)**2.)+(d2zdsdt(l)**2.))
      st(l)=st(l)+(dxdsp(l)*((gammafilament*0.5*(d5xds5(l+1)+d5xds5(l)))+dfxds(l))) &
            +(dydsp(l)*((gammafilament*0.5*(d5yds5(l+1)+d5yds5(l)))+dfyds(l)))&
            +(dzdsp(l)*((gammafilament*0.5*(d5zds5(l+1)+d5zds5(l)))+dfzds(l)))
      st(l)=st(l)+((ap(p)%fxl(l)/ds)*dxdsp(l))+((ap(p)%fyl(l)/ds)*dydsp(l))+((ap(p)%fzl(l)/ds)*dzdsp(l))

    elseif (fboundary=="f") then

      st(l)= d2dt2xs(l)-(dacc(l)*abs(inertia-1))
      st(l)= st(l)-((d2xdsdt(l)**2.)+(d2ydsdt(l)**2.)+(d2zdsdt(l)**2.))
      st(l)=st(l)+(dxdsp(l)*((gammafilament*0.5*(d5xds5(l+1)+d5xds5(l)))+dfxds(l))) &
            +(dydsp(l)*((gammafilament*0.5*(d5yds5(l+1)+d5yds5(l)))+dfyds(l)))&
            +(dzdsp(l)*((gammafilament*0.5*(d5zds5(l+1)+d5zds5(l)))+dfzds(l)))

    elseif (fboundary=="c") then

      st(l)= d2dt2xs(l)-(dacc(l)*abs(inertia-1))
      st(l)= st(l)-((d2xdsdt(l)**2.)+(d2ydsdt(l)**2.)+(d2zdsdt(l)**2.))
      st(l)=st(l)+(dxdsp(l)*((gammafilament*0.5*(d5xds5(l+1)+d5xds5(l)))+dfxds(l))) &
            +(dydsp(l)*((gammafilament*0.5*(d5yds5(l+1)+d5yds5(l)))+dfyds(l)))&
            +(dzdsp(l)*((gammafilament*0.5*(d5zds5(l+1)+d5zds5(l)))+dfzds(l)))
      st(l)=st(l)+((ap(p)%fxl(l)/ds)*dxdsp(l))+((ap(p)%fyl(l)/ds)*dydsp(l))+((ap(p)%fzl(l)/ds)*dzdsp(l))
      st(l)=st(l)+(((d4xds4(l+1)+d4xds4(l))*gammafilament/ds)*dxdsp(l)) &
                  +(((d4yds4(l+1)+d4yds4(l))*gammafilament/ds)*dydsp(l)) &
                  +(((d4zds4(l+1)+d4zds4(l))*gammafilament/ds)*dzdsp(l))

    endif

  else
    st(l)= d2dt2xs(l)-(dacc(l)*abs(inertia-1))   
    st(l)= st(l)-((d2xdsdt(l)**2.)+(d2ydsdt(l)**2.)+(d2zdsdt(l)**2.))
    st(l)=st(l)+(dxdsp(l)*((gammafilament*0.5*(d5xds5(l+1)+d5xds5(l)))+dfxds(l))) &
          +(dydsp(l)*((gammafilament*0.5*(d5yds5(l+1)+d5yds5(l)))+dfyds(l)))&
          +(dzdsp(l)*((gammafilament*0.5*(d5zds5(l+1)+d5zds5(l)))+dfzds(l)))
    st(l)=st(l)-(dxdsp(l)*dfcxds(l))-(dydsp(l)*dfcyds(l))-(dzdsp(l)*dfczds(l))
  endif
enddo



do l=1,NL-1
  if (l==1) then
    ft(l,l)=apx(l)
    ft(l,l+1)=aex(l)
    bt(l)=st(l)
  elseif (l==NL-1) then
    ft(l,l)=apx(l)
    ft(l,l-1)=awx(l)
    bt(l)=st(l)
  else
    ft(l,l)=apx(l)
    ft(l,l+1)=aex(l)
    ft(l,l-1)=awx(l)
    bt(l)=st(l)
  endif
enddo

do l=2,NL-1
  ft(l,l)=(ft(l-1,l)*(-(ft(l,l-1)/ft(l-1,l-1))))+ft(l,l)
  bt(l)=(bt(l-1)*(-(ft(l,l-1)/ft(l-1,l-1))))+bt(l)
  ft(l,l-1)=0
enddo

do l=NL-2,1,-1
  bt(l)=(bt(l+1)*(-(ft(l,l+1)/ft(l+1,l+1))))+bt(l)
  ft(l,l+1)=0
enddo

do l=1,NL-1
  Tension(l)=bt(l)/ft(l,l)   
enddo

!!*********************************inertialess filament***********************************
!elseif (inertia==0) then

!do l=1,NL-1
!dsc(l)=sqrt((dxds(l)**2.)+(dyds(l)**2.)+(dzds(l)**2.))!length of each element scaled with ds
!tension(l)=kappa*(dsc(l)-1.)
!enddo

!endif!if inertia

return
end subroutine comp_tension
end module mod_comp_tension
