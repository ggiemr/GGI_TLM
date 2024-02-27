! Plot the modal fields and also use the fields to calculate the 
! power carried by the mode, the dissipated power and hence the attenuation
!

! SET PLOT RANGE

h0=h                    ! centre of field calculation grid

xmin=h0-xrange_factor*h
xmax=h0+xrange_factor*h

ymin=-yrange_factor*h
ymax= yrange_factor*h

dx=(xmax-xmin)/dble(nx)
dy=(ymax-ymin)/dble(ny)

! END OF RANGE SETTING

xp0=0
dp0=dp0min

! set x values for field calculation

  xp0=0
  yp0=0         ! this is known due to the field symmetry
  dp0min=1e30

  dxmin=(xmax-xmin)/(2d0*nx*(1d0+(dxdy_range_factor-1d0)/2d0))
  dxmax=dxdy_range_factor*dxmin

  xp_field(0)=h0
  do i=1,nx
    dx=dxmin+dble(i-1)*(dxmax-dxmin)/dble(nx)
    xp_field(i) =xp_field(i-1)+dx    
    xp_field(-i)=2d0*xp_field(0)-xp_field(i)  
    
    dp0=abs(xp_field(i))
    if (dp0.LT.dp0min) then
      dp0min=dp0
      xp0=i
    end if
    
    dp0=abs(xp_field(-i))
    if (dp0.LT.dp0min) then
      dp0min=dp0
      xp0=-i
    end if
    
  end do
  
  do i=-nx,nx
! this could be improved...      
    if (i.EQ.-nx) then
      dxp_integral_field(i)=(xp_field(i+1)-xp_field(i))/2d0
    else if (i.EQ.nx) then
      dxp_integral_field(i)=(xp_field(i)-xp_field(i-1))/2d0
    else
      dxp_integral_field(i)=(xp_field(i+1)-xp_field(i-1))/2d0
    end if
  end do

! set y values for field calculation

  dymin=ymax/(ny*(1d0+(dxdy_range_factor-1d0)/2d0))
  dymax=dxdy_range_factor*dymin

  yp_field(0)=0d0
  do i=1,ny
    dy=dymin+dble(i-1)*(dymax-dymin)/dble(ny)
    yp_field(i) =yp_field(i-1)+dy    
    yp_field(-i)=2d0*yp_field(0)-yp_field(i)  
  end do
  
  do i=-ny,ny
! this could be improved...      
    if (i.EQ.-ny) then
      dyp_integral_field(i)=(yp_field(i+1)-yp_field(i))/2d0
    else if (i.EQ.ny) then
      dyp_integral_field(i)=(yp_field(i)-yp_field(i-1))/2d0
    else
      dyp_integral_field(i)=(yp_field(i+1)-yp_field(i-1))/2d0
    end if
  end do
  
  actual_xmin=xp_field(-nx)
  actual_xmax=xp_field(nx)
  actual_ymin=yp_field(-ny)
  actual_ymax=yp_field(ny)
  
  write(*,*)'Field evaluation grid:'
  write(*,*)'dxmin       =',real(dxmin),' dxmax       =',real(dxmax) 
  write(*,*)'xmin       =',real(xmin),' xmax       =',real(xmax) 
  write(*,*)'actual xmin=',real(actual_xmin),' actual xmax=',real(actual_xmax) 
  write(*,*)'dymin       =',real(dymin),' dymax       =',real(dymax) 
  write(*,*)'ymin       =',real(ymin),' ymax       =',real(ymax) 
  write(*,*)'actual ymin=',real(actual_ymin),' actual ymax=',real(actual_ymax) 
 
!! ************** TEMP ***********************
!  open(unit=55,file='xp.dat')
!  do ix=-nx,nx
!    xp=xp_field(ix)
!    dx=dxp_integral_field(ix)
!    write(55,*)xp,dx    
!  end do
!  close(unit=55)
!  open(unit=55,file='yp.dat')
!
!  do iy=-ny,ny  
!    yp=yp_field(iy)
!    dy=dyp_integral_field(iy)
!    write(55,*)yp,dy
!  end do
!  close(unit=55)
!! ************** TEMP ***********************

write(*,*)'Origin point for normalisation, xp0=',xp0,' yp0=',yp0
write(*,*)'coordinates:',xp_field(xp0),yp_field(yp0)

! position in space for Maxwell check and modified bessel function check
mcx = 1.4d0*h
mcy = 0.4d0*h    

! x y cell for Maxwell check
imcx=NINT(mcx-xmin)/dx
imcy=NINT(mcy-ymin)/dy 

! Calculate the unknowns R,T,N,M for the field calculation

if (plot_RTMN) then
  open(unit=80,file='Rl.dat')
  open(unit=81,file='Tl.dat')
  open(unit=82,file='Ml.dat')
  open(unit=83,file='Nl.dat')
end if

do i=-nlambda,nlambda

  l=lambda_field(i)

  u1l(i)=sqrt(l*l+beta*beta-k1*k1)
  if(dble(u1l(i)).LT.0d0) u1l(i)=u1l(i)*(-1d0)   ! sign check. Olsen paper section III
  
  u2l(i)=sqrt(l*l+beta*beta-k2*k2)
  if(dble(u2l(i)).LT.0d0) u2l(i)=u2l(i)*(-1d0)   ! sign check. Olsen paper section III

! Wait, equation 14
  Rl(i)=-1d0+(2d0*k1*k1/(k1*k1-beta*beta))*(l*l-u1l(i)*u2l(i))*u1l(i)/(k1*k1*u2l(i)+k2*k2*u1l(i))

! Wait, equation 13- gives the same solution as equation 14.
!  K13=(k1*k1-beta*beta)/(k2*k2-beta*beta)
!  Rl(i)=( l*l*beta*beta*(1d0-K13**3)+(eps0*w*u1l(i)-eps2*w*u2l(i)*K13)*(mu0*w*u1l(i)+mu0*w*u2l(i)*K13) )/ &
!        (-l*l*beta*beta*(1d0-K13**3)+(eps0*w*u1l(i)+eps2*w*u2l(i)*K13)*(mu0*w*u1l(i)+mu0*w*u2l(i)*K13) )
  
! Wait, equation 9
  Tl(i)=(k1*k1-beta*beta)*(1d0+Rl(i))/(k2*k2-beta*beta)

! Note that N and M are odd functions of l (lambda)
! Wait, equation 11 and 10
  Nl(i)=j*l*beta*(1d0+Rl(i)-Tl(i))/(w*mu0*(u1l(i)*(k2*k2-beta*beta)/(k1*k1-beta*beta)+u2l(i)))

! Wait, equation 10
  Ml(i)=Nl(i)*(k2*k2-beta*beta)/(k1*k1-beta*beta)
  
  if (plot_RTMN) then
    write(80,'(3ES12.4)')l,real(Rl(i)),imag(Rl(i))
    write(81,'(3ES12.4)')l,real(Tl(i)),imag(Tl(i))
    write(82,'(3ES12.4)')l,real(Ml(i)),imag(Ml(i))
    write(83,'(3ES12.4)')l,real(Nl(i)),imag(Nl(i))
  end if
  
end do

if (plot_RTMN) then
  close(unit=80)
  close(unit=81)
  close(unit=82)
  close(unit=83)
end if  

! ********** TEMP FOR TESTING ***********
!INCLUDE "check_modified_bessel.F90"

Ex(:,:)=(0d0,0d0)
Ey(:,:)=(0d0,0d0)
Ez(:,:)=(0d0,0d0)
Hx(:,:)=(0d0,0d0)
Hy(:,:)=(0d0,0d0)
Hz(:,:)=(0d0,0d0)

write(ch,'(I1)')n_modes

area_check=0d0

! loop over the spatial grid

do ix=-nx,nx

  xp=xp_field(ix)
  dx=dxp_integral_field(ix)

  do iy=-ny,ny
  
   yp=yp_field(iy)
   dy=dyp_integral_field(iy)
    
   if ( (ix.NE.0d0).OR.(iy.NE.0) ) then
    
    icount=0
    do i=-nlambda,nlambda

      l=lambda_field(i)

      dl=dl_integral_field(i)
      
      if (xp.LT.0d0) then 
! REGION 2, below the ground  ! no waves in +x direction
        
        Pe_px_my=( Tl(i)*exp(-u1l(i)*h)*exp(u2l(i)*xp) )*(exp(-j*l*yp)/u1l(i))*dl                
        Pm_px_my=( Nl(i)*exp(-u1l(i)*h)*exp(u2l(i)*xp) )*(exp(-j*l*yp)/u1l(i))*dl
        
 ! Wait, equation 1                                            
        Ex(ix,iy)= Ex(ix,iy)-j*beta*u2l(i)*(Pe_px_my)        &
                            -j*w*mu0*( (-j*l)*Pm_px_my ) 
        
        Ey(ix,iy)= Ey(ix,iy)-j*beta*((-j*l)*Pe_px_my)    &
                            +j*w*mu0*u2l(i)*(Pm_px_my )
        
        Ez(ix,iy)= Ez(ix,iy)+Pe_px_my*(k2*k2-beta*beta)   
        
 ! Wait, equation 2              
        Hx(ix,iy)= Hx(ix,iy)-j*beta*u2l(i)*(Pm_px_my)        &
                            +j*w*eps2*( (-j*l)*Pe_px_my ) 
        
        Hy(ix,iy)= Hy(ix,iy)-j*beta*((-j*l)*Pm_px_my)    &
                            -j*w*eps2*u2l(i)*(Pe_px_my)
        
        Hz(ix,iy)= Hz(ix,iy)+(Pm_px_my)*(k2*k2-beta*beta)   
        
        Field_Dx=Ex(ix,iy)*eps2
    
      else if (xp.LT.h) then 
! REGION 1A between the ground and the wire    

        if (Bessel_function_in_field_calc) then  
          Pe_px_my=(0d0,0d0)                                                    ! subtract K0r  (direct) term
          Pe_mx_my=( (1d0+Rl(i))*exp(-u1l(i)*(xp+h)) )*(exp(-j*l*yp)/u1l(i))*dl ! subtract K0r2 (image) term from R
        else
          Pe_px_my=(       exp( u1l(i)*(xp-h)) )*(exp(-j*l*yp)/u1l(i))*dl  
          Pe_mx_my=( Rl(i)*exp(-u1l(i)*(xp+h)) )*(exp(-j*l*yp)/u1l(i))*dl  
        end if
        
        Pm_mx_my=( Ml(i)*exp(-u1l(i)*(xp+h)) )*(exp(-j*l*yp)/u1l(i))*dl
        
 ! Wait, equation 1                
        Ex(ix,iy)= Ex(ix,iy)-j*beta*  u1l(i) *(Pe_px_my)        &
                            -j*beta*(-u1l(i))*(Pe_mx_my)        &
                            -j*w*mu0*( (-j*l)*Pm_mx_my ) 
        
        Ey(ix,iy)= Ey(ix,iy)-j*beta*(-j*l)*(Pe_px_my+Pe_mx_my)+ &
                            +j*w*mu0*(-u1l(i))*(Pm_mx_my )
        
        Ez(ix,iy)= Ez(ix,iy)+(Pe_px_my+Pe_mx_my)*(k1*k1-beta*beta)  
        
 ! Wait, equation 2              
        Hx(ix,iy)= Hx(ix,iy)-j*beta*(-u1l(i))*(Pm_mx_my)      &
                            +j*w*eps0*( (-j*l)*(Pe_px_my+Pe_mx_my) ) 
        
        Hy(ix,iy)= Hy(ix,iy)-j*beta*((-j*l)*Pm_mx_my)    &
                            -j*w*eps0*(-u1l(i))*(Pe_mx_my )  &
                            -j*w*eps0*(+u1l(i))*(Pe_px_my )
        
        Hz(ix,iy)= Hz(ix,iy)+(Pm_mx_my)*(k1*k1-beta*beta)   
        
        Field_Dx=Ex(ix,iy)*eps0
                
      else
! REGION 1B above the wire       ! no waves in -x direction
        
! Need to check signs here... 
        if (Bessel_function_in_field_calc) then
          Pe_mx_my=( (1d0+Rl(i))*exp(-u1l(i)*(xp+h)) )*(exp(-j*l*yp)/u1l(i))*dl ! subtract Kor and K0r2 (image) term from R
        else
          Pe_mx_my=( exp(-u1l(i)*(xp-h))+Rl(i)*exp(-u1l(i)*(xp+h)) )*(exp(-j*l*yp)/u1l(i))*dl
        end if
        Pm_mx_my=( Ml(i)*exp(-u1l(i)*(xp+h)) )*(exp(-j*l*yp)/u1l(i))*dl
        
 ! Wait, equation 1                
        Ex(ix,iy)= Ex(ix,iy)-j*beta*(-u1l(i))*(Pe_mx_my)        &
                            -j*w*mu0*( (-j*l)*Pm_mx_my ) 
        
        Ey(ix,iy)= Ey(ix,iy)-j*beta*(-j*l)*(Pe_mx_my)+          &
                            +j*w*mu0*(-u1l(i))*(Pm_mx_my )
        
        Ez(ix,iy)= Ez(ix,iy)+(Pe_mx_my)*(k1*k1-beta*beta)  
         
 ! Wait, equation 2              
        Hx(ix,iy)= Hx(ix,iy)-j*beta*(-u1l(i))*(Pm_mx_my)      &
                            +j*w*eps0*( (-j*l)*(Pe_mx_my) )
        
        Hy(ix,iy)= Hy(ix,iy)-j*beta*((-j*l)*Pm_mx_my)  &
                            -j*w*eps0*(-u1l(i))*(Pe_mx_my ) 
       
        Hz(ix,iy)= Hz(ix,iy)+(Pm_mx_my)*(k1*k1-beta*beta)   

        Field_Dx=Ex(ix,iy)*eps0
            
      end if ! region
            
    end do ! next l (lambda)
      
    if ((Bessel_function_in_field_calc).AND.(xp.GT.0d0)) then        ! in the air region, add in K0r contributions
              
      kr=sqrt(k1*k1-beta*beta)
      if(imag(kr).GT.0d0) kr=kr*(-1d0)   ! sign check for lossy propagation away from the wire
      r=sqrt( (xp-h)**2+yp**2)

      arg=(0d0,1d0)*r*kr ! argument of modified Bessel function
      CALL calc_modified_Bessel_function_and_derivative((0d0,0d0),arg,K0r,dK0r)
 
      partial_K0r_x=dK0r*(0d0,1d0)*kr*(xp-h)/r
      partial_K0r_y=dK0r*(0d0,1d0)*kr*yp/r
              
      r=sqrt( (xp+h)**2+yp**2)
      arg=(0d0,1d0)*r*kr ! argument of modified Bessel function
      CALL calc_modified_Bessel_function_and_derivative((0d0,0d0),arg,K0r2,dK0r2)

      partial_K0r2_x=dK0r2*(0d0,1d0)*kr*(xp+h)/r
      partial_K0r2_y=dK0r2*(0d0,1d0)*kr*yp/r
      
      Ex(ix,iy)=Ex(ix,iy)-j*beta*2d0*(partial_K0r_x-partial_K0r2_x)
      
      Ey(ix,iy)=Ey(ix,iy)-j*beta*2d0*(partial_K0r_y-partial_K0r2_y)
      
      Ez(ix,iy)=Ez(ix,iy)+2d0*(K0r-K0r2)*(k1*k1-beta*beta)  
      
      Hx(ix,iy)=Hx(ix,iy)+j*w*eps0*2d0*(partial_K0r_y-partial_K0r2_y)
      
      Hy(ix,iy)=Hy(ix,iy)-j*w*eps0*2d0*(partial_K0r_x-partial_K0r2_x)
      
!     Hz has no contribution from the primary and image fields
      
    end if
    
    area_check=area_check+dx*dy
    
  end if  ! not origin point on wire
    
 end do ! next y
    
end do ! next x

write(*,*)'Area check:',area_check
write(*,*)'Should be :',(actual_xmax-actual_xmin)*(actual_ymax-actual_ymin)
write(*,*)'ratio     :',area_check/((actual_xmax-actual_xmin)*(actual_ymax-actual_ymin))

! ************ NEED TO SORT OUT NORM VALUE **********

norm=Ez(xp0,yp0)

write(*,*)'|Field normalisation value|=',real(abs(norm))

Ex(:,:)= Ex(:,:)/norm
Ey(:,:)= Ey(:,:)/norm
Ez(:,:)= Ez(:,:)/norm
Hx(:,:)= Hx(:,:)/norm
Hy(:,:)= Hy(:,:)/norm
Hz(:,:)= Hz(:,:)/norm

if (show_field_plots) then

  open(unit=50,file='Ex.dat'//ch)
  open(unit=51,file='Ey.dat'//ch)
  open(unit=52,file='Ez.dat'//ch)

  open(unit=60,file='Hx.dat'//ch)
  open(unit=61,file='Hy.dat'//ch)
  open(unit=62,file='Hz.dat'//ch)

! loop over the spatial grid

  do ix=-nx,nx

    xp=xp_field(ix)
    dx=dxp_integral_field(ix)

    do iy=-ny,ny
  
    yp=yp_field(iy)
    dy=dyp_integral_field(iy)
  
      if ( (mod(ix,plot_everyx).Eq.0).AND.(mod(iy,plot_everyy).Eq.0) ) then
        write(50,'(4ES12.4)')xp,yp,real(Ex(ix,iy)),imag(Ex(ix,iy))
        write(51,'(4ES12.4)')xp,yp,real(Ey(ix,iy)),imag(Ey(ix,iy))
        write(52,'(4ES12.4)')xp,yp,real(Ez(ix,iy)),imag(Ez(ix,iy))
  
        write(60,'(4ES12.4)')xp,yp,real(Hx(ix,iy)),imag(Hx(ix,iy))
        write(61,'(4ES12.4)')xp,yp,real(Hy(ix,iy)),imag(Hy(ix,iy))
        write(62,'(4ES12.4)')xp,yp,real(Hz(ix,iy)),imag(Hz(ix,iy))
      end if
    
    end do ! next y
  
    if (mod(ix,plot_everyx).Eq.0) then
      write(50,*)
      write(51,*)
      write(52,*)
  
      write(60,*)
      write(61,*)
      write(62,*)
    end if
 
  end do ! next x

  close(unit=50)
  close(unit=51)
  close(unit=52)

  close(unit=60)
  close(unit=61)
  close(unit=62)
  
end if

if (maxwell_check) then
! Numerical check that Maxwell's equations are satisfied by the calculated fields
  
  INCLUDE 'Maxwell_check.F90'

end if

INCLUDE 'attenuation_calc.F90'

!STOP 1


  
