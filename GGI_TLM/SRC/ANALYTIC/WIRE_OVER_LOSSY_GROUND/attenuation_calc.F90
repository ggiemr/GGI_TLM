
! loop over the field mesh to calculate the power carried by the mode
! and the power lost and hence the mode attenuation

Pec=(0d0,0d0)
Pzc=(0d0,0d0)

! loop over the spatial grid to calculate the propagating power
do ix=-nx,nx

  xp=xp_field(ix)
  dx=dxp_integral_field(ix)

  do iy=-ny,ny
  
   yp=yp_field(iy)
   dy=dyp_integral_field(iy)

    Pzc=Pzc+0.5d0*( Ex(ix,iy)*conjg(Hy(ix,iy)) - Ey(ix,iy)*conjg(Hx(ix,iy)) )*dx*dy
    
  end do ! next y
 
end do  ! next x

! work out an analytic contribution from the central cell containing the wire
Hwire=(Hy(1,0)+Hx(0,-1)-Hy(-1,0)-Hx(0,1))/4d0
r_Hwire=sqrt(  ( (xp_field(1)-xp_field(-1))/2d0 )**2+ ( (yp_field(1)-yp_field(-1))/2d0 )**2  )/sqrt(2d0)
I_wire=Hwire*2d0*pi*r_Hwire

write(*,*)'r_Hwire=',r_Hwire
write(*,*)'Hwire=',Hwire
write(*,*)'Iwire=',I_wire

rs=r_Hwire
rw=a

C_wire=2.0*pi*eps0/log(rs/rw)
L_wire=(mu0/(2.0*pi))*log(rs/rw)
Z0_wire=sqrt(L_wire/C_wire)

write(*,*)'Z0_wire=',Z0_wire
 
Pzwc=0.5*abs(I_wire**2)*Z0_wire

eps_ratio=eps0/eps2

! loop over the spatial grid to calculate the dissipated power
do ix=-nx,nx

  xp=xp_field(ix)
  dx=dxp_integral_field(ix)

! limits of this integration cell in x  
  xpp=xp+dx/2d0
  xpm=xp-dx/2d0
  
  if (xpm.LT.0d0) then    ! there should be some contribution from the lossy ground material
    
    if (xpp.GT.0d0) then  ! the contribution is smaller than dx
      dx=0d0-xpm
    end if

    do iy=-ny,ny
  
     yp=yp_field(iy)
     dy=dyp_integral_field(iy)
    
      if (xp.GE.0d0) then  ! the field evaluation point is in the air so we need to recognise that Ex is discontinuous

       Ex_in=Ex(ix,iy)*eps_ratio
       Pec=Pec+0.5d0*sigma*(Ex_in*conjg(Ex_in) + Ey(ix,iy)*conjg(Ey(ix,iy)) +Ez(ix,iy)*conjg(Ez(ix,iy)) )*dx*dy

      else
   
       Pec=Pec+0.5d0*sigma*(Ex(ix,iy)*conjg(Ex(ix,iy)) + Ey(ix,iy)*conjg(Ey(ix,iy)) +Ez(ix,iy)*conjg(Ez(ix,iy)) )*dx*dy
      
      end if
 
    end do ! next y
 
  end if
 
end do  ! next x

Pe=dble(Pec)
Pz=dble(Pzc)
Pzw=dble(Pzwc)

alpha_f=Pe/(2d0*(Pz+Pzw))

write(*,*)
write(*,*)'Attenuation calculation from fields'
write(*,*)'lmax=kxmax=',lmax
write(*,*)'lambda_x max=',6.28315d0/lmax
write(*,*)'dl=dkx=',dl
write(*,*)'lambda_x min=',6.28315d0/dl
write(*,*)'Complex Poynting vector:',Pzc
write(*,*)
write(*,*)'***Power carried by mode field *** :',real(Pz)
write(*,*)'***Power carried close to wire***  :',real(Pzw)
write(*,*)'***Dissipated power     ***  :',real(Pe)
write(*,*)
write(*,*)'***Attenuation constant from field***   :',alpha_f
write(*,*)'***Attenuation over specified propagation distance***   :',&
                                                  real(20d0*log10( exp(-alpha_f*prop_dist) )),' dB'
write(*,*)'***Attenuation over 5km***   :',&
                                                  real(20d0*log10( exp(-alpha_f*5000.0) )),' dB'

!write(*,*)'beta=',cmplx(beta),' 9km attenuation=',real(20*log10(abs(exp(-gammaw*prop_dist)))),'dB'

write(*,*)'_________________________________________________________________________________'
write(*,'(A)')'    nx    ny  xrng  yrng   nl    lrng    Pz          Pe         alpha  attn(9km)'
write(*,'(2I6,2F6.1,I6,F6.1,2ES12.3,F10.5,F8.1)')nx,ny,xrange_factor,yrange_factor,&
          nlambda,lrange_factor,real(Pz),real(Pe),real(alpha_f),real(20*log10(abs(exp(-alpha_f*prop_dist))))
write(*,*)'_________________________________________________________________________________'
write(*,*)
  
beta_f=dble(beta)
write(26,'(3ES12.4)',ADVANCE='NO')1000.0*real(alpha_f),real(beta_f),20*log10(abs(exp(-alpha_f*prop_dist)))

