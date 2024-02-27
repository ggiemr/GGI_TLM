
dz=sqrt(dx*dx+dy*dy)

ejbetaz=exp(-j*beta*dz)

if (mcx.GT.0d0) then
  mc_epsilon=eps0     ! in air
else  
  mc_epsilon=eps2     ! in ground plane
end if
! Get fields on local grid

do ix=-1,1
  do iy=-1,1
  
    Exmc(ix,iy,0)=Ex(imcx+ix,imcy+iy)
    Exmc(ix,iy,1)=Ex(imcx+ix,imcy+iy)*ejbetaz  ! forward propagation
    Exmc(ix,iy,-1)=Ex(imcx+ix,imcy+iy)/ejbetaz  ! backward propagation
  
    Eymc(ix,iy,0)=Ey(imcx+ix,imcy+iy)
    Eymc(ix,iy,1)=Ey(imcx+ix,imcy+iy)*ejbetaz  ! forward propagation
    Eymc(ix,iy,-1)=Ey(imcx+ix,imcy+iy)/ejbetaz  ! backward propagation
  
    Ezmc(ix,iy,0)=Ez(imcx+ix,imcy+iy)
    Ezmc(ix,iy,1)=Ez(imcx+ix,imcy+iy)*ejbetaz  ! forward propagation
    Ezmc(ix,iy,-1)=Ez(imcx+ix,imcy+iy)/ejbetaz  ! backward propagation
  
    Hxmc(ix,iy,0)=Hx(imcx+ix,imcy+iy)
    Hxmc(ix,iy,1)=Hx(imcx+ix,imcy+iy)*ejbetaz  ! forward propagation
    Hxmc(ix,iy,-1)=Hx(imcx+ix,imcy+iy)/ejbetaz  ! backward propagation
  
    Hymc(ix,iy,0)=Hy(imcx+ix,imcy+iy)
    Hymc(ix,iy,1)=Hy(imcx+ix,imcy+iy)*ejbetaz  ! forward propagation
    Hymc(ix,iy,-1)=Hy(imcx+ix,imcy+iy)/ejbetaz  ! backward propagation
  
    Hzmc(ix,iy,0)=Hz(imcx+ix,imcy+iy)
    Hzmc(ix,iy,1)=Hz(imcx+ix,imcy+iy)*ejbetaz  ! forward propagation
    Hzmc(ix,iy,-1)=Hz(imcx+ix,imcy+iy)/ejbetaz  ! backward propagation

  end do
end do

! fields at centre

Exc=Exmc(0,0,0)
Eyc=Eymc(0,0,0)
Ezc=Ezmc(0,0,0)

Hxc=Hxmc(0,0,0)
Hyc=Hymc(0,0,0)
Hzc=Hzmc(0,0,0)

! E field derivatives

dEx_dy=(Exmc(0,1,0)-Exmc(0,-1,0))/(dy*2d0)
dEx_dz=(Exmc(0,0,1)-Exmc(0,0,-1))/(dz*2d0)

dEy_dx=(Eymc(1,0,0)-Eymc(-1,0,0))/(dx*2d0)
dEy_dz=(Eymc(0,0,1)-Eymc(0,0,-1))/(dz*2d0)

dEz_dx=(Ezmc(1,0,0)-Ezmc(-1,0,0))/(dx*2d0)
dEz_dy=(Ezmc(0,1,0)-Ezmc(0,-1,0))/(dy*2d0)

! H field derivatives

dHx_dy=(Hxmc(0,1,0)-Hxmc(0,-1,0))/(dy*2d0)
dHx_dz=(Hxmc(0,0,1)-Hxmc(0,0,-1))/(dz*2d0)

dHy_dx=(Hymc(1,0,0)-Hymc(-1,0,0))/(dx*2d0)
dHy_dz=(Hymc(0,0,1)-Hymc(0,0,-1))/(dz*2d0)

dHz_dx=(Hzmc(1,0,0)-Hzmc(-1,0,0))/(dx*2d0)
dHz_dy=(Hzmc(0,1,0)-Hzmc(0,-1,0))/(dy*2d0)

! Check components of curl H=epsilon_0 dE/dt

write(*,*)'Check components of curl H=epsilon dE/dt'

LHS=dHz_dy-dHy_dz
RHS=mc_epsilon*j*w*Exc

write(*,*)'x:',LHS
write(*,*)'x:',RHS
write(*,*)'  ',LHS/RHS,abs(LHS/RHS)
write(*,*)

LHS=dHx_dz-dHz_dx
RHS=mc_epsilon*j*w*Eyc

write(*,*)'y:',LHS
write(*,*)'y:',RHS
write(*,*)'  ',LHS/RHS,abs(LHS/RHS)
write(*,*)

LHS=dHy_dx-dHx_dy
RHS=mc_epsilon*j*w*Ezc

write(*,*)'z:',LHS
write(*,*)'z:',RHS
write(*,*)'  ',LHS/RHS,abs(LHS/RHS)
write(*,*)

! Check components of curl E=-mu_0 dH/dt

write(*,*)'Check components of curl E=-mu_0 dH/dt'

LHS=dEz_dy-dEy_dz
RHS=-mu0*j*w*Hxc

write(*,*)'x:',LHS
write(*,*)'x:',RHS
write(*,*)'  ',LHS/RHS,abs(LHS/RHS)
write(*,*)

LHS=dEx_dz-dEz_dx
RHS=-mu0*j*w*Hyc

write(*,*)'y:',LHS
write(*,*)'y:',RHS
write(*,*)'  ',LHS/RHS,abs(LHS/RHS)
write(*,*)

LHS=dEy_dx-dEx_dy
RHS=-mu0*j*w*Hzc

write(*,*)'z:',LHS
write(*,*)'z:',RHS
write(*,*)'  ',LHS/RHS,abs(LHS/RHS)
write(*,*)
