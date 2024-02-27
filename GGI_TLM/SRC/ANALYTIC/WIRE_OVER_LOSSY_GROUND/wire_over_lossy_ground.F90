
PROGRAM propagation_model

IMPLICIT NONE

! ****************************************************************************
! ****************************************************************************
! ****************************************************************************

! The following parameters and flags can be used to provide more information
! from the solution but do require the software to be recompiled
! In order to do this, set the flag run_from_internal_setup=.TRUE. to run from an internal setup.
! You can then specify the geometry, field calculation parameters etc using the parameters below.

logical :: run_from_internal_setup=.FALSE.  

! Field calculation parameters
real*8,parameter  :: xrange_factor=20d0 !  10 multiplier for h to give the xmin and xmax values
real*8,parameter  :: yrange_factor=40d0 !  20 multiplier for h to give the ymin and ymax values
real*8,parameter  :: lrange_factor=10d0   ! multiplier for lambda to give lmax value
real*8,parameter  :: dlrange_factor=10d0  ! range of delta_lambda
integer,parameter :: nlambda=1000
integer,parameter :: nx=50          !  50 full x range is -nx to nx
integer,parameter :: ny=100         ! 100 full y range is -ny to ny
integer,parameter :: plot_everyx=4 !4 
integer,parameter :: plot_everyy=4 !4
real*8,parameter  :: dxdy_range_factor=8.0d0   ! range of dx and dy in field calculation

logical :: plot_integrand=.FALSE.
logical :: show_plots=.FALSE.
logical :: calculate_fields=.FALSE.
logical :: show_field_plots=.FALSE.
logical :: maxwell_check=.FALSE.
logical :: plot_RTMN=.FALSE.
logical :: Bessel_function_in_field_calc=.TRUE.  ! Should not be changed

! ****************************************************************************
! ****************************************************************************
! ****************************************************************************

real*8 :: fmin,fmax,f,w
real*8 :: logfmin,logfmax,logfstep,logf

real*8 :: sigma,a,h,epsr
complex*16 :: eps2

real*8 :: test1,test2,test3,test4,test5

complex*16 :: Z,Y,gamma,beta
complex*16 :: Zc,Yc,gammac,betac
complex*16 :: Zc2,Yc2,gammac2,betac2
complex*16 :: Zw,Yw,gammaw,betaw
complex*16 :: delta,k1,k2,ka,Jc,ALPHA,QmjP,NmjM,K13

integer :: nf,floop

! complex plane search stuff for Wait's solution

integer    :: itest,min_test_point
real*8     :: dbeta_r,dbeta_i
complex*16 :: beta_test(0:4),new_beta,beta_save
real*8     :: error_test(0:4),error,min_error,initial_error,final_error,spatial_avg_final_error

real*8     :: prop_dist

integer    :: iteration
integer    :: max_iterations=1000

! Wait model mode search stuff

integer,parameter :: max_modes=10000

integer :: n_trial_modes,n_modes,mode

complex*16 :: beta_err_min(max_modes)
complex*16 :: beta_found(max_modes)

integer    :: n_re_beta,ire
real*8     :: re_beta_min,re_beta_max,d_re_beta,re_beta

integer    :: n_im_beta,iim
real*8     :: im_beta_min,im_beta_max,d_im_beta,im_beta

real*8     :: beta_diff,min_beta_diff
integer    :: mode_test
logical    :: already_found_this_mode

logical :: plot_beta_map

real*8,allocatable     :: err_array(:,:)
complex*16,allocatable :: beta_array(:,:)

real*8 :: tw_a,tw_r,tw_L,tw_C
complex*16:: Ztlm,Ytlm,Ztot,Ytot
complex*16:: gammatlm

! Field calculation stuff

real*8 :: l,lmax,dl,dlmin,dlmax
real*8 :: fimin
complex*16 :: K0a,K0h,arg
complex*16 :: u1,u2
complex*16 :: integrand,integrand2
complex*16 :: exp_term

real*8     :: r

logical :: Olsen_eqn=.TRUE.
!logical :: Olsen_eqn=.FALSE.
complex*16 :: Kterm,Integral_term

integer :: nlambda_local
integer :: icount

! expansion coefficients which are a function of l (lambda)
real*8     :: lambda_field(-nlambda:nlambda)
real*8     :: dl_integral_field(-nlambda:nlambda)

real*8     :: lambda(0:nlambda)
real*8     :: dl_integral(0:nlambda)

! discretisation in x and y for field evaluation and attenuation calculation
real*8     :: xp_field(-nx:nx)
real*8     :: dxp_integral_field(-nx:nx)
real*8     :: dxmin,dxmax,h0,xpp,xpm

real*8     :: yp_field(-ny:ny)
real*8     :: dyp_integral_field(-ny:ny)
real*8     :: dymin,dymax

real*8     :: actual_xmin,actual_xmax
real*8     :: actual_ymin,actual_ymax

integer :: xp0,yp0
real*8  :: dp0,dp0min
complex*16 :: eps_ratio,Ex_in

complex*16 :: u1l(-nlambda:nlambda),u2l(-nlambda:nlambda)
complex*16 :: Rl(-nlambda:nlambda),Tl(-nlambda:nlambda),Nl(-nlambda:nlambda),Ml(-nlambda:nlambda)
real*8     :: delta_l

real*8     :: area_check

complex*16 :: kx,ky,kz,kr,K0ra,K0r,K0r2,dK0r,dK0r2,K0rp,K0rm

complex*16 :: partial_K0r_x ,partial_K0r_y
complex*16 :: partial_K0r2_x,partial_K0r2_y

complex*16 :: K0r_i
complex*16 :: partial_K0r_x_i ,partial_K0r_y_i

real*8 :: beta_f

complex*16 :: Hwire,I_wire
real*8     :: r_Hwire
real*8     :: Z0_wire,L_wire,C_wire,rs,rw

complex*16 :: Pe_px_my,Pe_mx_my
complex*16 :: Pm_px_my,Pm_mx_my

! complex electromagnetic fields
complex*16 :: Ex(-nx:nx,-ny:ny),Ey(-nx:nx,-ny:ny),Ez(-nx:nx,-ny:ny),Field_Dx,norm
complex*16 :: Hx(-nx:nx,-ny:ny),Hy(-nx:nx,-ny:ny),Hz(-nx:nx,-ny:ny)

complex*16 :: Pec,Pzc,Pzwc
real*8     :: Pe,Pz,Pzw,alpha_f

real*8 :: xmin,xmax,ymin,ymax,dx,dy,dz,xp,yp,rp

character :: ch

! Maxwell check stuff...

real*8 :: mcx,mcy           ! position in space for Maxwell check
integer :: imcx,imcy        ! x y cell for Maxwell check

complex*16 :: Exmc(-1:1,-1:1,-1:1),Eymc(-1:1,-1:1,-1:1),Ezmc(-1:1,-1:1,-1:1)
complex*16 :: Hxmc(-1:1,-1:1,-1:1),Hymc(-1:1,-1:1,-1:1),Hzmc(-1:1,-1:1,-1:1)

integer :: iz

complex*16 :: Exc,Eyc,Ezc,Hxc,Hyc,Hzc
complex*16 :: ejbetaz

complex*16 :: mc_epsilon

complex*16 :: dEx_dy,dEx_dz,dEy_dx,dEy_dz,dEz_dx,dEz_dy
complex*16 :: dHx_dy,dHx_dz,dHy_dx,dHy_dz,dHz_dx,dHz_dy

complex*16 :: LHS,RHS

integer :: i,ix,iy

character(LEN=80) :: FH2_DIRECTORY

! constants

real*8,parameter :: small=1d-12

real*8,parameter :: large=1d12

real*8,parameter :: pi=3.141592653589793d0

real*8,parameter :: Z0=376.73031346177d0

real*8,parameter :: Y0=1D0/Z0

complex*16,parameter :: j=(0d0,1d0)

real*8,parameter :: eps0=8.8541878176D-12

real*8,parameter :: mu0=3.141592653589793d0*4D-7

real*8,parameter :: c0=2.99792458D8


! START

if (run_from_internal_setup) then

! Set up default test case values

  a=0.013      ! wire radius
  h=20.0       ! wire height over ground
  sigma=0.001  ! ground electrical conductivity
  epsr=1.0     ! ground relative permittivity

  fmin=2e6
  fmax=2e6
  nf=1

  prop_dist=9000d0

  n_re_beta=101
  n_im_beta=101
  
  calculate_fields=.TRUE.

else  ! Get the problem specification from user inputs

  write(*,*)'GGI_wire_over_lossy_ground'
  write(*,*)'Calculate the complex propagation constant of modes on a wire over lossy ground'
  write(*,*)'using the following solution methods:'
  write(*,*)'Full method due to Carson'
  write(*,*)'Approximate method due to Carson'
  write(*,*)'Rigorous method due to Wait (multi-mode)'
  write(*,*)'This software uses the Fortran complex bessel function module mod_zbes.f90'
  write(*,*)'downloaded from: dl.acm.org/doi/10.1145/1916461.1916471'
  write(*,*)
  write(*,*)'Enter the wire radius (m)'
  read(*,*)a
  write(*,*)'Enter the wire height over the ground (m)'
  read(*,*)h
  write(*,*)'Enter the ground electrical conductivity (S/m)'
  read(*,*)sigma
  write(*,*)'Enter the ground relative permittivity'
  read(*,*)epsr
  write(*,*)'Enter the minimum frequency for analysis (Hz)'
  read(*,*)fmin
  write(*,*)'Enter the maximum frequency for analysis (Hz)'
  read(*,*)fmax
  write(*,*)'Enter the number of frequencies to evaluate (log frequency scale)'
  read(*,*)nf

  write(*,*)'Writing configuration to file: wire_over_lossy_ground_configuration.dat'

  open(unit=10,file='wire_over_lossy_ground_configuration.dat')

  write(10,*)real(a),'  wire radius (m)'
  write(10,*)real(h),'  wire height over the ground (m)'
  write(10,*)real(sigma),'  ground electrical conductivity (S/m)'
  write(10,*)real(epsr),'  ground relative permittivity'
  write(10,*)real(fmin),'  minimum frequency for analysis'
  write(10,*)real(fmax),'  maximum frequency for analysis'
  write(10,*)nf,'  number of frequencies to evaluate (log frequency scale)'

  close(unit=10)

  n_re_beta=101
  n_im_beta=101

  prop_dist=9000d0
  
end if

!! specification for TLM thin wire model propagation constant estimate 
!tw_a=0.013                         ! TLM thin wire model radius
!tw_r=a                             ! TLM thin wire model reference radius
!tw_L=(mu0/(2.0*pi))*log(tw_r/tw_a)
!tw_C=2.0*pi*eps0/log(tw_r/tw_a)

ALLOCATE( err_array(1:n_re_beta,1:n_im_beta) )
ALLOCATE( beta_array(1:n_re_beta,1:n_im_beta) )

open(unit=10,file='k0.out')
open(unit=12,file='Approx_Carson_propagation_model.out')
open(unit=14,file='Full_Carson_propagation_model.out')
open(unit=16,file='Wait_propagation_model.out')
open(unit=18,file='Wait_propagation_model_all_modes.out')

!open(unit=20,file='TLM_thin_wire_model_estimate.out')  ! removed

if (calculate_fields) then
  open(unit=26,file='Wait_propagation_model_from_fields.out')
end if

logfmin=log10(fmin)
logfmax=log10(fmax)
if (nf.NE.1) then
  logfstep=(logfmax-logfmin)/dble(nf-1)
  plot_beta_map=.FALSE.
else
  logfstep=0d0
  plot_beta_map=.TRUE.
end if

do floop=1,nf

  logf=logfmin+dble(floop-1)*logfstep
  f=10.0**(logf)
  w=2d0*pi*f
  
  write(*,*)
  write(*,*)'Frequency=',real(f)
  
  write(10,'(2ES12.4)')f,w*sqrt(mu0*eps0)
  
! set up the intrgration discretisation for integrals over lambda

! get the discretisation for l (lambda in Wait paper)

  fimin=1d-4

! integration sampling for propagation constant solution
  lmax=-lrange_factor*log(fimin)/(2d0*h)  ! increase range of integral over lambda for field calculation

  dlmin=lmax/(nlambda*(1d0+(dlrange_factor-1d0)/2d0))
  dlmax=dlrange_factor*dlmin
    
  write(*,*)'lmax=',lmax
  write(*,*)'dlrange_factor=',dlrange_factor
  write(*,*)'lmax/nlambda=',lmax/dble(nlambda)
  write(*,*)'dlmin=',dlmin,' dlmax=',dlmax

  lambda(0)=0d0
  do i=1,nlambda
    dl=dlmin+dble(i-1)*(dlmax-dlmin)/dble(nlambda)
    lambda(i)=lambda(i-1)+dl
  end do
  do i=0,nlambda
    l=lambda(i)
! this could be improved...      
    if (i.EQ.0) then
      dl_integral(i)=(lambda(1))/2d0
    else if (i.EQ.nlambda) then
      dl_integral(i)=(lambda(i)-lambda(i-1))/2d0
    else
      dl_integral(i)=(lambda(i+1)-lambda(i-1))/2d0
    end if
  end do

! integration sampling for field evaluation
  lmax=-lrange_factor*log(fimin)/(2d0*h)  ! increase range of integral over lambda for field calculation

  dlmin=lmax/(nlambda*(1d0+(dlrange_factor-1d0)/2d0))
  dlmax=dlrange_factor*dlmin

  write(*,*)'lmax_field=',lmax
  write(*,*)'lmax_field/nlambda=',lmax/dble(nlambda)
  write(*,*)'dlmin_field=',dlmin,' dlmax=',dlmax

  lambda_field(0)=0d0
  do i=1,nlambda
    dl=dlmin+dble(i-1)*(dlmax-dlmin)/dble(nlambda)
    lambda_field(i)=lambda_field(i-1)+dl
    lambda_field(-i)=-lambda_field(i)  
  end do
  do i=-nlambda,nlambda
    l=lambda_field(i)
! this could be improved...      
    if (i.EQ.-nlambda) then
      dl_integral_field(i)=(lambda_field(i+1)-lambda_field(i))/2d0
    else if (i.EQ.nlambda) then
      dl_integral_field(i)=(lambda_field(i)-lambda_field(i-1))/2d0
    else
      dl_integral_field(i)=(lambda_field(i+1)-lambda_field(i-1))/2d0
    end if
  end do
  
  if (plot_RTMN) then
  
    open(unit=70,file='lambda.dat')
    open(unit=71,file='lambda_field.dat')
    do i=0,nlambda
      write(70,'(2ES12.4)')lambda(i),dl_integral(i)
    end do
    do i=-nlambda,nlambda
      write(71,'(2ES12.4)')lambda_field(i),dl_integral_field(i)
    end do
    close(unit=70)
    close(unit=71)
    
  end if

! Approximate Carson model (no magnetic materials)

  delta=sqrt(2.0/(w*mu0*sigma))
  Zc=w*mu0/8.0+(j*w*mu0/(2.0*pi))*( log(sqrt(2.0)*delta/a)+log(0.926) )
  Yc=2.0*pi*j*w*eps0/log(2.0*h/a)    ! admittance for wire over perfect ground plane
  gammac=sqrt(Zc*Yc)
  k2=w*sqrt(mu0*eps0*epsr)
  
  write(*,*)'Zc =',cmplx(Zc),' Yc =',cmplx(Yc),' gammac =',cmplx(gammac)
    
  test1=real(2.0*k2*h)
  betac=gammac/j
  write(12,'(4ES12.4)')f,1000.0*real(gammac),imag(gammac),20*log10(abs(exp(-gammac*prop_dist)))
  
  
! Full Carson model 

  eps2=(eps0*epsr-j*sigma/w)
  
  k1=w*sqrt(mu0*eps0)                    ! this is real
  k2=w*sqrt(mu0*(eps0*epsr-j*sigma/w))   ! this is complex. 
  if(imag(k2).GT.0d0) k2=k2*(-1d0)       ! Imaginary part should be negative (Olsen, last parageraph of section III)
  
  CALL calc_carson_integral(k2,f,a,h,sigma,epsr,Jc,nlambda,lambda,dl_integral)

  Zc2=(j*w*mu0/(2.0*pi))*(log(2.0*h/a)+Jc)
  Yc2=2.0*pi*j*w*eps0/log(2.0*h/a)
  betac2=k1*sqrt(1d0+Jc/log(2.0*h/a))
  
  gammac2=j*betac2
  write(*,*)'Zc2=',cmplx(Zc2),' Yc2=',cmplx(Yc2),' gammac2=',cmplx(gammac2)
  gammac2=sqrt(Zc2*Yc2)  
  write(*,*)'Zc2=',cmplx(Zc2),' Yc2=',cmplx(Yc2),' gammac2=',cmplx(gammac2)
    
  test1=abs(a*sqrt(k1*k1-betac2*betac2))
  test2=abs(2d0*h*sqrt(k1*k1-betac2*betac2))
  test3=a/(2.0*h)
  test4=abs(k1*h)
  test5=abs(k1*k1/(k2*k2))
  
  write(*,*)'--------------------------------------------------'
  write(*,*)'Carson approximation conditions:'
  
  write(*,*)'test1:',test1,' <<1'
  write(*,*)'test2:',test2,' <<1'
  write(*,*)'test3:',test3,' <<1'
  write(*,*)'test4:',test4,' <<1'
  write(*,*)'test5:',test5,' <<1'
  
  write(*,*)'--------------------------------------------------'
  
  write(14,'(4ES12.4)')f,1000.0*real(gammac2),imag(gammac2),20*log10(abs(exp(-gammac2*prop_dist)))
  
!  GOTO 8888  ! uncomment to skip Wait's method
  
! Wait model: starting beta from the Carson model 

  ka=k1*k2/sqrt(k2*k2+k1*k1)
  if(imag(ka).GT.0d0) ka=-ka    ! Olsen, section VII, top of p 1417

  write(*,*)
  write(*,*)'Approx Carson model, beta=',betac
  write(*,*)'Full Carson model,   beta=',betac2
  write(*,*)
  write(*,*)'Branch points:'
  write(*,*)'beta=k0=',k1
  write(*,*)'beta=k2=',k2
  write(*,*)'beta=ka=',ka
  write(*,*)"Calculating resuts of Wait's model:"
  
  err_array(:,:)=0d0
  beta_array(:,:)=(0d0,0d0)
  
! Set search range for beta

  re_beta_min=k1*0.95
  re_beta_max=k1*1.1
  d_re_beta=(re_beta_max-re_beta_min)/dble(n_re_beta-1)
  
  write(*,*)'re_beta_min=',re_beta_min,' re_beta_max=',re_beta_max,' d_re_beta',d_re_beta

  im_beta_min=imag(sqrt(Zc2*Yc2)/j) *5d0                  ! range derived from Carson model
  im_beta_min=min(im_beta_min,-0.5e-3)
  
  im_beta_max=0d0
  d_im_beta=(im_beta_max-im_beta_min)/dble(n_im_beta-1)
  
  write(*,*)'im_beta_min=',im_beta_min,' im_beta_max=',im_beta_max,' d_im_beta',d_im_beta
  
  min_beta_diff=sqrt(d_re_beta**2+d_im_beta**2)/10.0   ! measure for testing if modes are equivalent
  
! evaluate error on grid

  write(*,*)'evaluate error on grid'
  
  if (plot_beta_map) then
    open(unit=50,file='beta_err_map.dat')
  end if
  
  do ire=1,n_re_beta
  
    re_beta=re_beta_min+dble(ire-1)*d_re_beta

    do iim=1,n_im_beta
  
      im_beta=im_beta_min+dble(iim-1)*d_im_beta
      
      beta=dcmplx(re_beta,im_beta)
            
      beta_array(ire,iim)=beta
      
      plot_integrand=.FALSE.
      INCLUDE 'calc_error_Wait.F90'  
        
      err_array(ire,iim)=error
      
      if (plot_beta_map) then
        write(50,'(3ES12.4)')real(beta_array(ire,iim)),imag(beta_array(ire,iim)),real(log10(err_array(ire,iim)))
      end if
      
    end do
    
    if (plot_beta_map) then
      write(50,*)
    end if
    
  end do
  
  if (plot_beta_map) then
    close(unit=50)
  end if

! find local minima in error

  n_trial_modes=0
  
  do iim=2,n_im_beta-1
  
    do ire=2,n_re_beta-1
    
      error=err_array(ire,iim)
          
      if ( (error.LT.err_array(ire+1,iim  )).AND.   &
           (error.LT.err_array(ire-1,iim  )).AND.   &
           (error.LT.err_array(ire  ,iim+1)).AND.   &
           (error.LT.err_array(ire  ,iim-1))      ) then

! found a local minimum indicating a mode

        n_trial_modes=n_trial_modes+1 
        beta_err_min(n_trial_modes)=beta_array(ire,iim)
           
      end if
      
    end do
    
  end do

  write(16,'(ES12.4)',ADVANCE='NO')f
  if (calculate_fields) then
    write(26,'(ES12.4)',ADVANCE='NO')f
  end if

! loop over local minima in error and accurately calculate the mode

  n_modes=0
  
  do mode=1,n_trial_modes

! set starting point for mode search

    beta=beta_err_min(mode)
  
    k1=w*sqrt(mu0*eps0)
    k2=w*sqrt(mu0*(eps0*epsr-j*sigma/w))
    if(imag(k2).GT.0d0) k2=k2*(-1d0)       ! Imaginary part should be negative (Olsen, last parageraph of section III)
  
    beta_test(0)=beta
    
    iteration=0
    initial_error=0d0
  
    dbeta_r=dble(beta)/100d0
    dbeta_i=dble(beta)/100d0
  
5000 CONTINUE   ! start of a new loop of the optimisation process

    iteration=iteration+1

! set the points of the simplex
    beta_test(1)=beta_test(0)-dbeta_r
    beta_test(2)=beta_test(0)+dbeta_r
    beta_test(3)=beta_test(0)-j*dbeta_i
    beta_test(4)=beta_test(0)+j*dbeta_i
  
    min_error=1e30
    min_test_point=-1
  
    do itest=0,4
  
      beta=beta_test(itest)
      plot_integrand=.FALSE.
      INCLUDE 'calc_error_Wait.F90'    
      error_test(itest)=error
          
      if (iteration.EQ.1) then                   ! NEW: include all 5 points in initial error calculation
        initial_error=initial_error+error/5d0
      end if
    
      if (error.LT.min_error) then
        min_error=error
        min_test_point=itest
      end if
  
    end do
  
    if (min_test_point.EQ.-1) then
      write(*,*)'ERROR evaulating error in Wait solution:'
      STOP 1
    end if
  
    if (min_test_point.EQ.0) then
!   reduce the size of the simplex
      dbeta_r=dbeta_r/2d0
      dbeta_i=dbeta_i/2d0
    else
!   put the minimum error solution at the centre of the simplex  
      beta_test(0)=beta_test(min_test_point)
    end if
  
    if (iteration.LT.max_iterations) GOTO 5000
 
! recreate the simplex around the minimum point 
    betaw=beta_test(0)
    
! set the points of the simplex surrounding the minimum point but increase the simplex size to 
! detect minima near branch lines in the solution which are not true solutions

    dbeta_r=dble(beta)/1000d0
    dbeta_i=dble(beta)/1000d0
    
    beta_test(1)=beta_test(0)-dbeta_r
    beta_test(2)=beta_test(0)+dbeta_r
    beta_test(3)=beta_test(0)-j*dbeta_i
    beta_test(4)=beta_test(0)+j*dbeta_i

! evaluate the error on all points of the simplex and calculate the final error

    spatial_avg_final_error=0d0
    do itest=4,0,-1
      beta=beta_test(itest)
      plot_integrand=.FALSE.
      INCLUDE 'calc_error_Wait.F90'    
      error_test(itest)=error
      spatial_avg_final_error=spatial_avg_final_error+error/5d0
    end do
    
    final_error=error_test(0)
       
    if ( (real(final_error/initial_error).LT.1d-6).AND. &
         (real(spatial_avg_final_error/initial_error).LT.1d-1).AND. &
         (real(betaw).GE.0d0).AND. (imag(betaw).LT.0d0)  ) then
         
      write(*,'(A,I4,A,F10.6,A,F10.6,A,ES10.2,A,ES10.2)')'    mode=',mode, &
              ' betaw =',real(betaw),'+j',imag(betaw),'  point_err=',real(final_error/initial_error), &
              ' avg_err=',real(spatial_avg_final_error/initial_error)
         
! assume that the error is tending to zero and this is a real forward propagating mode with attenuation
    
      if (mode.EQ.1) then  ! this must be a new mode
        n_modes=1
        beta_found(n_modes)=betaw
        write(16,'(3ES12.4)',ADVANCE='NO')1000.0*real(gammaw),imag(gammaw),20*log10(abs(exp(-gammaw*prop_dist)))
        write(18,'(4ES12.4)')f,1000.0*real(gammaw),imag(gammaw),20*log10(abs(exp(-gammaw*prop_dist)))
       
        if (plot_beta_map) then
          write(*,*)'***************************************************'
          write(*,*)'Plotting Wait integrand for mode',n_modes
          write(*,*)'***************************************************'
          beta=betaw
          write(*,*)'beta=',cmplx(beta),' 9km attenuation=',real(20*log10(abs(exp(-gammaw*prop_dist)))),'dB'
          plot_integrand=.TRUE.
          
          INCLUDE 'calc_error_Wait.F90'    
          
        end if

! MODE 1 is sometimes the guided Zenneck wave so may choose not to plot this one          
        if (calculate_fields) then
          INCLUDE 'plot_fields.F90'    
        end if
        
!! TLM thin wire model estimation...  REMOVED
!        Ztlm=j*w*tw_L
!        Ytlm=j*w*tw_C
!        Ztot=Zw+Ztlm
!        Ytot=Yw*Ytlm/(Yw+Ytlm)
!        gammatlm=sqrt(Ztot*Ytot)
!        write(20,'(4ES12.4)')f,1000.0*real(gammatlm),imag(gammatlm),20*log10(abs(exp(-gammatlm*prop_dist)))
!         
!        write(*,*)
!        write(*,*)'mode=',n_modes
!        write(*,*)'Z=',Zw
!        write(*,*)'Y=',Yw
!        write(*,*)'R=',real(Zw)
!        write(*,*)'L=',real(Zw/(j*w))
!        write(*,*)'G=',real(Yw)
!        write(*,*)'C=',real(Yw/(j*w))
          
       else
    
! check to see if we have already found this mode
        already_found_this_mode=.FALSE.
        do mode_test=1,n_modes
          beta_diff=abs(betaw-beta_found(mode_test))
          if (beta_diff.LT.min_beta_diff) then
            already_found_this_mode=.TRUE.
          end if
        end do
    
        if(.NOT.already_found_this_mode) then
          n_modes=n_modes+1
          beta_found(n_modes)=betaw    
          write(16,'(3ES12.4)',ADVANCE='NO')1000.0*real(gammaw),imag(gammaw),20*log10(abs(exp(-gammaw*prop_dist)))
          write(18,'(4ES12.4)')f,1000.0*real(gammaw),imag(gammaw),20*log10(abs(exp(-gammaw*prop_dist)))
    
          if (plot_beta_map) then
            write(*,*)'***************************************************'
            write(*,*)'Plotting Wait integrand for mode',n_modes
            write(*,*)'***************************************************'
            beta=betaw
            write(*,*)'beta=',cmplx(beta),' 9km attenuation=',real(20*log10(abs(exp(-gammaw*prop_dist)))),'dB'
!            plot_integrand=.TRUE.
                                    
            INCLUDE 'calc_error_Wait.F90'    
            
          end if
            
          if (calculate_fields) then
            INCLUDE 'plot_fields.F90'    
          end if
          
!! TLM thin wire model estimation... REMOVED
!          Ztlm=j*w*tw_L
!          Ytlm=j*w*tw_C
!          Ztot=Zw+Ztlm
!          Ytot=Yw*Ytlm/(Yw+Ytlm)
!          gammatlm=sqrt(Ztot*Ytot)
!          write(20,'(4ES12.4)')f,1000.0*real(gammatlm),imag(gammatlm),20*log10(abs(exp(-gammatlm*prop_dist)))
                    
        end if
      
      end if  ! mode=1 or not...
      
    else
    
      write(*,'(A,I4,A,F10.6,A,F10.6,A,ES10.2,A,ES10.2)')'*** mode=',mode, &
              ' betaw =',real(betaw),'+j',imag(betaw),'  point_err=',real(final_error/initial_error), &
              ' avg_err=',real(spatial_avg_final_error/initial_error)

    end if    ! relative error is small
    
  end do ! next mode
  
  write(16,*)
  
  if (calculate_fields) then
    write(26,*)
  end if
  
8888 CONTINUE

  flush(unit=10)
  flush(unit=12)
  flush(unit=14)
  flush(unit=16)
  flush(unit=18)
  
!  flush(unit=20)

end do

close(unit=10) 
close(unit=12) 
close(unit=14) 
close(unit=16)
close(unit=18)
close(unit=20)
 
DEALLOCATE( err_array )
DEALLOCATE( beta_array )

STOP

END PROGRAM propagation_model
!
! ________________________________________________________________________________
!
!
SUBROUTINE calc_modified_Bessel_function(znu,arg,K0)

USE mod_zbes

IMPLICIT NONE

complex*16 :: znu,arg,K0

! local variables

complex*16 :: zz,H1n,H2n

integer    :: info

! Evaluate modified Bessel function

  if (imag(arg).GE.0d0) then
    zz=arg*exp((0d0,-3.1415926535d0)/2d0)
    CALL hankel2(znu,zz,H2n,info)
    if (info.NE.0) then
      write(*,*)'ERROR in hankel2:'
      write(*,*)'znu=',znu
      write(*,*)'zz =',zz
      write(*,*)'info=',info
    end if
    K0=((0d0,-3.1415926535d0)/2d0)*exp(znu*(0d0,-3.1415926535d0)/2d0)*H2n
  else
    zz=arg*exp((0d0,3.1415926535d0)/2d0)
    CALL hankel1(znu,zz,H1n,info)
    if (info.NE.0) then
      write(*,*)'ERROR in hankel2:'
      write(*,*)'znu=',znu
      write(*,*)'zz =',zz
      write(*,*)'info=',info
    end if
    K0=((0d0,3.1415926535d0)/2d0)*exp(znu*(0d0,3.1415926535d0)/2d0)*H1n
  end if
  
  RETURN

END SUBROUTINE calc_modified_Bessel_function

!
! ________________________________________________________________________________
!
!
SUBROUTINE calc_modified_Bessel_function_and_derivative(znu,arg,K0,dK0)

USE mod_zbes

IMPLICIT NONE

complex*16 :: znu,arg,K0,dK0

! local variables

complex*16 :: zz,H1n,H2n,H1nm,H2nm

complex*16 :: znum,K0nm
integer    :: info

! Evaluate modified Bessel function
 
  if (imag(arg).GE.0d0) then
  
    zz=arg*exp((0d0,-3.1415926535d0)/2d0)
    CALL hankel2(znu,zz,H2n,info)
    if (info.NE.0) then
      write(*,*)'ERROR in hankel2:'
      write(*,*)'znu=',znu
      write(*,*)'zz =',zz
      write(*,*)'info=',info
    end if
    K0=((0d0,-3.1415926535d0)/2d0)*exp(znu*(0d0,-3.1415926535d0)/2d0)*H2n
    
  else
  
    zz=arg*exp((0d0,3.1415926535d0)/2d0)
    CALL hankel1(znu,zz,H1n,info)
    if (info.NE.0) then
      write(*,*)'ERROR in hankel2:'
      write(*,*)'znu=',znu
      write(*,*)'zz =',zz
      write(*,*)'info=',info
    end if
    K0=((0d0,3.1415926535d0)/2d0)*exp(znu*(0d0,3.1415926535d0)/2d0)*H1n
    
  end if

! Evaluate modified Bessel function of order znu-1
  znum=znu-(1d0,0d0)
  
  if (imag(arg).GE.0d0) then
  
    zz=arg*exp((0d0,-3.1415926535d0)/2d0)
    CALL hankel2(znum,zz,H2nm,info)
    if (info.NE.0) then
      write(*,*)'ERROR in hankel2:'
      write(*,*)'znu=',znum
      write(*,*)'zz =',zz
      write(*,*)'info=',info
    end if
    K0nm=((0d0,-3.1415926535d0)/2d0)*exp(znum*(0d0,-3.1415926535d0)/2d0)*H2nm
    
  else
  
    zz=arg*exp((0d0,3.1415926535d0)/2d0)
    CALL hankel1(znum,zz,H1nm,info)
    if (info.NE.0) then
      write(*,*)'ERROR in hankel2:'
      write(*,*)'znu=',znum
      write(*,*)'zz =',zz
      write(*,*)'info=',info
    end if
    K0nm=((0d0,3.1415926535d0)/2d0)*exp(znum*(0d0,3.1415926535d0)/2d0)*H1nm
    
  end if

  dK0=exp(dcmplx(0d0,-3.1415926535d0))*K0nm-(znu/zz)*K0

  RETURN

END SUBROUTINE calc_modified_Bessel_function_and_derivative

