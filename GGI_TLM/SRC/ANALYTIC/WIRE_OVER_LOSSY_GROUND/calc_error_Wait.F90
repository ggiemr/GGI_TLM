  if (Olsen_eqn) then

! use the form from the Olsen paper  
    CALL Wait_integral(k1,k2,f,a,h,beta,ALPHA,QmjP,NmjM,plot_integrand, &
                       show_plots,Jc,nlambda,lambda,dl_integral)

    Zw=(j*w*mu0/(2.0*pi))*(ALPHA+2d0*QmjP)
    Yw=(2.0*pi*j*w*eps0 )/(ALPHA+2d0*NmjM)
  
    gammaw=sqrt(Zw*Yw)
    betaw=gammaw/j
    new_beta=gammaw/j
  
    error=abs(beta-new_beta)

  else
! use equation 17 in the Wait paper  

    CALL Wait_integral2(k1,k2,f,a,h,beta,Kterm,Integral_term,ALPHA,QmjP,NmjM, &
                        plot_integrand,show_plots,Jc,nlambda,lambda,dl_integral)

    Zw=(j*w*mu0/(2.0*pi))*(ALPHA+2d0*QmjP)
    Yw=(2.0*pi*j*w*eps0 )/(ALPHA+2d0*NmjM)
  
    gammaw=sqrt(Zw*Yw)
    betaw=gammaw/j
    new_beta=gammaw/j
   
    error=abs(Kterm+Integral_term)
  
  end if
