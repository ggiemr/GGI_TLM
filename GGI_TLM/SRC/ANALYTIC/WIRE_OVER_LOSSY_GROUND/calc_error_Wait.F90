  
  CALL wiat_integral(k1,k2,f,a,h,beta,ALPHA,QmjP,NmjM)

  Zw=(j*w*mu0/(2.0*pi))*(ALPHA+2d0*QmjP)
  Yw=2.0*pi*j*w*eps0/   (ALPHA+2d0*NmjM)
  
  gammaw=sqrt(Zw*Yw)
  betaw=gammaw/j
  new_beta=gammaw/j
  
  error=abs(beta-new_beta)
