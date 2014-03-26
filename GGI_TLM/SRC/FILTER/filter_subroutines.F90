!
!    GGI_TLM Time domain electromagnetic field solver based on the TLM method
!    Copyright (C) 2013  Chris Smartt
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.   
!
!  SUBROUTINE read_Sfilter(s1,ip_unit)
!  SUBROUTINE write_Sfilter(s1,op_unit)
!  SUBROUTINE write_S_PZ_filter(s1,op_unit)
!  SUBROUTINE write_S_PR_filter(s1,op_unit)
!  SUBROUTINE read_Zfilter(Z1,ip_unit)
!  SUBROUTINE write_Zfilter(z1,op_unit)
!  SUBROUTINE output_Sfilter_frequency_response(s1,fmin,fmax,fstep,op_unit)
!  SUBROUTINE output_Zfilter_frequency_response(z1,fmin,fmax,fstep,op_unit)
!  SUBROUTINE output_Zfilter_frequency_response_fast_slow(Zf,Zs,fmin,fmax,fstep,op_unit)
!  SUBROUTINE evaluate_Zfilter(Z1,ft,G)
!  SUBROUTINE timeshift_Zfilter(ft)
!  SUBROUTINE bilinear_s_to_z( filter_s_in , dt , filter_z )
!  SUBROUTINE bilinear_z_to_s( filter_z  , filter_s )
!  SUBROUTINE reciprocal_Sfilter(a,res) 
!  SUBROUTINE Z_fast_slow_decomposition( Z_in, f_fast  , Z_slow )
!  SUBROUTINE copy_Zfilter_response(ft1,ft2)
!  SUBROUTINE subtract_1_Sfilter(a,res) 
!  SUBROUTINE calculate_magnetic_susceptibility_impedance_Sfilter(a,K,res) 
!  SUBROUTINE calculate_electric_susceptibility_impedance_Sfilter(a,K,res) 
!  SUBROUTINE bicubic_s_to_z
!  SUBROUTINE deallocate_Zfilter_data(Zfilter_data)
!  SUBROUTINE deallocate_Sfilter
!
! NAME
!     read_Sfilter
!
! DESCRIPTION
!       read Laplace domain filter coefficients from file
!
! SEE ALSO
!
!
! HISTORY
!
!     started 02/03/09 CJS
!
  SUBROUTINE read_Sfilter(s1,ip_unit)

USE filter_types
  
  type(Sfilter)	:: s1
  integer	:: ip_unit
  
  integer i

!START

! read filter coefficients from the given file unit

  read(ip_unit,*)s1%wnorm
  read(ip_unit,*)s1%a%order
  if (allocated (s1%a%coeff) ) deallocate(s1%a%coeff)
  allocate( s1%a%coeff(0:s1%a%order) )
  read(ip_unit,*)(s1%a%coeff(i),i=0,s1%a%order)  
  read(ip_unit,*)s1%b%order
  if (allocated (s1%b%coeff) ) deallocate(s1%b%coeff)
  allocate( s1%b%coeff(0:s1%b%order) )
  read(ip_unit,*)(s1%b%coeff(i),i=0,s1%b%order)  
  
  RETURN
  END SUBROUTINE read_Sfilter
!
! NAME
!     write_Sfilter
!
! DESCRIPTION
!       write Laplace domain filter coefficients to screen or file 
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
  SUBROUTINE write_Sfilter(s1,op_unit)

USE filter_types
  
  type(Sfilter)	:: s1
  integer	:: op_unit 
  
  integer i

!START

  if (op_unit.eq.0) then
! write to screen
    write(*,*)'Laplace domain filter'
    write(*,*)'wnorm=',s1%wnorm
    write(*,*)'a order=',s1%a%order
    write(*,*)(s1%a%coeff(i),i=0,s1%a%order)  
    write(*,*)'b order=',s1%b%order
    write(*,*)(s1%b%coeff(i),i=0,s1%b%order)  
    write(*,*)' '
  else
    write(op_unit,*)s1%wnorm,'  # w normalisation constant'
    write(op_unit,*)s1%a%order,'  # a order, a coefficients follow below:'
    write(op_unit,*)(s1%a%coeff(i),i=0,s1%a%order)  
    write(op_unit,*)s1%b%order,'  # b order, b coefficients follow below:'
    write(op_unit,*)(s1%b%coeff(i),i=0,s1%b%order)  
  end if
  
  RETURN
  END SUBROUTINE write_Sfilter
!
! NAME
!     write_S_PZ_filter
!
! DESCRIPTION
!       write Laplace domain filter in pole/ zero format to screen or file 
!
! SEE ALSO
!
!
! HISTORY
!
!     started 4/01/13 CJS
!
  SUBROUTINE write_S_PZ_filter(s1,op_unit)

USE filter_types
  
  type(Sfilter_PZ)	:: s1
  integer		:: op_unit 
  
  integer i

!START

  if (op_unit.eq.0) then
! write to screen
    write(*,*)'Laplace domain filter, Pole Zero format'
    write(*,*)'wnorm=',s1%wnorm
    write(*,*)'order=',s1%order
    write(*,*)'G    =',s1%G
    write(*,*)'n_real_poles        =',s1%n_real_poles
    write(*,*)'n_complex_poles     =',s1%n_complex_poles
    write(*,*)'n_complex_pole_pairs=',s1%n_complex_pole_pairs
    do i=1,s1%order
      write(*,*)'Pole',i,' ',s1%poles(i)
    end do
    write(*,*)'n_real_zeros        =',s1%n_real_zeros
    write(*,*)'n_complex_zeros     =',s1%n_complex_zeros
    write(*,*)'n_complex_zero_pairs=',s1%n_complex_zero_pairs
    do i=1,s1%order
      write(*,*)'zero',i,' ',s1%zeros(i)
    end do
  else
    write(op_unit,*)s1%wnorm,' wnorm'
    write(op_unit,*)s1%order,' order'
    write(op_unit,*)s1%G,' G '
    write(op_unit,*)s1%n_real_poles,' n_real_poles'
    write(op_unit,*)s1%n_complex_poles,' n_complex_poles'
    write(op_unit,*)s1%n_complex_pole_pairs,' n_complex_pole_pairs'
    do i=1,s1%order
      write(op_unit,*)s1%poles(i),' Pole',i
    end do
    write(op_unit,*)s1%n_real_zeros,' n_real_zeros'
    write(op_unit,*)s1%n_complex_zeros,' n_complex_zeros'
    write(op_unit,*)s1%n_complex_zero_pairs,' n_complex_zero_pairs'
    do i=1,s1%order
      write(op_unit,*)s1%zeros(i),' zero',i
    end do
  end if
  
  RETURN
  END SUBROUTINE write_S_PZ_filter
!
! NAME
!     write_S_PR_filter
!
! DESCRIPTION
!       write Laplace domain filter in pole/ residue format to screen or file 
!
! SEE ALSO
!
!
! HISTORY
!
!     started 4/01/13 CJS
!
  SUBROUTINE write_S_PR_filter(s1,op_unit)

USE filter_types
  
  type(Sfilter_PR)	:: s1
  integer		:: op_unit 
  
  integer i

!START
  REAL*8  			:: wnorm
  INTEGER			:: order
  REAL*8			:: C
  INTEGER			:: n_complex_poles
  INTEGER			:: n_complex_pole_pairs
  INTEGER			:: n_real_poles
  LOGICAL,allocatable		:: complex_pole(:)
  COMPLEX*16,allocatable	:: poles(:)
  COMPLEX*16,allocatable	:: residues(:)

  if (op_unit.eq.0) then
! write to screen
    write(*,*)'Laplace domain filter, Pole Zero format'
    write(*,*)'wnorm=',s1%wnorm
    write(*,*)'order=',s1%order
    write(*,*)'C    =',s1%C
    write(*,*)'n_complex_poles     =',s1%n_complex_poles
    write(*,*)'n_complex_pole_pairs=',s1%n_complex_pole_pairs
    write(*,*)'n_real_poles        =',s1%n_real_poles
    do i=1,s1%order
      write(*,*)'Pole',i,' ',s1%poles(i)
    end do
    do i=1,s1%order
      write(*,*)'residue',i,' ',s1%residues(i)
    end do
  else
    write(op_unit,*)s1%wnorm,' wnorm'
    write(op_unit,*)s1%order,' order'
    write(op_unit,*)s1%C,' C '
    write(op_unit,*)s1%n_real_poles,' n_real_poles'
    write(op_unit,*)s1%n_complex_poles,' n_complex_poles'
    write(op_unit,*)s1%n_complex_pole_pairs,' n_complex_pole_pairs'
    do i=1,s1%order
      write(op_unit,*)s1%poles(i),' Pole',i
    end do
    do i=1,s1%order
      write(op_unit,*)s1%residues(i),' residue',i
    end do
  end if
  
  RETURN
  END SUBROUTINE write_S_PR_filter
!
! NAME
!     read_Zfilter
!
! DESCRIPTION
!       read Z domain filter coefficients from file 
!
! SEE ALSO
!
!
! HISTORY
!
!     started 14/12/2012 CJS
!
  SUBROUTINE read_Zfilter(Z1,ip_unit)

USE filter_types
  
  type(Zfilter)	:: Z1
  integer	:: ip_unit
  
  integer i

!START

! read filter coefficients from the given file unit

  read(ip_unit,*)Z1%wnorm
  read(ip_unit,*)Z1%T
  read(ip_unit,*)Z1%a%order
  if (allocated (Z1%a%coeff) ) deallocate(Z1%a%coeff)
  allocate( Z1%a%coeff(0:Z1%a%order) )
  read(ip_unit,*)(Z1%a%coeff(i),i=0,Z1%a%order)  
  read(ip_unit,*)Z1%b%order
  if (allocated (Z1%b%coeff) ) deallocate(Z1%b%coeff)
  allocate( Z1%b%coeff(0:Z1%b%order) )
  read(ip_unit,*)(Z1%b%coeff(i),i=0,Z1%b%order)  
  
  RETURN
  END SUBROUTINE read_Zfilter
!
! NAME
!     write_Zfilter
!
! DESCRIPTION
!       write Z transform filter coefficients to screen or file 
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
  SUBROUTINE write_Zfilter(z1,op_unit)

USE filter_types
  
  type(Zfilter) :: z1
  integer 	:: op_unit
  
  integer i

!START

  if (op_unit.eq.0) then
! write to screen
    write(*,*)'Z transform filter'
    write(*,*)'wnorm   =',z1%wnorm
    write(*,*)'Timestep=',z1%T
    write(*,*)'a order   =',z1%a%order
    write(*,*)(z1%a%coeff(i),i=0,z1%a%order)  
    write(*,*)'b order   =',z1%b%order
    write(*,*)(z1%b%coeff(i),i=0,z1%b%order)  
    write(*,*)' '
  else
    write(op_unit,*)z1%wnorm,' wnorm '
    write(op_unit,*)z1%T,' Timestep'
    write(op_unit,*)z1%a%order,' a order'
    write(op_unit,*)(z1%a%coeff(i),i=0,z1%a%order)
    write(op_unit,*)z1%b%order,' b order'
    write(op_unit,*)(z1%b%coeff(i),i=0,z1%b%order)
  end if
  
  RETURN
  END SUBROUTINE write_Zfilter
  
!
! NAME
!      output_Sfilter_frequency_response
!     
! DESCRIPTION
!       Evaluate the frequency resaponse of a Laplace domain filter and
! 	write the response to file
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
  SUBROUTINE output_Sfilter_frequency_response(s1,fmin,fmax,fstep,op_unit)

USE filter_types
USE filter_functions
  
  type(Sfilter)	:: s1
  real*8	:: fmin,fmax,fstep
  integer	:: op_unit 

! local variables  
  real*8 f,w
  complex*16 jw,num,den,response
  integer n
  integer n_frequencies,frequency_loop

!START

!  do f=fmin,fmax,fstep

  n_frequencies=int( (fmax-fmin)/fstep )+1

  do frequency_loop=1,n_frequencies

    f=fmin+(frequency_loop-1)*fstep
  
    response=evaluate_Sfilter_frequency_response(s1,f)
        
    write(op_unit,1000)f,real(response),imag(response)
1000 format(3E16.7)

  end do
  
  RETURN
  END SUBROUTINE output_Sfilter_frequency_response
  
!
! NAME
!      output_Zfilter_frequency_response
!     
! DESCRIPTION
!       Evaluate the frequency resaponse of a Z domain filter and
! 	write the response to file
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/01/09 CJS
!
  SUBROUTINE output_Zfilter_frequency_response(z1,fmin,fmax,fstep,op_unit)

USE filter_types
USE filter_functions
  
  type(Zfilter)	:: z1
  real*8	:: fmin,fmax,fstep
  integer	:: op_unit 

! local variables  
  real*8 f,w
  complex*16 jw,num,den,response
  integer n
  integer n_frequencies,frequency_loop

!START

!  do f=fmin,fmax,fstep

  n_frequencies=int( (fmax-fmin)/fstep )+1

  do frequency_loop=1,n_frequencies

    f=fmin+(frequency_loop-1)*fstep
  
    response=evaluate_Zfilter_frequency_response(z1,f)

    write(op_unit,1000)f,real(response),imag(response)
1000 format(3E16.7)

  end do
  
  RETURN
  END SUBROUTINE output_Zfilter_frequency_response
  
  
!
! NAME
!      output_Zfilter_frequency_response_fast_slow
!     
! DESCRIPTION
!       Evaluate the frequency resaponse of a Z domain filter which
!       has been decomposed into the fast/ slow format and
! 	write the response to file
!
! SEE ALSO
!
!
! HISTORY
!
!     started 21/02/11 CJS
!
  SUBROUTINE output_Zfilter_frequency_response_fast_slow(Zf,Zs,fmin,fmax,fstep,op_unit)

USE filter_types
USE filter_functions
  
  real*8	:: Zf
  type(Zfilter)	:: Zs
  real*8	:: fmin,fmax,fstep
  integer	:: op_unit 

! local variables  
  real*8 f,w
  complex*16 jw,num,den,response
  integer n
  integer n_frequencies,frequency_loop

!START

!  do f=fmin,fmax,fstep

  n_frequencies=int( (fmax-fmin)/fstep )+1

  do frequency_loop=1,n_frequencies

    f=fmin+(frequency_loop-1)*fstep
  
    response=evaluate_Zfilter_frequency_response(Zs,f)
    
    response=response+Zf

    write(op_unit,1000)f,real(response),imag(response)
1000 format(3E16.7)

  end do
  
  RETURN
  END SUBROUTINE output_Zfilter_frequency_response_fast_slow
  
!
! NAME
!      evaluate_Zfilter
!      
! DESCRIPTION
!       Evaluate the filter response for the current timestep with the
!       forcing function value =G, filter coeffs =Z1, filter store, ft
! 	
!
! SEE ALSO
!
!
! HISTORY
!
!     started 28/01/09 CJS
!
  SUBROUTINE evaluate_Zfilter(Z1,ft,G)

USE filter_types
USE filter_functions
  
  type(Zfilter)		  :: z1
  type(Zfilter_response)  :: ft
  real*8		  :: G

! local variables  
  integer i
  real*8 W0,F0

!START

  W0=G
  do i=1,ft%order
    W0=W0-z1%b%coeff(i)*ft%w(i)
  end do
  W0=W0/z1%b%coeff(0)
  ft%w(0)=W0

  F0=0d0
  do i=0,ft%order
    F0=F0+z1%a%coeff(i)*ft%w(i)
  end do
  ft%f=F0

  
  RETURN
  END SUBROUTINE evaluate_Zfilter
  
!
! NAME
!      timeshift_Zfilter(ft)
!      
! DESCRIPTION
!       time shift filter backstore variables for the next timestep
! 	
!
! SEE ALSO
!
!
! HISTORY
!
!     started 29/01/09 CJS
!
  SUBROUTINE timeshift_Zfilter(ft)

USE filter_types
USE filter_functions
  
  type(Zfilter_response)  :: ft

! local variables  
  integer i

!START

! shift W
  do i=ft%order,1,-1
    ft%w(i)=ft%w(i-1)
  end do
  
  RETURN
  END SUBROUTINE timeshift_Zfilter
!
! NAME
!       bilinear_s_to_z
!
! DESCRIPTION
!       Bilinear transformation from s plane to z plane 
!       i.e. calculate the digital filter coefficients from the 
!       s plane representaion of the filter function and the
!       calculation timestep
!
!       The filter function is found as a polynomial in z^(-1) 
!       This is the revised material format   
!
!
! SEE ALSO
!     polynomial_operators.f90
!
! HISTORY
!
!     started 17/05/05 CJS
!
! COMMENTS
!     makes use of operator overloading for polynomial types     
!

SUBROUTINE bilinear_s_to_z( filter_s_in , dt , filter_z )
  
  ! Modules used
  
  USE filter_types
  USE polynomial_types
  USE polynomial_operators
  USE polynomial_functions
  
  IMPLICIT NONE

  !Variables passed as arguments to subroutine

  type(Sfilter)       :: filter_s_in
  type(Zfilter)       :: filter_z
  real*8              :: dt
  
  ! Local variables
  type(Sfilter)       :: filter_s
        
  integer order,max_order
  integer aorder,border
  
  real*8 gain
  
  type(polynomial)     :: a,b1,b2,c
  
  type(polynomial)     :: snum,sden
  type(polynomial)     :: znum,zden
  type(polynomial)     :: znum2,zden2
  type(polynomial)     :: ztemp
  
  complex*16,allocatable:: poles(:)
  complex*16,allocatable:: rpoles(:),cpoles(:)
  complex*16,allocatable:: roots(:)
  complex*16,allocatable:: rroots(:),croots(:)
  integer nreal,ncomplex
  
  real*8 rec2,modc,b0
  
  real*8 dt2
  integer i,j
  integer loop
    
  ! START

!  print*,'Called Bilinear s to z'

! copy filter, reducing the order if highest order terms are zero  

  filter_s%a%order=0
  do i=0,filter_s_in%a%order
    if (filter_s_in%a%coeff(i).NE.0d0) then
      filter_s%a%order=filter_s_in%a%order
    end if
  end do

  filter_s%b%order=0
  do i=0,filter_s_in%b%order
    if (filter_s_in%b%coeff(i).NE.0d0) then
      filter_s%b%order=filter_s_in%b%order
    end if
  end do

  filter_s%wnorm=filter_s_in%wnorm
  allocate (filter_s%a%coeff(0:filter_s%a%order))
  allocate (filter_s%b%coeff(0:filter_s%b%order))
  
  do i=0,filter_s%a%order
    filter_s%a%coeff(i)=filter_s_in%a%coeff(i)
  end do
  
  do i=0,filter_s%b%order
    filter_s%b%coeff(i)=filter_s_in%b%coeff(i)
  end do    
  
!  write(*,*)'Bilinear_s_to_z'
!  write(*,*)'filter_in'
!  CALL write_Sfilter(filter_s_in,0)
!  write(*,*)'filter'
!  CALL write_Sfilter(filter_s,0)

  max_order=filter_s%a%order
  if (filter_s%b%order.gt.max_order) max_order=filter_s%b%order
  
  aorder=filter_s%a%order
  border=filter_s%b%order
  
  dt2=2d0/(dt*filter_s%wnorm)

  allocate (poles(1:border))
  allocate (rpoles(1:border))
  allocate (cpoles(1:border))
  
  allocate (roots(1:aorder))
  allocate (rroots(1:aorder))
  allocate (croots(1:aorder))

!find poles and roots of transfer function and sort into 
!real and complex

!copy numerator and denominator of transfer function into polynomial types 
  
  snum=filter_s%a
  sden=filter_s%b

!find roots of the polynomials ! note: roots and poles are already allocated
 call findroots(snum,roots,aorder)
 call findroots(sden,poles,border)

!calculate the gain term
if (sden%coeff(border).ne.0d0) then
  gain=snum%coeff(aorder)/sden%coeff(border)
else
  write(*,*)'Problem calculating gain, sden%coeff(max_order)=0'
  write(*,*)'max_order=',max_order
  write(*,*)'snum=',snum%coeff(0:aorder)
  write(*,*)'sden=',sden%coeff(0:border)
  stop
end if

! loop over numerator then denominator and calculate 
! {1+z^(-1)}a(z) for numerator terms and 
! {1+z^(-1)}b(z) for denominator terms

do loop=1,2
  
  if(loop.eq.1) then
!do numerator terms 
    call rootsort(aorder,roots,rroots,		            &
                    croots,nreal,ncomplex,aorder)
    order=aorder
  else 
!do denominator terms 
    call rootsort(border,poles,rpoles,		            &
                    cpoles,nreal,ncomplex,border)
    order=border
  end if

! initialise work polynomials
  a=1.0D0
  b2%order=2
  ALLOCATE (b2%coeff(0:b2%order))
  b1%order=1
  ALLOCATE (b1%coeff(0:b1%order))
  c=0.0D0

!transform complex roots

  do i=1,ncomplex
  
    if(loop.eq.1) then
      modc=dble(croots(i*2-1)*croots(i*2))
      rec2=dble(croots(i*2-1)+croots(i*2))
    else
      modc=dble(cpoles(i*2-1)*cpoles(i*2))
      rec2=dble(cpoles(i*2-1)+cpoles(i*2)) 
    end if
    
! revised format
    b2%coeff(0)=dt2*dt2-rec2*dt2+modc
    b2%coeff(1)=-2d0*dt2*dt2+2d0*modc
    b2%coeff(2)=dt2*dt2+rec2*dt2+modc
        
    c=a*b2
    a=c
    
  end do

!transform real roots

  do i=1,nreal
      
! revised format...    
    
    if(loop.eq.1) then
      b1%coeff(0)=dt2-rroots(i)
      b1%coeff(1)=-dt2-rroots(i)
    else
      b1%coeff(0)=dt2-rpoles(i)
      b1%coeff(1)=-dt2-rpoles(i)   
    end if
    
    c=a*b1
    a=c
    
  end do
  
  DEALLOCATE (b2%coeff)
  DEALLOCATE (b1%coeff)
  
  if(loop.eq.1) then
!put numerator terms, scaled by the gain into znum 

    znum=allocate_polynomial(aorder)
    do i=0,aorder
      znum%coeff(i)=a%coeff(i)*gain
    end do
  else   
!put denominator terms, into zden 
    zden=allocate_polynomial(border)
    do i=0,border
      zden%coeff(i)=a%coeff(i)
    end do
  end if

! return to do denominator terms
end do

! multiply top and bottom by (1+z^(-1))

ztemp=allocate_polynomial(1)

ztemp%coeff(0)=1d0
ztemp%coeff(1)=1d0

if (aorder.gt.border) then
  order=aorder-border
  znum2=znum
  zden2=zden
  do i=1,order      
    a=zden2*ztemp
    zden2=a
  end do
else if (border.gt.aorder) then
  order=border-aorder
  znum2=znum
  zden2=zden
  do i=1,order      
    a=znum2*ztemp
    znum2=a
  end do
else ! numerator and denominator orders are the same so complete cancellation.
  znum2=znum
  zden2=zden
end if

!scale so that b0=1.0

  b0=zden2%coeff(0)
  if (b0.ne.0d0) then
    do i=0,max_order
      znum2%coeff(i)=znum2%coeff(i)/b0
    end do
    do i=0,max_order
      zden2%coeff(i)=zden2%coeff(i)/b0
    end do
  end if
  
!copy numerator and denominator of transfer function into z_filter type 
  filter_z%wnorm     =filter_s%wnorm 
  filter_z%T         =dt*filter_s%wnorm
  filter_z%a     =znum2
  filter_z%b     =zden2

deallocate (poles)
deallocate (rpoles)
deallocate (cpoles)
deallocate (roots)
deallocate (rroots)
deallocate (croots)

RETURN

END
!
! NAME
!       bilinear_z_to_s
!
! DESCRIPTION
!       Bilinear transformation from z plane to s plane 
!       i.e. calculate the s plane representation from the
!       digital filter coefficients 
!
!
! SEE ALSO
!     polynomial_operators.f90
!
! HISTORY
!
!     started 05/02/08 CJS
!
! COMMENTS
!     makes use of operator overloading for polynomial types     
!

SUBROUTINE bilinear_z_to_s( filter_z  , filter_s )
  
  ! Modules used
  
  USE filter_types
  USE polynomial_types
  USE polynomial_operators
  USE polynomial_functions
  
  IMPLICIT NONE

  !Variables passed as arguments to subroutine

  type(Zfilter)       :: filter_z
  type(Sfilter)       :: filter_s
  
  ! Local variables
        
  integer order,max_order
  
  real*8 gain
  
  type(polynomial)     :: a,b1,b2,c
  
  type(polynomial)     :: snum,sden
  type(polynomial)     :: znum,zden
  
  complex*16,allocatable:: poles(:),roots(:)
  complex*16,allocatable:: rroots(:),croots(:)
  integer nreal,ncomplex
  
  real*8 rec2,modc,b0
  
  real*8 dt2,dt
  integer i,j
  integer loop
  
  ! START

  print*,'Called Bilinear z to s'

! copy filter, reducing the order if highest order terms are zero  

  max_order=filter_z%a%order
  if (filter_z%b%order.gt.max_order) max_order=filter_z%b%order
  
  dt=filter_z%T
  dt2=dt/2d0

  allocate (poles(1:max_order))
  allocate (roots(1:max_order))
  allocate (rroots(1:max_order))
  allocate (croots(1:max_order))

!find poles and roots of transfer function and sort into 
!real and complex

!copy numerator and denominator of transfer function into polynomial types 
  
  znum=filter_z%a
  zden=filter_z%b

print*,'Order=',max_order

print*,'A'
print*,znum%coeff(0:max_order)

print*,'B'
print*,zden%coeff(0:max_order)

!find roots of the polynomials ! note: roots and poles are already allocated
 call findroots(znum,roots,max_order)
 call findroots(zden,poles,max_order)

print*,'roots'
do i=1,max_order
  print*,roots(i)
end do

print*,'poles'
do i=1,max_order
  print*,poles(i)
end do

!calculate the gain term
gain=znum%coeff(max_order)/zden%coeff(max_order)

print*,'gain=',gain

! loop over numerator then denominator and calculate 
! {1+z^(-1)}a(z) for numerator terms and 
! {1+z^(-1)}b(z) for denominator terms

do loop=1,2
  
  if(loop.eq.1) then
!do numerator terms 
    call rootsort(max_order,roots,rroots,		            &
                    croots,nreal,ncomplex,max_order)
  else 
!do denominator terms 
    call rootsort(max_order,poles,rroots,		            &
                    croots,nreal,ncomplex,max_order)
  end if

! initialise work polynomials
  a=1.0D0
  b2%order=2
  ALLOCATE (b2%coeff(0:b2%order))
  b1%order=1
  ALLOCATE (b1%coeff(0:b1%order))
  c=0.0D0

!transform complex roots

  do i=1,ncomplex
  
    modc=dble(croots(i*2-1)*croots(i*2))
    rec2=dble(croots(i*2-1)+croots(i*2))

! revised format

    b2%coeff(0)=1d0-rec2+modc
    b2%coeff(1)=2d0*dt2*(-1d0+modc)
    b2%coeff(2)=dt2*dt2*(1d0+rec2+modc)
        
    c=a*b2
    a=c
    
  end do

!transform real roots

  do i=1,nreal
      
! revised format...    
    b1%coeff(0)=1d0-rroots(i)
    b1%coeff(1)=-dt2*(1d0+rroots(i))
    
    c=a*b1
    a=c
    
  end do
  
  DEALLOCATE (b2%coeff)
  DEALLOCATE (b1%coeff)

  print*,'A=',a%coeff(0:max_order)

  if(loop.eq.1) then
!put numerator terms, scaled by the gain into znum 

    snum=allocate_polynomial(max_order)
    do i=0,max_order
      snum%coeff(i)=a%coeff(i)*gain
    end do
  else   
!put denominator terms, into zden 
    sden=allocate_polynomial(max_order)
    do i=0,max_order
      sden%coeff(i)=a%coeff(i)
    end do
  end if

! return to do denominator terms
end do

!scale so that b0=1.0

  b0=sden%coeff(0)
  if (b0.ne.0d0) then
    do i=0,max_order
      snum%coeff(i)=snum%coeff(i)/b0
    end do
    do i=0,max_order
      sden%coeff(i)=sden%coeff(i)/b0
    end do
  end if
  
!  print*,'Snum=',snum
!  print*,'Sden=',sden
  
!copy numerator and denominator of transfer function into z_filter type 
  filter_s%wnorm     =filter_z%wnorm 
  filter_s%a     =snum
  filter_s%b     =sden

deallocate (poles)
deallocate (roots)
deallocate (rroots)
deallocate (croots)

RETURN

END
!
! NAME
!     reciprocal_Sfilter
!
! DESCRIPTION
!    input an S plane filter, A(s)/B(s) and return 
!    the reciprocal filter i.e. B(s)/A(s)
!
! HISTORY
!
!     started 22/01/09 CJS
!

SUBROUTINE reciprocal_Sfilter(a,res) 

USE filter_types

IMPLICIT NONE

! argument types
  type(Sfilter) :: a
  
! Result type
  type(Sfilter) :: ans
  type(Sfilter) :: res

! function definition

  ans%wnorm=a%wnorm
  ans%a%order=a%b%order
  ans%b%order=a%a%order
  allocate (ans%a%coeff(0:ans%a%order))
  allocate (ans%b%coeff(0:ans%b%order))
  
  ans%a%coeff(0:ans%a%order)=a%b%coeff(0:ans%a%order)
  ans%b%coeff(0:ans%b%order)=a%a%coeff(0:ans%b%order)
    
  if (allocated(res%a%coeff) ) deallocate(res%a%coeff)
  if (allocated(res%b%coeff) ) deallocate(res%b%coeff)
  res=ans
  
END SUBROUTINE reciprocal_Sfilter

!
! NAME
!       Z_fast_slow_decomposition
!
! DESCRIPTION
!       decompose a Z domain filter into an instantaneous (fast) and
!       delayed (slow) system
!
!
! SEE ALSO
!     
!
! HISTORY
!
!     started 07/04/09 CJS
!
! COMMENTS
!     makes use of operator overloading for polynomial types     
!

SUBROUTINE Z_fast_slow_decomposition( Z_in, f_fast  , Z_slow )
  
  ! Modules used
  
  USE filter_types 
  USE filter_functions
  USE filter_operators
  USE polynomial_types
  USE polynomial_operators
  USE polynomial_functions
  
  IMPLICIT NONE

  !Variables passed as arguments to subroutine

  type(Zfilter)       :: Z_in
  real*8 	      :: f_fast
  type(Zfilter)       :: Z_slow
  
  ! Local variables
  
  type(Zfilter)       :: Z_out

  integer aorder,border,max_order
  
  integer i
    
! START

!  print*,'Z_fast_slow_decomposition'

  aorder=Z_in%a%order
  border=Z_in%b%order

  if (aorder.ne.border) then
    write(*,*)'Error in Z_fast_slow_decomposition'
    write(*,*)'aorder.ne.border'
    write(*,*)aorder,border
    stop
  end if
  
  f_fast=Z_in%a%coeff(0)/Z_in%b%coeff(0)
    
  Z_out=Z_in
    
  Z_out%a%coeff(0)=0d0
  do i=1,border
    Z_out%a%coeff(i)=Z_out%a%coeff(i)-Z_in%b%coeff(i)*f_fast
  end do
    
  Z_slow=Z_out
  
  RETURN

END
!
! NAME
!    copy_Zfilter_response  
!      
! DESCRIPTION
!       
! 	
!
! SEE ALSO
!
!
! HISTORY
!
!     started 23/11/10 CJS
!
  SUBROUTINE copy_Zfilter_response(ft1,ft2)

USE filter_types
USE filter_functions
  
  type(Zfilter_response)  :: ft1
  type(Zfilter_response)  :: ft2
  
! local variables  
  integer i,order

!START

! check order of filter functions
  if (ft1%order.ne.ft2%order) then
    write(*,*)'Error in copy_Zfilter_response'
    write(*,*)'filter order discrepancy'
    write(*,*)'ft1 order=',ft1%order
    write(*,*)'ft2 order=',ft2%order
    stop
  end if
  
  order=ft1%order
  ft2%f=ft1%f
  ft2%w(0:order)=ft1%w(0:order)
  
  RETURN
  END SUBROUTINE copy_Zfilter_response
  
!
! NAME
!     subtract_1_Sfilter
!
! DESCRIPTION
!    input an S plane filter, A(s)/B(s) and return 
!    the filter function A(s)/B(s)-1=(A(s)-B(s))/B(s)
!
! HISTORY
!
!     started 22/01/09 CJS
!

SUBROUTINE subtract_1_Sfilter(a,res) 

USE filter_types

IMPLICIT NONE

! argument types
  type(Sfilter) :: a
  
! Result type
  type(Sfilter) :: ans
  type(Sfilter) :: res
  
  integer i

! function definition

  ans%wnorm=a%wnorm
  ans%a%order=a%a%order
  ans%b%order=a%b%order
  allocate (ans%a%coeff(0:ans%a%order))
  allocate (ans%b%coeff(0:ans%b%order))
  
  ans%a%coeff(0:ans%a%order)=a%a%coeff(0:ans%a%order)
  ans%b%coeff(0:ans%b%order)=a%b%coeff(0:ans%b%order)
  
  do i=0,a%b%order
    ans%a%coeff(i)=ans%a%coeff(i)-ans%b%coeff(i)
  end do
  
  if (allocated(res%a%coeff) ) deallocate(res%a%coeff)
  if (allocated(res%b%coeff) ) deallocate(res%b%coeff)
  res=ans
  
END SUBROUTINE subtract_1_Sfilter

  
!
! NAME
!     calculate_magnetic_susceptibility_impedance_Sfilter
!
! DESCRIPTION
!    input an S plane filter, mur(s)=A(s)/B(s) and return the stub impedance filter Z(s)=s*Lstub
!    where Lstub=K*(mur(s)-1)
!
!   Note: we must take into account the filter frequency normalisation 
!
!
! HISTORY
!
!     started 06/09/12 CJS
!

SUBROUTINE calculate_magnetic_susceptibility_impedance_Sfilter(a,K,res) 

USE filter_types

IMPLICIT NONE

! argument types
  type(Sfilter) :: a
  real*8	:: K
  
! Result type
  type(Sfilter) :: ans
  type(Sfilter) :: res
  
  integer i

! function definition

  ans%wnorm=a%wnorm
  ans%a%order=a%a%order+1
  ans%b%order=a%b%order
  allocate (ans%a%coeff(0:ans%a%order))
  allocate (ans%b%coeff(0:ans%b%order))

! numerator terms  
  ans%a%coeff(0)=0d0
  do i=0,a%a%order
    ans%a%coeff(i+1)=K*a%wnorm*a%a%coeff(i)
  end do
  do i=0,a%b%order
    ans%a%coeff(i+1)=ans%a%coeff(i+1)-K*a%wnorm*a%b%coeff(i)
  end do
  
  ans%b%coeff(0:ans%b%order)=a%b%coeff(0:ans%b%order)
  
  if (allocated(res%a%coeff) ) deallocate(res%a%coeff)
  if (allocated(res%b%coeff) ) deallocate(res%b%coeff)
  res=ans
  
END SUBROUTINE calculate_magnetic_susceptibility_impedance_Sfilter
!
! NAME
!     calculate_electric_susceptibility_impedance_Sfilter
!
! DESCRIPTION
!    input an S plane filter, epsr(s)=A(s)/B(s) and return the stub impedance filter Z(s)=1/(s*Cstub(s))
!    where Cstub(s)=K*(epsr-1)
!
!   Note: we must take into account the filter frequency normalisation 
!
!
! HISTORY
!
!     started 06/09/12 CJS
!

SUBROUTINE calculate_electric_susceptibility_impedance_Sfilter(a,K,res) 

USE filter_types

IMPLICIT NONE

! argument types
  type(Sfilter) :: a
  real*8	:: K
  
! Result type
  type(Sfilter) :: ans
  type(Sfilter) :: res
  
  integer i

! function definition

  ans%wnorm=a%wnorm
  ans%b%order=a%a%order+1
  ans%a%order=a%b%order
  allocate (ans%a%coeff(0:ans%a%order))
  allocate (ans%b%coeff(0:ans%b%order))

! denominator terms  
  ans%b%coeff(0)=0d0
  do i=0,a%a%order
    ans%b%coeff(i+1)=K*a%wnorm*a%a%coeff(i)
  end do
  do i=0,a%b%order
    ans%b%coeff(i+1)=ans%b%coeff(i+1)-K*a%wnorm*a%b%coeff(i)
  end do
  
  ans%a%coeff(0:a%b%order)=a%b%coeff(0:a%b%order)
  
  if (allocated(res%a%coeff) ) deallocate(res%a%coeff)
  if (allocated(res%b%coeff) ) deallocate(res%b%coeff)
  res=ans
  
END SUBROUTINE calculate_electric_susceptibility_impedance_Sfilter


!
! NAME
!       bicubic_s_to_z
!
! DESCRIPTION
!       Bicubic transformation from s plane to z plane 
!       i.e. calculate the digital filter coefficients from the 
!       s plane representaion of the filter function and the
!       calculation timestep
!
!       The filter function is found as a polynomial in z^(-3) 
!       This is the revised material format   
!
!
! SEE ALSO
!     polynomial_operators.f90
!
! HISTORY
!
!     started 17/04/12 CJS
!
! COMMENTS
!     makes use of operator overloading for polynomial types     
!

SUBROUTINE bicubic_s_to_z( filter_s_in , dt , filter_z )
  
  ! Modules used
  
  USE filter_types
  USE polynomial_types
  USE polynomial_operators
  USE polynomial_functions
  
  IMPLICIT NONE

  !Variables passed as arguments to subroutine

  type(Sfilter)       :: filter_s_in
  type(Zfilter)       :: filter_z
  real*8              :: dt
  
  ! Local variables
        
  type(Sfilter)       :: filter_s
  integer order,max_order
  integer aorder,border
  
  real*8 gain
  
  type(polynomial)     :: a,b1,b2,c
  
  type(polynomial)     :: snum,sden
  type(polynomial)     :: znum,zden
  type(polynomial)     :: znum2,zden2
  type(polynomial)     :: ztemp
  
  complex*16,allocatable:: poles(:)
  complex*16,allocatable:: rpoles(:),cpoles(:)
  complex*16,allocatable:: roots(:)
  complex*16,allocatable:: rroots(:),croots(:)
  integer nreal,ncomplex
  
  real*8 rec2,modc,b0
  
  real*8 dt23
  integer i,j
  integer loop
    
  ! START

! copy filter, reducing the order if highest order terms are zero  
  filter_s%a%order=0
  do i=0,filter_s_in%a%order
    if (filter_s_in%a%coeff(i).NE.0d0) then
      filter_s%a%order=filter_s_in%a%order
    end if
  end do

  filter_s%b%order=0
  do i=0,filter_s_in%b%order
    if (filter_s_in%b%coeff(i).NE.0d0) then
      filter_s%b%order=filter_s_in%b%order
    end if
  end do

  filter_s%wnorm=filter_s_in%wnorm
  allocate (filter_s%a%coeff(0:filter_s%a%order))
  allocate (filter_s%b%coeff(0:filter_s%b%order))
  
  do i=0,filter_s%a%order
    filter_s%a%coeff(i)=filter_s_in%a%coeff(i)
  end do
  
  do i=0,filter_s%b%order
    filter_s%b%coeff(i)=filter_s_in%b%coeff(i)
  end do    
  
  max_order=filter_s%a%order
  if (filter_s%b%order.gt.max_order) max_order=filter_s%b%order
  
  aorder=filter_s%a%order
  border=filter_s%b%order
  
  dt23=2d0/(3d0*dt*filter_s%wnorm)

  allocate (poles(1:border))
  allocate (rpoles(1:border))
  allocate (cpoles(1:border))
  
  allocate (roots(1:aorder))
  allocate (rroots(1:aorder))
  allocate (croots(1:aorder))

!find poles and roots of transfer function and sort into 
!real and complex

!copy numerator and denominator of transfer function into polynomial types 
  
  snum=filter_s%a
  sden=filter_s%b

!find roots of the polynomials ! note: roots and poles are already allocated
 call findroots(snum,roots,aorder)
 call findroots(sden,poles,border)

!calculate the gain term
if (sden%coeff(border).ne.0d0) then
  gain=snum%coeff(aorder)/sden%coeff(border)
else
  write(*,*)'Problem calculating gain, sden%coeff(max_order)=0'
  write(*,*)'max_order=',max_order
  write(*,*)'snum=',snum%coeff(0:aorder)
  write(*,*)'sden=',sden%coeff(0:border)
  stop
end if

! loop over numerator then denominator and calculate 
! {1+z^(-3)}a(z) for numerator terms and 
! {1+z^(-3)}b(z) for denominator terms

do loop=1,2
  
  if(loop.eq.1) then
!do numerator terms 
    call rootsort(aorder,roots,rroots,		            &
                    croots,nreal,ncomplex,aorder)
    order=aorder
  else 
!do denominator terms 
    call rootsort(border,poles,rpoles,		            &
                    cpoles,nreal,ncomplex,border)
    order=border
  end if

! initialise work polynomials; b2 for complex roots, b1 for real roots
  a=1.0D0
  
  b2%order=6
  ALLOCATE (b2%coeff(0:b2%order))
  
  b1%order=3
  ALLOCATE (b1%coeff(0:b1%order))
  
  c=0.0D0

!transform complex roots

  do i=1,ncomplex
  
    if(loop.eq.1) then
      modc=dble(croots(i*2-1)*croots(i*2))
      rec2=dble(croots(i*2-1)+croots(i*2))
    else
      modc=dble(cpoles(i*2-1)*cpoles(i*2))
      rec2=dble(cpoles(i*2-1)+cpoles(i*2)) 
    end if
    
! revised format
    b2%coeff(0)=dt23*dt23-rec2*dt23+modc
    b2%coeff(1)=0d0
    b2%coeff(2)=0d0
    b2%coeff(3)=-2d0*dt23*dt23+2d0*modc
    b2%coeff(4)=0d0
    b2%coeff(5)=0d0
    b2%coeff(6)=dt23*dt23+rec2*dt23+modc
        
    c=a*b2
    a=c
    
  end do

!transform real roots

  do i=1,nreal
      
! revised format...    
    
    if(loop.eq.1) then
      b1%coeff(0)=dt23-rroots(i)
      b1%coeff(1)=0d0
      b1%coeff(2)=0d0
      b1%coeff(3)=-dt23-rroots(i)
    else
      b1%coeff(0)=dt23-rpoles(i)
      b1%coeff(1)=0d0
      b1%coeff(2)=0d0
      b1%coeff(3)=-dt23-rpoles(i)   
    end if
    
    c=a*b1
    a=c
    
  end do
  
  DEALLOCATE (b2%coeff)
  DEALLOCATE (b1%coeff)
  
  if(loop.eq.1) then
!put numerator terms, scaled by the gain into znum 

    znum=allocate_polynomial(3*aorder)
    do i=0,3*aorder
      znum%coeff(i)=a%coeff(i)*gain
    end do
  else   
!put denominator terms, into zden 
    zden=allocate_polynomial(3*border)
    do i=0,3*border
      zden%coeff(i)=a%coeff(i)
    end do
  end if

! return to do denominator terms
end do

! multiply top and bottom by (1+z^(-3))

ztemp=allocate_polynomial(3)

ztemp%coeff(0)=1d0
ztemp%coeff(1)=0d0
ztemp%coeff(2)=0d0
ztemp%coeff(3)=1d0

! need to check this... is the order determinaton correct?

if (aorder.gt.border) then
  order=aorder-border
  znum2=znum
  zden2=zden
  do i=1,order      
    a=zden2*ztemp
    zden2=a
  end do
else if (border.gt.aorder) then
  order=border-aorder
  znum2=znum
  zden2=zden
  do i=1,order      
    a=znum2*ztemp
    znum2=a
  end do
else ! numerator and denominator orders are the same so complete cancellation.
  znum2=znum
  zden2=zden
end if

!scale so that b0=1.0

  b0=zden2%coeff(0)
  if (b0.ne.0d0) then
    do i=0,znum2%order
      znum2%coeff(i)=znum2%coeff(i)/b0
    end do
    do i=0,zden2%order
      zden2%coeff(i)=zden2%coeff(i)/b0
    end do
  end if
  
!copy numerator and denominator of transfer function into z_filter type 
  filter_z%wnorm     =filter_s%wnorm 
  filter_z%T         =dt*filter_s%wnorm
  filter_z%a     =znum2
  filter_z%b     =zden2

deallocate (poles)
deallocate (rpoles)
deallocate (cpoles)
deallocate (roots)
deallocate (rroots)
deallocate (croots)

RETURN

END SUBROUTINE bicubic_s_to_z 
!
! NAME
!      deallocate_Zfilter(ft)
!      
! DESCRIPTION
!       deallocate Zfilter type
! 	
!
! SEE ALSO
!
!
! HISTORY
!
!     started 26/09/12 CJS
!
  SUBROUTINE deallocate_Zfilter(Zfilter_in)

USE filter_types
USE filter_functions
  
  type(Zfilter)  :: Zfilter_in

! local variables  
  integer i

!START

  if (allocated(  Zfilter_in%a%coeff )) then
    DEALLOCATE( Zfilter_in%a%coeff )
  end if 
  if (allocated(  Zfilter_in%b%coeff )) then
    DEALLOCATE( Zfilter_in%b%coeff )
  end if 
  
  RETURN
  END SUBROUTINE deallocate_Zfilter
!
! NAME
!      deallocate_Sfilter(ft)
!      
! DESCRIPTION
!       deallocate Sfilter type
! 	
!
! SEE ALSO
!
!
! HISTORY
!
!     started 26/09/12 CJS
!
  SUBROUTINE deallocate_Sfilter(Sfilter_in)

USE filter_types
USE filter_functions
  
  type(Sfilter)  :: Sfilter_in

! local variables  
  integer i

!START

  if (allocated(  Sfilter_in%a%coeff )) then
    DEALLOCATE( Sfilter_in%a%coeff )
  end if 
  if (allocated(  Sfilter_in%b%coeff )) then
    DEALLOCATE( Sfilter_in%b%coeff )
  end if 
  
  RETURN
  END SUBROUTINE deallocate_Sfilter
!
! NAME
!      deallocate_Zfilter_data
!      
! DESCRIPTION
!       deallocate Zfilter_data 
! 	
!
! SEE ALSO
!
!
! HISTORY
!
!     started 26/09/12 CJS
!
  SUBROUTINE deallocate_Zfilter_data(Zfilter_data)

USE filter_types
USE filter_functions
  
  type(Zfilter_response)  :: Zfilter_data

! local variables  

!START

  if (allocated(  Zfilter_data%w )) then
    DEALLOCATE( Zfilter_data%w )
  end if 
   
  RETURN
  END SUBROUTINE deallocate_Zfilter_data
