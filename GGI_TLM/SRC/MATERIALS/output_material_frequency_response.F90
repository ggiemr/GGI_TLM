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
!
! NAME
!      output_material_frequency_response
!     
! DESCRIPTION
!       Evaluate the frequency response of a Laplace domain material filter and
! 	write the response to file
!
! SEE ALSO
!
!
! HISTORY
!
!     started 29/11/10 CJS
!
  subroutine output_material_frequency_response(s1,sigma,fmin,fmax,fstep,op_unit,filename)

USE filter_types
USE filter_functions
USE constants
  
  type(Sfilter)	:: s1
  real*8	:: sigma
  real*8	:: fmin,fmax,fstep
  integer	:: op_unit 
  character*256  :: filename

! local variables  
  real*8 f,w
  complex*16 jw,response
  integer n
  integer n_frequencies,frequency_loop

!START
  open(unit=op_unit,FILE=filename)
  
  write(*,*)'Name:',filename
  write(*,*)'A order=',s1%a%order
  write(*,*)'A coefficients'
  do n=0,s1%a%order
    write(*,8000)n,s1%a%coeff(n)
  end do
  write(*,*)'B order=',s1%b%order
  write(*,*)'B coefficients'
  do n=0,s1%b%order
    write(*,8000)n,s1%b%coeff(n)
  end do
8000 format(I5,E14.6)

! loop over frequency
  n_frequencies=int( (fmax-fmin)/fstep )+1

  do frequency_loop=1,n_frequencies

    f=fmin+(frequency_loop-1)*fstep

    jw=j*2d0*pi*f  
    response=evaluate_Sfilter_frequency_response(s1,f)
    
    if (abs(jw).gt.1e-10) then  
      response=response+sigma/(jw)
    end if
        
    write(op_unit,1000)f,dble(response),dimag(response)
1000 format(3E16.7)

  end do
  
  close(unit=op_unit)
  
  return
  
  end subroutine output_material_frequency_response
