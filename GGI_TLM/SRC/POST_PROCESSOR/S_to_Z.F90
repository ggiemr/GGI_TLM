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
! SUBROUTINE S_to_Z_symmetric
!
! NAME
!    S_to_Z_symmetric
!
! DESCRIPTION
!     S parameter to impedance parameter transformation for symmetric structures
!     i.e. S11=S22, Z11=Z22
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 19/09/2013 CJS
!
!
SUBROUTINE S_to_Z_symmetric


USE post_process
USE file_information

IMPLICIT NONE

! local variables
  
  real*8	:: f
  complex*16	:: S11,S21,Z11,Z21
  real*8	:: Z0
   
  complex 	:: S(2,2),IPS(2,2),IMS(2,2),IMSinv(2,2),Z(2,2),DET
 
  integer	:: n_frequencies
  integer	:: frequency_loop
  integer	:: function_number
  
! START

!  write(*,*)'S_to_Z transformation for symmetric structures'
  
  n_functions_of_time=0
  n_functions_of_frequency=4
  
  CALL Allocate_post_data()
  
  CALL read_frequency_domain_S_parameter_data(1,2)
  
  write(*,*)'Enter the reference wave impedance, Z0'
  read(*,*)Z0
  
  write(record_user_inputs_unit,*)Z0,' Reference wave impedance, Z0'

! Allocate memory for result

  n_frequencies=function_of_frequency(1)%n_frequencies

  function_number=3  
  function_of_frequency(function_number)%n_frequencies=n_frequencies  
  ALLOCATE ( function_of_frequency(function_number)%frequency(1:n_frequencies) )
  function_of_frequency(function_number)%frequency(1:n_frequencies)=		&
                function_of_frequency(1)%frequency(1:n_frequencies)     
  ALLOCATE ( function_of_frequency(function_number)%value(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%magnitude(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%phase(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%dB(1:n_frequencies) )

  function_number=4 
  function_of_frequency(function_number)%n_frequencies=n_frequencies  
  ALLOCATE ( function_of_frequency(function_number)%frequency(1:n_frequencies) )
  function_of_frequency(function_number)%frequency(1:n_frequencies)=		&
                function_of_frequency(1)%frequency(1:n_frequencies)     
  ALLOCATE ( function_of_frequency(function_number)%value(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%magnitude(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%phase(1:n_frequencies) )
  ALLOCATE ( function_of_frequency(function_number)%dB(1:n_frequencies) )
  
  do frequency_loop=1,n_frequencies
  
    S11=function_of_frequency(1)%value(frequency_loop)
    S21=function_of_frequency(2)%value(frequency_loop)
    
    IMS(1,1)= ((1.0,0.0)-S11)/Z0
    IMS(1,2)= (	       -S21)/Z0
    IMS(2,1)=-(	       -S21)/Z0
    IMS(2,2)=-((1.0,0.0)-S11)/Z0
   
    IPS(1,1)=(1.0,0.0)+S11
    IPS(1,2)=	      +S21
    IPS(2,1)=	      +S21
    IPS(2,2)=(1.0,0.0)+S11
   
    DET=IMS(1,1)*IMS(2,2)-IMS(1,2)*IMS(2,1)
   
    IMSinv(1,1)= IMS(2,2)/DET
    IMSinv(1,2)=-IMS(1,2)/DET
    IMSinv(2,1)=-IMS(2,1)/DET
    IMSinv(2,2)= IMS(1,1)/DET
   
    Z(1,1)=IMSinv(1,1)*IPS(1,1)+IMSinv(1,2)*IPS(2,1)
    Z(1,2)=IMSinv(1,1)*IPS(1,2)+IMSinv(1,2)*IPS(2,2)
    Z(2,1)=IMSinv(2,1)*IPS(1,1)+IMSinv(2,2)*IPS(2,1)
    Z(2,2)=IMSinv(2,1)*IPS(1,2)+IMSinv(2,2)*IPS(2,2)
       
    function_number=3    
    function_of_frequency(function_number)%value(frequency_loop)=Z(1,1)
    function_of_frequency(function_number)%magnitude(frequency_loop)=	&
                    abs(function_of_frequency(function_number)%value(frequency_loop))
    function_of_frequency(function_number)%phase(frequency_loop)=	&
                    atan2( imag(function_of_frequency(function_number)%value(frequency_loop)), &
                           dble(function_of_frequency(function_number)%value(frequency_loop))   )
    function_of_frequency(function_number)%dB(frequency_loop)=	&
                    20d0*log10(function_of_frequency(function_number)%magnitude(frequency_loop))
    
    function_number=4 
    function_of_frequency(function_number)%value(frequency_loop)=Z(2,1)
    function_of_frequency(function_number)%magnitude(frequency_loop)=	&
                    abs(function_of_frequency(function_number)%value(frequency_loop))
    function_of_frequency(function_number)%phase(frequency_loop)=	&
                    atan2( imag(function_of_frequency(function_number)%value(frequency_loop)), &
                           dble(function_of_frequency(function_number)%value(frequency_loop))   )
    function_of_frequency(function_number)%dB(frequency_loop)=	&
                    20d0*log10(function_of_frequency(function_number)%magnitude(frequency_loop))

  end do ! next frequency value
  
  function_number=3
  write(*,*)'Write Z11 data to file'
  CALL write_Frequency_Domain_Data(function_number)
  
  function_number=4
  write(*,*)'Write Z12 data to file'
  CALL write_Frequency_Domain_Data(function_number)
 
  function_number=4
  write(*,*)'Write Z21 data to file'
  CALL write_Frequency_Domain_Data(function_number)
  
  function_number=3
  write(*,*)'Write Z22 data to file'
  CALL write_Frequency_Domain_Data(function_number)

  CALL Deallocate_post_data()

  RETURN
  

  
END SUBROUTINE S_to_Z_symmetric
