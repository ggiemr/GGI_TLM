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
! SUBROUTINE S_to_Z
!
! NAME
!    S_to_Z
!
! DESCRIPTION
!     S parameter to impedance parameter transformation
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 19/09/2013 CJS as S_to_Z_symmetric
!     generalised 4/12/2013 becomming S_to_Z
!
SUBROUTINE S_to_Z


USE post_process
USE file_information

IMPLICIT NONE

! local variables
  
  real*8	:: f
  complex*16	:: S11,S12,S21,S22
  complex*16	:: Z11,Z12,Z21,Z22
  real*8	:: Z0,Z1,Z2
   
  complex*16 	:: S(2,2),IPS(2,2),IMS(2,2),IMSinv(2,2),Z(2,2),DET
 
  integer	:: n_frequencies
  integer	:: frequency_loop
  integer	:: function_number
  
  character	:: ch
  logical	:: symmetric
  integer	:: opformat
  
  integer	:: nports
  
! START

!  write(*,*)'S_to_Z transformation'

  write(*,*)'How many ports? (1/2)'
  read(*,*)nports
  write(record_user_inputs_unit,*)nports,'  ! number of ports'
  
  if (nports.EQ.1) then
  
    n_functions_of_time=0
    n_functions_of_frequency=2
  
    CALL Allocate_post_data()
  
    write(*,*)'Stage 1: Read S11 data'
    CALL read_one_port_frequency_domain_S_parameter_data(1)
    
    write(*,*)'Enter the reference wave impedance, Z0'
    read(*,*)Z0  
    write(record_user_inputs_unit,*)Z0,' Reference wave impedance, Z0'
       
    write(*,*)'Enter the output format (1: nomal frequency domain function or 2: fdmat format) '
    read(*,*)opformat
    write(record_user_inputs_unit,'(I10,A)')opformat,' output format (1: nomal frequency domain function or 2: fdmat format)'

! Allocate memory for result

    n_frequencies=function_of_frequency(1)%n_frequencies

    function_number=2
    
    function_of_frequency(function_number)%n_frequencies=n_frequencies  
    
    ALLOCATE ( function_of_frequency(function_number)%frequency(1:n_frequencies) )
    function_of_frequency(function_number)%frequency(1:n_frequencies)=  	      &
    		  function_of_frequency(1)%frequency(1:n_frequencies)	  
    ALLOCATE ( function_of_frequency(function_number)%value(1:n_frequencies) )
    ALLOCATE ( function_of_frequency(function_number)%magnitude(1:n_frequencies) )
    ALLOCATE ( function_of_frequency(function_number)%phase(1:n_frequencies) )
    ALLOCATE ( function_of_frequency(function_number)%dB(1:n_frequencies) )
  
    do frequency_loop=1,n_frequencies
    
      S11=function_of_frequency(1)%value(frequency_loop)
      Z(1,1)= Z0*(1d0+S11)/(1d0-S11)
        
      function_of_frequency(2)%value(frequency_loop)=Z(1,1)
      function_of_frequency(2)%magnitude(frequency_loop)=	    &
    		    abs(function_of_frequency(2)%value(frequency_loop))
      function_of_frequency(2)%phase(frequency_loop)=   &
    		      atan2( imag(function_of_frequency(2)%value(frequency_loop)), &
    			     dble(function_of_frequency(2)%value(frequency_loop))	)
      function_of_frequency(2)%dB(frequency_loop)=      &
    		      20d0*log10(function_of_frequency(2)%magnitude(frequency_loop))

    end do ! next frequency value
    
    if (opformat.eq.1) then
  
      write(*,*)'Write Z11 data to file'
      CALL write_Frequency_Domain_Data(2)
  
    else
  
      write(*,*)'Write Z11 data to file'
      CALL write_Frequency_Domain_Data_simple_format(2)
  
    end if
  
  else
! two ports
  
    n_functions_of_time=0
    n_functions_of_frequency=8
  
    CALL Allocate_post_data()
  
    write(*,*)'Stage 1: Read S11, S21 data'
    CALL read_frequency_domain_S_parameter_data(1,3)
  
    write(*,*)'Is the structure symmetrical ? (y/n)'
    read(*,'(A)')ch
    write(record_user_inputs_unit,'(A)')ch
    if ( (ch.eq.'Y').OR.(ch.eq.'y') ) then
      symmetric=.TRUE.
    else if ( (ch.eq.'N').OR.(ch.eq.'n') ) then
      symmetric=.FALSE.
    else
      write(*,*)'Response should be y/n'
      STOP
    end if  
  
    if (symmetric) then
! Allocate memory for S12, S22 and copy over from the already read S11 and S21

! S12=S21
      n_frequencies=function_of_frequency(3)%n_frequencies
      function_of_frequency(2)%n_frequencies=n_frequencies  
      ALLOCATE ( function_of_frequency(2)%frequency(1:n_frequencies) )
      ALLOCATE ( function_of_frequency(2)%value(1:n_frequencies) )
      ALLOCATE ( function_of_frequency(2)%magnitude(1:n_frequencies) )
      ALLOCATE ( function_of_frequency(2)%phase(1:n_frequencies) )
      ALLOCATE ( function_of_frequency(2)%dB(1:n_frequencies) )
    
      function_of_frequency(2)%frequency(1:n_frequencies)=function_of_frequency(3)%frequency(1:n_frequencies)     
      function_of_frequency(2)%value(1:n_frequencies)    =function_of_frequency(3)%value(1:n_frequencies)
      function_of_frequency(2)%magnitude(1:n_frequencies)=function_of_frequency(3)%magnitude(1:n_frequencies)
      function_of_frequency(2)%phase(1:n_frequencies)    =function_of_frequency(3)%phase(1:n_frequencies)    
      function_of_frequency(2)%dB(1:n_frequencies)       =function_of_frequency(3)%dB(1:n_frequencies)	   

! S22=S11
      n_frequencies=function_of_frequency(1)%n_frequencies
      function_of_frequency(4)%n_frequencies=n_frequencies  
      ALLOCATE ( function_of_frequency(4)%frequency(1:n_frequencies) )
      ALLOCATE ( function_of_frequency(4)%value(1:n_frequencies) )
      ALLOCATE ( function_of_frequency(4)%magnitude(1:n_frequencies) )
      ALLOCATE ( function_of_frequency(4)%phase(1:n_frequencies) )
      ALLOCATE ( function_of_frequency(4)%dB(1:n_frequencies) )
    
      function_of_frequency(4)%frequency(1:n_frequencies)=function_of_frequency(1)%frequency(1:n_frequencies)     
      function_of_frequency(4)%value(1:n_frequencies)    =function_of_frequency(1)%value(1:n_frequencies)
      function_of_frequency(4)%magnitude(1:n_frequencies)=function_of_frequency(1)%magnitude(1:n_frequencies)
      function_of_frequency(4)%phase(1:n_frequencies)    =function_of_frequency(1)%phase(1:n_frequencies)    
      function_of_frequency(4)%dB(1:n_frequencies)       =function_of_frequency(1)%dB(1:n_frequencies)	   
  
    else
      write(*,*)'Stage 2: Read S12, S22 data'
      CALL read_frequency_domain_S_parameter_data(2,4)
    end if
  
    write(*,*)'Enter the reference wave impedance, Z0'
    read(*,*)Z0

! could be generalised but assume reference impedance is the same on both sides and constant with frequency  
    Z1=Z0
    Z2=Z0
  
    write(record_user_inputs_unit,*)Z0,' Reference wave impedance, Z0'
       
    write(*,*)'Enter the output format (1: nomal frequency domain function or 2: fdmat format) '
    read(*,*)opformat
    write(record_user_inputs_unit,'(I10,A)')opformat,' output format (1: nomal frequency domain function or 2: fdmat format)'

! Allocate memory for result

    n_frequencies=function_of_frequency(1)%n_frequencies

    do function_number=5,8
      function_of_frequency(function_number)%n_frequencies=n_frequencies  
      ALLOCATE ( function_of_frequency(function_number)%frequency(1:n_frequencies) )
      function_of_frequency(function_number)%frequency(1:n_frequencies)=		&
                    function_of_frequency(1)%frequency(1:n_frequencies)     
      ALLOCATE ( function_of_frequency(function_number)%value(1:n_frequencies) )
      ALLOCATE ( function_of_frequency(function_number)%magnitude(1:n_frequencies) )
      ALLOCATE ( function_of_frequency(function_number)%phase(1:n_frequencies) )
      ALLOCATE ( function_of_frequency(function_number)%dB(1:n_frequencies) )
    end do
  
    do frequency_loop=1,n_frequencies
  
      S11=function_of_frequency(1)%value(frequency_loop)
      S12=function_of_frequency(2)%value(frequency_loop)
      S21=function_of_frequency(3)%value(frequency_loop)
      S22=function_of_frequency(4)%value(frequency_loop)
    
! Impose reciprocity i.e. make S12=S21

      S12=(S12+S21)/2d0
      S21=S12
    
      IMS(1,1)= ((1D0,0D0)-S11)
      IMS(1,2)= (	        -S12)
      IMS(2,1)=-(	        -S21)
      IMS(2,2)=-((1D0,0D0)-S22)
   
      IPS(1,1)=(1D0,0D0)+S11
      IPS(1,2)=	      +S12
      IPS(2,1)=	      +S21
      IPS(2,2)=(1D0,0D0)+S22
   
      DET=IMS(1,1)*IMS(2,2)-IMS(1,2)*IMS(2,1)
   
      IMSinv(1,1)= IMS(2,2)/DET
      IMSinv(1,2)=-IMS(1,2)/DET
      IMSinv(2,1)=-IMS(2,1)/DET
      IMSinv(2,2)= IMS(1,1)/DET
#
! Note the sign convention for the impedance parameters...   
      Z(1,1)= (IPS(1,1)*IMSinv(1,1)+IPS(1,2)*IMSinv(2,1))*Z1
      Z(1,2)=-(IPS(1,1)*IMSinv(1,2)+IPS(1,2)*IMSinv(2,2))*Z2
      Z(2,1)= (IPS(2,1)*IMSinv(1,1)+IPS(2,2)*IMSinv(2,1))*Z1
      Z(2,2)=-(IPS(2,1)*IMSinv(1,2)+IPS(2,2)*IMSinv(2,2))*Z2
       
      function_of_frequency(5)%value(frequency_loop)=Z(1,1)
      function_of_frequency(6)%value(frequency_loop)=Z(1,2)
      function_of_frequency(7)%value(frequency_loop)=Z(2,1)
      function_of_frequency(8)%value(frequency_loop)=Z(2,2)   
    
      do function_number=5,8
        function_of_frequency(function_number)%magnitude(frequency_loop)=	&
                      abs(function_of_frequency(function_number)%value(frequency_loop))
        function_of_frequency(function_number)%phase(frequency_loop)=	&
                        atan2( imag(function_of_frequency(function_number)%value(frequency_loop)), &
                               dble(function_of_frequency(function_number)%value(frequency_loop))   )
        function_of_frequency(function_number)%dB(frequency_loop)=	&
                        20d0*log10(function_of_frequency(function_number)%magnitude(frequency_loop))
      end do

    end do ! next frequency value
    
    if (opformat.eq.1) then
  
      function_number=5
      write(*,*)'Write Z11 data to file'
      CALL write_Frequency_Domain_Data(function_number)
  
      function_number=6
      write(*,*)'Write Z12 data to file'
      CALL write_Frequency_Domain_Data(function_number)
 
      function_number=7
      write(*,*)'Write Z21 data to file'
      CALL write_Frequency_Domain_Data(function_number)
  
      function_number=8
      write(*,*)'Write Z22 data to file'
      CALL write_Frequency_Domain_Data(function_number)

    else
  
      function_number=5
      write(*,*)'Write Z11 data to file'
      CALL write_Frequency_Domain_Data_simple_format(function_number)
  
      function_number=6
      write(*,*)'Write Z12 data to file'
      CALL write_Frequency_Domain_Data_simple_format(function_number)
 
      function_number=7
      write(*,*)'Write Z21 data to file'
      CALL write_Frequency_Domain_Data_simple_format(function_number)
  
      function_number=8
      write(*,*)'Write Z22 data to file'
      CALL write_Frequency_Domain_Data_simple_format(function_number)
  
    end if

  end if ! one or two ports

  CALL Deallocate_post_data()

  RETURN
  
END SUBROUTINE S_to_Z

