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
! SUBROUTINE twoport
!
! NAME
!    twoport
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 18/01/2013 CJS
!     4/9/2014 CJS   Add an absorptance column to output data i.e. A=sqrt(S11**2-S21**2)
!
!
SUBROUTINE twoport


USE post_process
USE file_information
USE constants

IMPLICIT NONE

! local variables

  integer	:: function_number
  integer	:: n_frequencies
  integer	:: frequency_loop
  
  real*8	:: frequency,w
  
  real*8	:: f1,f2,f3,f4,f5,f6,f7,f8
  
  character(len=256)	:: filename
  
  real*8	:: d_port
  real*8	:: d_interface
  
  integer	:: Excitation_port
       
  real*8	:: d1,l1,d2,l2
  
  complex*16	:: FT(8)
  
  complex*16	:: beta1,beta2,Z1,Z2,El1,Er1,El2,Er2
  
  complex*16	:: S11,S12,S21,S22
  
  real*8	:: Pr,Pt,Pa
    
  integer,parameter :: E11=1
  integer,parameter :: E12=2
  integer,parameter :: H11=3
  integer,parameter :: H12=4
  integer,parameter :: E21=5
  integer,parameter :: E22=6
  integer,parameter :: H21=7
  integer,parameter :: H22=8

! START

  write(*,*)'Two port analysis'

  write(*,*)'Port configuration:'
  write(*,*)''
  write(*,*)'  op11   PORT_1  op12            | scatterer |         op22   PORT_2  op21  '
  write(*,*)''
       
! Allocate frequency domain data   
  
  n_functions_of_time=0
  n_functions_of_frequency=8
  
  CALL Allocate_post_data()
  
! read input files     
  write(*,*)'Read the file for E field at port 1, op11'
  function_number=E11
  CALL read_frequency_domain_data(function_number)
  
  write(*,*)'Read the file for E field at port 1, op12'
  function_number=E12
  CALL read_frequency_domain_data(function_number)
  
  write(*,*)'Read the file for H field at port 1, op11'
  function_number=H11
  CALL read_frequency_domain_data(function_number)
  
  write(*,*)'Read the file for H field at port 1, op12'
  function_number=H12
  CALL read_frequency_domain_data(function_number)
  
  write(*,*)'Read the file for E field at port 2, op21'
  function_number=E21
  CALL read_frequency_domain_data(function_number)
  
  write(*,*)'Read the file for E field at port 2, op22'
  function_number=E22
  CALL read_frequency_domain_data(function_number)
  
  write(*,*)'Read the file for H field at port 2, op21'
  function_number=H21
  CALL read_frequency_domain_data(function_number)
  
  write(*,*)'Read the file for H field at port 2, op22'
  function_number=H22
  CALL read_frequency_domain_data(function_number)
  
! check that the frequencies match...

  if ( (function_of_frequency(1)%n_frequencies.NE.function_of_frequency(2)%n_frequencies).OR.	&
       (function_of_frequency(1)%n_frequencies.NE.function_of_frequency(3)%n_frequencies).OR.	&
       (function_of_frequency(1)%n_frequencies.NE.function_of_frequency(4)%n_frequencies).OR.	&
       (function_of_frequency(1)%n_frequencies.NE.function_of_frequency(5)%n_frequencies).OR.	&
       (function_of_frequency(1)%n_frequencies.NE.function_of_frequency(6)%n_frequencies).OR.	&
       (function_of_frequency(1)%n_frequencies.NE.function_of_frequency(7)%n_frequencies).OR.	&
       (function_of_frequency(1)%n_frequencies.NE.function_of_frequency(8)%n_frequencies) ) then
  
    write(*,*)'Frequency mismatch between functions'
    write(*,*)'n_frequencies, f1=',function_of_frequency(1)%n_frequencies
    write(*,*)'n_frequencies, f2=',function_of_frequency(2)%n_frequencies
    write(*,*)'n_frequencies, f3=',function_of_frequency(3)%n_frequencies
    write(*,*)'n_frequencies, f4=',function_of_frequency(4)%n_frequencies
    write(*,*)'n_frequencies, f5=',function_of_frequency(5)%n_frequencies
    write(*,*)'n_frequencies, f6=',function_of_frequency(6)%n_frequencies
    write(*,*)'n_frequencies, f7=',function_of_frequency(7)%n_frequencies
    write(*,*)'n_frequencies, f8=',function_of_frequency(8)%n_frequencies
    STOP
    
  end if
  
  do frequency_loop=1,function_of_frequency(1)%n_frequencies  
    
    f1=function_of_frequency(1)%frequency(frequency_loop)
    f2=function_of_frequency(2)%frequency(frequency_loop)
    f3=function_of_frequency(3)%frequency(frequency_loop)
    f4=function_of_frequency(4)%frequency(frequency_loop)
    f5=function_of_frequency(5)%frequency(frequency_loop)
    f6=function_of_frequency(6)%frequency(frequency_loop)
    f7=function_of_frequency(7)%frequency(frequency_loop)
    f8=function_of_frequency(8)%frequency(frequency_loop)
  
    if ( (f1.NE.f2).OR.(f1.NE.f3).OR.(f1.NE.f4).OR.(f1.NE.f5).OR.(f1.NE.f6).OR.(f1.NE.f7).OR.(f1.NE.f8) ) then	 
      write(*,*)'Frequency mismatch in functions f1, f2, f3, f4, f5, f6, f7 and f8'
      write(*,*)'frequency number',frequency_loop
      write(*,*)'frequency, f1=',f1
      write(*,*)'frequency, f2=',f2
      write(*,*)'frequency, f3=',f3
      write(*,*)'frequency, f4=',f4
      write(*,*)'frequency, f5=',f5
      write(*,*)'frequency, f6=',f6
      write(*,*)'frequency, f7=',f7
      write(*,*)'frequency, f8=',f8
      STOP     
    end if  
    
  end do ! next frequency to check
  
! Read other solution data  
        
  write(*,*)'Enter port separation at output plane 1  (m)'
  read(*,*)d_port
  write(record_user_inputs_unit,*)d_port,' port separation at output plane 1 (m)'
  d1=d_port/2d0
  
  write(*,*)'Enter port midpoint to interface distamce at output plane 1 (m)'
  read(*,*)d_interface
  write(record_user_inputs_unit,*)d_interface,' port midpoint to interface distamce at output plane 1  (m)'
  l1=d_interface
        
  write(*,*)'Enter port separation at output plane 2  (m)'
  read(*,*)d_port
  write(record_user_inputs_unit,*)d_port,' port separation at output plane 2 (m)'
  d2=d_port/2d0
  
  write(*,*)'Enter port midpoint to interface distamce at output plane 2 (m)'
  read(*,*)d_interface
  write(record_user_inputs_unit,*)d_interface,' port midpoint to interface distamce at output plane 2 (m)'
  l2=d_interface
       
  write(*,*)'Enter the port number (1 or 2) from which the problem is excited'
  read(*,*)Excitation_port
  write(record_user_inputs_unit,'(I10,A)')Excitation_port,' Excitation_port '
 
! File for results

  write(*,*)'Enter the frequency domain data output filename'
  read(*,*)filename
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  OPEN(unit=local_file_unit,file=filename)
  
!  OPEN(unit=local_file_unit2,file='temp.dat')
  
  if (excitation_port.eq.1) then
    write(local_file_unit,*)	&
    '# f,dble(S11),dimag(S11),abs(S11),20.0*log10(abs(S11)),dble(S21),dimag(S21),abs(S21),20.0*log10(abs(S21)), Pr, Pt, Pa'
  else if (excitation_port.eq.2) then  
    write(local_file_unit,*)	&
    '# f,dble(S12),dimag(S12),abs(S12),20.0*log10(abs(S12)),dble(S22),dimag(S22),abs(S22),20.0*log10(abs(S22)), Pr, Pt, Pa'
  end if
  
  n_frequencies=function_of_frequency(1)%n_frequencies
  
  do frequency_loop=1,n_frequencies
  
    frequency=function_of_frequency(1)%frequency(frequency_loop)
    w=2D0*pi*frequency
 
    FT(E11)=function_of_frequency(E11)%value(frequency_loop)
    FT(E12)=function_of_frequency(E12)%value(frequency_loop)
    FT(H11)=function_of_frequency(H11)%value(frequency_loop)
    FT(H12)=function_of_frequency(H12)%value(frequency_loop)
    FT(E21)=function_of_frequency(E21)%value(frequency_loop)
    FT(E22)=function_of_frequency(E22)%value(frequency_loop)
    FT(H21)=function_of_frequency(H21)%value(frequency_loop)
    FT(H22)=function_of_frequency(H22)%value(frequency_loop)
    
    CALL oneport_calc(FT(E11),FT(E12),FT(H11),FT(H12),d1,beta1,Z1,Er1,El1)
    CALL oneport_calc(FT(E22),FT(E21),FT(H22),FT(H21),d2,beta2,Z2,Er2,El2)
   
! we now have the incident and scattered wave amplitudes either side of the interface
! use these to calculate the appropriate S parameters for the excited port 
! then write the S parameters and the scattered power to files	 
    if (excitation_port.eq.1) then
    
      S11=El1*exp(j*beta1*l1)/( Er1*exp(-j*beta1*l1) )
      S21=Er2*exp(j*beta2*l2)/( Er1*exp(-j*beta1*l1) )
      Pr=abs(s11)**2
      Pt=abs(s21)**2
      Pa=1d0-Pr-Pt
      
      write(local_file_unit,8000)frequency,dble(S11),dimag(S11),abs(S11),20.0*log10(abs(S11)),	&
                                           dble(S21),dimag(S21),abs(S21),20.0*log10(abs(S21)),Pr,Pt,Pa

!      write(local_file_unit2,8000)frequency,dble(beta1),dimag(beta1),dble(beta2),dimag(beta2)
      
    else if (excitation_port.eq.2) then
    
      S12=El1*exp(j*beta1*l1)/( El2*exp(-j*beta2*l2) )
      S22=Er2*exp(j*beta2*l2)/( El2*exp(-j*beta2*l2) )
      Pr=abs(s22)**2
      Pt=abs(s12)**2
      Pa=1d0-Pr-Pt
      
      write(local_file_unit,8000)frequency,dble(S12),dimag(S12),abs(S12),20.0*log10(abs(S12)),	&
                                           dble(S22),dimag(S22),abs(S22),20.0*log10(abs(S22)),Pr,Pt,Pa

!      write(local_file_unit2,8000)frequency,dble(beta1),dimag(beta1),dble(beta2),dimag(beta2)
      
    end if

8000     format(12E14.5)
    
  end do ! next frequency value
  
  CLOSE(unit=local_file_unit)
  
!  CLOSE(unit=local_file_unit2)

  CALL Deallocate_post_data()

  RETURN
  
END SUBROUTINE twoport
