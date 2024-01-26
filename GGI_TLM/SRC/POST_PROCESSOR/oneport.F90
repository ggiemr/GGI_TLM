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
! SUBROUTINE oneport
! SUBROUTINE oneport_calc
!
! NAME
!    oneport
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
!
!
SUBROUTINE oneport


USE post_process
USE file_information
USE constants

IMPLICIT NONE

! local variables

  integer	:: function_number
  integer	:: n_frequencies
  integer	:: frequency_loop
  
  real*8	:: frequency,w
  
  real*8	:: f1,f2,f3,f4
  
  character(len=256)	:: filename
  
  real*8	:: d_port
  real*8	:: d_interface
       
  real*8	:: d,l
  
  complex*16	:: FT(4)
  
  complex*16	:: beta,Z,R,Ei,Er,Zs
  
  integer,parameter :: E1=1
  integer,parameter :: E2=2
  integer,parameter :: H1=3
  integer,parameter :: H2=4
  
! START

  write(*,*)'Oneport analysis'

  write(*,*)'Port configuration:'
  write(*,*)''
  write(*,*)'  op1   PORT_1  op2            | scatterer | '
  write(*,*)''
       
! Allocate frequency domain data   
  
  n_functions_of_time=0
  n_functions_of_frequency=4
  
  CALL Allocate_post_data()
  
! read input files     
  write(*,*)'Read the file for E field at port 1, op1'
  function_number=E1
  CALL read_frequency_domain_data(function_number)
  
  write(*,*)'Read the file for E field at port 1, op2'
  function_number=E2
  CALL read_frequency_domain_data(function_number)
  
  write(*,*)'Read the file for H field at port 1, op1'
  function_number=H1
  CALL read_frequency_domain_data(function_number)
  
  write(*,*)'Read the file for H field at port 1, op2'
  function_number=H2
  CALL read_frequency_domain_data(function_number)
  
! check that the frequencies match...

  if ( (function_of_frequency(1)%n_frequencies.NE.function_of_frequency(2)%n_frequencies).OR.	&
       (function_of_frequency(1)%n_frequencies.NE.function_of_frequency(3)%n_frequencies).OR.	&
       (function_of_frequency(1)%n_frequencies.NE.function_of_frequency(4)%n_frequencies) ) then
  
    write(*,*)'Frequency mismatch between functions'
    write(*,*)'n_frequencies, f1=',function_of_frequency(1)%n_frequencies
    write(*,*)'n_frequencies, f2=',function_of_frequency(2)%n_frequencies
    write(*,*)'n_frequencies, f3=',function_of_frequency(3)%n_frequencies
    write(*,*)'n_frequencies, f4=',function_of_frequency(4)%n_frequencies
    STOP
    
  end if
  
  do frequency_loop=1,function_of_frequency(1)%n_frequencies  
    
    f1=function_of_frequency(1)%frequency(frequency_loop)
    f2=function_of_frequency(2)%frequency(frequency_loop)
    f3=function_of_frequency(3)%frequency(frequency_loop)
    f4=function_of_frequency(4)%frequency(frequency_loop)
  
    if ( (f1.NE.f2).OR.(f1.NE.f3).OR.(f1.NE.f4) ) then	 
      write(*,*)'Frequency mismatch in functions f1, f2, f3 and f4'
      write(*,*)'frequency number',frequency_loop
      write(*,*)'frequency, f1=',f1
      write(*,*)'frequency, f2=',f2
      write(*,*)'frequency, f3=',f3
      write(*,*)'frequency, f4=',f4
      STOP     
    end if  
    
  end do ! next frequency to check
  
! Read other solution data  
        
  write(*,*)'Enter port separation (m)'
  read(*,*)d_port
  write(record_user_inputs_unit,*)d_port,' port separation (m)'
  d=d_port/2d0
  
  write(*,*)'Enter port midpoint to interface distamce (m)'
  read(*,*)d_interface
  write(record_user_inputs_unit,*)d_interface,' port midpoint to interface distamce (m)'
  l=d_interface
 
! File for results

  write(*,*)'Enter the frequency domain data output filename'
  read(*,*)filename
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  OPEN(unit=local_file_unit,file=trim(filename))
  OPEN(unit=local_file_unit2,file=trim(filename)//'.Zs')

  n_frequencies=function_of_frequency(1)%n_frequencies
  
  write(local_file_unit,'(A,A)')	&
	'    frequency       Re{Beta}     Im{Beta}       Re{Z}          Im{Z}',	&
	'         Re{R}          Im{R}        |R|          |R|(dB)      Re{w/Beta}    1-|R|**2'
 
  write(local_file_unit2,'(A)')'    frequency        Re{Zs}       Im{Zs}'
  
  do frequency_loop=1,n_frequencies
  
    frequency=function_of_frequency(1)%frequency(frequency_loop)
    w=2D0*pi*frequency
 
    FT(1)=function_of_frequency(1)%value(frequency_loop)
    FT(2)=function_of_frequency(2)%value(frequency_loop)
    FT(3)=function_of_frequency(3)%value(frequency_loop)
    FT(4)=function_of_frequency(4)%value(frequency_loop)
    
    CALL oneport_calc(FT(E1),FT(E2),FT(H1),FT(H2),d,beta,Z,Ei,Er)
        
    R=(Er/Ei)*exp(2d0*j*beta*l)
    
    Zs= Z*(1d0+R)/(1d0-R)
    
    write(local_file_unit,8000)frequency,dble(beta),dimag(beta),dble(Z),dimag(Z),	&
                                         dble(R),dimag(R),abs(R),20d0*log10(abs(R)),    &
					 dble(w/beta),1d0-abs(R)**2
    
    write(local_file_unit2,8010)frequency,dble(zs),dimag(Zs),abs(Zs) 

8000     format(11E14.5)
8010     format(4E14.5)
    
  end do ! next frequency value
  
  CLOSE(unit=local_file_unit)
  CLOSE(unit=local_file_unit2)

  CALL Deallocate_post_data()

  RETURN
  

  
END SUBROUTINE oneport
!
!
!
  SUBROUTINE oneport_calc(E1,E2,H1,H2,d,beta,Z,Er,El)
 
USE constants 
  
  complex*16 E1,E2,H1,H2,beta,Z,Er,El
  real*8 d

! local variables

  complex*16 a,b,c,s,beta1,beta2,P,P1,P2,Eit,Ert

! START
  
  a=E2*H1+E1*H2
  b=-2d0*(E1*H1+E2*H2)
  c=E2*H1+E1*H2

  s=sqrt(b*b-4d0*a*c)

  p1=(-b+s)/(2d0*a)
  beta1=log(p1)/(2d0*j*d)

  p2=(-b-s)/(2d0*a)
  beta2=log(p2)/(2d0*j*d)

  beta=beta1
  
  if (dble(beta).lt.0d0) beta=-beta1
  
  P=exp(j*beta*d)
  Eit=(E2*P-E1*P*P*P)/(1d0-P*P*P*P)
  Z=(2d0*Eit*P-E1)/H1
  Ert=-Eit*P*P+E1*P

! check sign of Z and reverse if necessary

  if (dble(Z).lt.0d0) Z=-Z

  El=Ert
  Er=Eit

  return
  
  END SUBROUTINE oneport_calc
