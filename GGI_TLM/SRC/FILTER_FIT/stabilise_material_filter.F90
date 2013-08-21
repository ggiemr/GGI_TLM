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
!     stabilise_dielectric_or_magnetic_material
!
! DESCRIPTION
!     Given the permittivity / permeability filter function in pole-residue form:
!     1. Ensure that all poles are stable (LHS of s plane)
!     2. Ensure that Im{epsr}<=0 for all testing frequencies
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 13/12/2012 CJS
!
!
SUBROUTINE stabilise_dielectric_or_magnetic_material

USE FF_parameters
USE FF_general
USE FF_input_data
USE FF_filters
USE FF_file_stuff

IMPLICIT NONE
 
! local variables

  integer	:: pole
  real*8	:: ap,app
  real*8	:: cp,cpp
  real*8	:: test

! START

!  write(*,*)'Stabilising Material model'
       
! check whether the conductivity is positive
  if(filter_sigma(1).lt.0d0) then	
!    write(*,*)'Stabilising sigma',filter_sigma(1)
! new estimate for conductivity
    filter_sigma(1)=-filter_sigma(1)/10d0
  end if
		
! check whether the constant term (value as w-> infinity) is 
! less than 1. If so set it to 1.1 as a start for the optimisation. 
  if(filter_S_PR(1)%C.lt.1d0) then    
!    print*,'Stabilising C:',filter_S_PR(1)%C
    filter_S_PR(1)%C=1.1d0
  end if
          
! loop over poles     
  pole=1
10  CONTINUE       

    if (.NOT.filter_S_PR(1)%complex_pole(pole)) then
	 
! real pole. We require that the residues of real poles are positive
      if (dble(filter_S_PR(1)%residues(pole)).lt.0d0) then
!	print*,'Stabilising residue:',filter_S_PR(1)%residues(pole)
! unstable residue so make positive to stabilise the system
        filter_S_PR(1)%residues(pole)=-filter_S_PR(1)%residues(pole)
      end if	
	   
      pole=pole+1   
	   
    else
! complex pole pair 

      ap=  dble(filter_S_PR(1)%poles(pole))
      app=dimag(filter_S_PR(1)%poles(pole))
      cp=  dble(filter_S_PR(1)%residues(pole))
      cpp=dimag(filter_S_PR(1)%residues(pole))

! 1. ensure that real part of zero>0
      if (cp.lt.0D0) then
!        print*,'Stabilising pole, 1:',cp
        cp=-cp
      end if	      

! 2. condition on imaginary part of zero...
      test=(cp*(app*app-ap*ap))/(2D0*ap*app)
      if (app.lt.0D0) then
        if (cpp.lt.test) then
!          print*,'Stabilising pole, 2:',cpp,test
          cpp=test
        end if
      else
        if (cpp.gt.test) then
!          print*,'Stabilising pole, 2:',cpp,test
          cpp=test
        end if
      end if	      

      filter_S_PR(1)%poles(pole  )=cmplx(ap, app)
      filter_S_PR(1)%poles(pole+1)=cmplx(ap,-app)
      filter_S_PR(1)%residues(pole  )=cmplx(cp, cpp)
      filter_S_PR(1)%residues(pole+1)=cmplx(cp,-cpp)

      pole=pole+2  

    end if ! complex pole

  if (pole.le.order) goto 10

  RETURN

END SUBROUTINE stabilise_dielectric_or_magnetic_material
