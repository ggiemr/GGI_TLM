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
!     SUBROUTINE add_dielectric_layers
!
! DESCRIPTION
!	
!     
! COMMENTS
!     
!
!
!
! HISTORY
!
!     started 14/3/19 CJS
!     
!

SUBROUTINE add_dielectric_layers

USE PCB_simulation

IMPLICIT NONE

! local variables

character(LEN=256) :: material_name

real :: dxmin,dxmax,dymin,dymax,dzmin,dzmax

integer :: i

! START

  read(10,*)n_dielectrics
  
  do i=1,n_dielectrics
    
    read(10,*)dxmin,dymin,dzmin,dxmax,dymax,dzmax
  
    read(10,'(A)')material_name

! Add the dielectric block to the volume geometry list

    n_volumes=n_volumes+1
    
    if(n_volumes.GT.max_volumes) then
    
      write(*,*)'ERROR in GGI_TLM_create_PCB_simulation_model: maximum number of volumes exceeded'
      write(*,*)'Maximum number of volumes is set to ',max_volumes
      write(*,*)'in /GGI_TLM/SRC/TLM_MODULES/PCB_simulation_setup_modules.F90'
      STOP 1
    
    end if

    volume_type(n_volumes)=volume_type_rectangular_block2
    
    volume_parameters(n_volumes,1)=dxmin
    volume_parameters(n_volumes,2)=dymin
    volume_parameters(n_volumes,3)=dzmin
    volume_parameters(n_volumes,4)=dxmax
    volume_parameters(n_volumes,5)=dymax
    volume_parameters(n_volumes,6)=dzmax
    
! Volume material stuff

    n_volume_materials=n_volume_materials+1
    volume_material_type(n_volume_materials)=volume_material_type_DISPERSIVE
    volume_material_name(n_volume_materials)=trim(material_name)
    volume_material_to_volume_list(n_volume_materials)=n_volumes         

  end do

RETURN  
  
END SUBROUTINE add_dielectric_layers
