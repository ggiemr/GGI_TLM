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
! SUBROUTINE check_solver_input_data
!
! NAME
!     check_solver_input_data
!
! DESCRIPTION
!     check that all the data required for the problem is present in the input file
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 12/10/2012 CJS
!
!
SUBROUTINE check_solver_input_data

USE TLM_general
USE TLM_volume_materials
USE TLM_surface_materials
USE TLM_excitation
USE TLM_output
USE Cables
USE geometry
USE file_information

IMPLICIT NONE

! local variables

  integer	:: i,j,k
  integer	:: volume
  integer	:: surface
  integer	:: line
  integer	:: point
  integer	:: excitation_function
  integer	:: cable
  integer	:: cable_end
  integer	:: cable_geometry

! START  

  CALL write_line('CALLED: check_solver_input_data',0,output_to_screen_flag)

! Check volume material data i.e. do the volumes in the volume list exist?
  do i=1,n_volume_materials
  
    do j=1,volume_material_list(i)%n_volumes
    
      volume=volume_material_list(i)%volume_list(j)
      if ( (volume.gt.n_volumes).OR.(volume.lt.1) )then
      
        CALL write_line('ERROR in volume_material_list: Volume does not exist',0,.TRUE.)
        CALL write_line_integer('Material number:',i,0,.TRUE.)
        CALL write_line_integer('Volume number:',volume,0,.TRUE.)
        CALL write_line_integer('Number of volumes=',n_volumes,0,.TRUE.)
	STOP
	
      end if
      
    end do  ! next volume in volume list
    
  end do ! next volume material

! Check surface material data i.e. do the surfaces in the surface list exist?
  do i=1,n_surface_materials
  
    do j=1,surface_material_list(i)%n_surfaces
    
      surface=surface_material_list(i)%surface_list(j)
      if ( (surface.gt.n_surfaces).OR.(surface.lt.1) )then
      
        CALL write_line('ERROR in surface_material_list: surface does not exist',0,.TRUE.)
        CALL write_line_integer('Material number:',i,0,.TRUE.)
        CALL write_line_integer('surface number:',surface,0,.TRUE.)
        CALL write_line_integer('Number of surfaces=',n_surfaces,0,.TRUE.)
	STOP
	
      end if
      
    end do ! next surface in surface list
    
  end do ! next surface material
  
! Check excitation surfaces i.e. do the surfaces exist?  
  do i=1,n_excitation_surfaces
  
    excitation_function=excitation_surfaces(i)%excitation_function_number
! note excitaion function can be 0 i.e. no excitation
    if ( (excitation_function.gt.n_excitation_functions).OR.(excitation_function.lt.0) )then 
    
      CALL write_line('ERROR in excitation_surfaces: excitation_function does not exist',0,.TRUE.)
      CALL write_line_integer('Excitation surface number:',i,0,.TRUE.)
      CALL write_line_integer('excitation function number:',excitation_function,0,.TRUE.)
      CALL write_line_integer('Number of excitation functions=',n_excitation_functions,0,.TRUE.)
      STOP

    end if
    
    surface=excitation_surfaces(i)%surface_number
    if ( (surface.gt.n_surfaces).OR.(surface.lt.1) )then
    
      CALL write_line('ERROR in excitation_surfaces: surface does not exist',0,.TRUE.)
      CALL write_line_integer('Excitation surface number:',i,0,.TRUE.)
      CALL write_line_integer('surface number:',surface,0,.TRUE.)
      CALL write_line_integer('Number of surfaces=',n_surfaces,0,.TRUE.)
      STOP

    end if
    
  end do ! next excitation surface 

! Check excitation points i.e. do the points exist?  
  do i=1,n_excitation_points
  
    excitation_function=excitation_points(i)%excitation_function_number
! note excitaion function can be 0 i.e. no excitation
    if ( (excitation_function.gt.n_excitation_functions).OR.(excitation_function.lt.0) )then 
    
      CALL write_line('ERROR in excitation_points: excitation_function does not exist',0,.TRUE.)
      CALL write_line_integer('Excitation point number:',i,0,.TRUE.)
      CALL write_line_integer('excitation function number:',excitation_function,0,.TRUE.)
      CALL write_line_integer('Number of excitation functions=',n_excitation_functions,0,.TRUE.)
      STOP

    end if
  
    point=excitation_points(i)%point_number
    if ( (point.gt.n_points).OR.(point.lt.1) )then
    
      CALL write_line('ERROR in excitation_points: point does not exist',0,.TRUE.)
      CALL write_line_integer('Excitation point number:',i,0,.TRUE.)
      CALL write_line_integer('point number:',point,0,.TRUE.)
      CALL write_line_integer('Number of points=',n_points,0,.TRUE.)
      STOP

    end if
    
  end do ! next excitation point
  
! Check Huygens surface  
  if (n_huygens_surfaces.ne.0) then
  
    surface=huygens_surface%surface_number
    if ( (surface.gt.n_surfaces).OR.(surface.lt.0) )then
    
      CALL write_line('ERROR in Huygens_surface: surface does not exist',0,.TRUE.)
      CALL write_line_integer('surface number:',surface,0,.TRUE.)
      CALL write_line_integer('Number of surfaces=',n_surfaces,0,.TRUE.)
      STOP

    end if
 
  end if ! n_huygens_surfaces.ne.0
  
! Check Far field surface  
  if (n_far_field_surfaces.ne.0) then
  
    surface=far_field_surface%surface_number
    if ( (surface.gt.n_surfaces).OR.(surface.lt.0) )then
    
      CALL write_line('ERROR in Far Field_surface: surface does not exist',0,.TRUE.)
      CALL write_line_integer('surface number:',surface,0,.TRUE.)
      CALL write_line_integer('Number of surfaces=',n_surfaces,0,.TRUE.)
      STOP

    end if
 
  end if ! n_far_field_surfaces.ne.0
  
! Check RCS surface  
  if (n_rcs_surfaces.ne.0) then
  
    surface=rcs_surface%surface_number
    if ( (surface.gt.n_surfaces).OR.(surface.lt.0) )then
    
      CALL write_line('ERROR in RCS_surface: surface does not exist',0,.TRUE.)
      CALL write_line_integer('surface number:',surface,0,.TRUE.)
      CALL write_line_integer('Number of surfaces=',n_surfaces,0,.TRUE.)
      STOP

    end if

! we require a single Huygens surface to provide the illuminating field for RCS calculation    
    if (n_huygens_surfaces.ne.1) then
    
      CALL write_line('ERROR in RCS_surface: We require a single Huygens surface definition',0,.TRUE.)
      CALL write_line_integer('Number of Huygens surfaces=',n_huygens_surfaces,0,.TRUE.)
      STOP

    end if

  end if ! n_rcs_surfaces.ne.0
  
! Check Frequency output surface  
  if (n_frequency_output_surfaces.ne.0) then
  
    do i=1,n_frequency_output_surfaces
    
      surface=frequency_output_surface(i)%surface_number
      if ( (surface.gt.n_surfaces).OR.(surface.lt.0) )then
    
        CALL write_line('ERROR in : frequency_output_surface does not exist',0,.TRUE.)
        CALL write_line_integer('surface number:',surface,0,.TRUE.)
        CALL write_line_integer('Number of surfaces=',n_surfaces,0,.TRUE.)
        STOP

      end if
      
    end do ! next frequency_output_surface
    
  end if ! n_frequency_output_surfaces.ne.0
  
! Check Frequency output surface  
  if (n_frequency_domain_power_surfaces.ne.0) then
  
    do i=1,n_frequency_domain_power_surfaces
    
      surface=frequency_domain_power_surface(i)%surface_number
      if ( (surface.gt.n_surfaces).OR.(surface.lt.0) )then
    
        CALL write_line('ERROR in : frequency_domain_power_surface does not exist',0,.TRUE.)
        CALL write_line_integer('surface number:',surface,0,.TRUE.)
        CALL write_line_integer('Number of surfaces=',n_surfaces,0,.TRUE.)
        STOP

      end if
      
    end do ! next frequency_domain_power_surface
    
  end if ! n_frequency_domain_power_surfaces.ne.0
  
! Check output surfaces i.e. do the surfaces exist?  
  do i=1,n_output_surfaces
  
    surface=output_surfaces(i)%surface_number
    if ( (surface.gt.n_surfaces).OR.(surface.lt.1) )then
    
      CALL write_line('ERROR in output_surfaces: surface does not exist',0,.TRUE.)
      CALL write_line_integer('output surface number:',i,0,.TRUE.)
      CALL write_line_integer('surface number:',surface,0,.TRUE.)
      CALL write_line_integer('Number of surfaces=',n_surfaces,0,.TRUE.)
      STOP

    end if
    
  end do ! next output surface

! Check output points i.e. do the points exist?  
  do i=1,n_output_points
  
    point=output_points(i)%point_number
    if ( (point.gt.n_points).OR.(point.lt.1) )then
    
      CALL write_line('ERROR in output_points: point does not exist',0,.TRUE.)
      CALL write_line_integer('output point number:',i,0,.TRUE.)
      CALL write_line_integer('point number:',point,0,.TRUE.)
      CALL write_line_integer('Number of points=',n_points,0,.TRUE.)
      STOP

    end if
    
  end do ! next output point

! check that cable junction excitation functions exist
  do i=1,n_cable_junctions
      
    do j=1,cable_junction_list(i)%number_of_cables
      
      do k=1,cable_junction_list(i)%n_external_conductors(j)
    
        excitation_function=cable_junction_list(i)%excitation_function(k)%value(k)
! note excitaion function can be 0 i.e. no excitation
        if ( (excitation_function.gt.n_excitation_functions).OR.(excitation_function.lt.0) )then 
    
          CALL write_line('ERROR in cable_junction_list: excitation_function does not exist',0,.TRUE.)
          CALL write_line_integer('Cable junction number:',i,0,.TRUE.)
          CALL write_line_integer('Cable number:',j,0,.TRUE.)
          CALL write_line_integer('excitation function number:',excitation_function,0,.TRUE.)
          CALL write_line_integer('Number of excitation functions=',n_excitation_functions,0,.TRUE.)
          STOP
       
        end if
	
      end do ! next conductor on this cable
      
    end do ! next cable in cable list
    
  end do ! next cable junction
  
! check frequency warping
  if (bicubic_warp_flag.AND.frequency_scale_flag) then
  
    CALL write_line('Error: We cannont have both bicubic_warp_flag and frequency_scale_flag set',0,.TRUE.)
    STOP
   
  end if
  
! check outer boundary wrapping - reflection coefficient should be 1
  if (wrap_x) then
    if ((R_xmin.lt.0.9999d0).OR.(R_xmin.gt.1.0001d0).OR.(R_xmax.lt.0.9999d0).OR.(R_xmax.gt.1.0001d0)) then
      CALL write_line('Error X boundary: Reflection coefficient must be 1 with outer boundary wrapping on',0,.TRUE.)
      STOP
    end if
  end if

  if (wrap_y) then
    if ((R_ymin.lt.0.9999d0).OR.(R_ymin.gt.1.0001d0).OR.(R_ymax.lt.0.9999d0).OR.(R_ymax.gt.1.0001d0)) then
      CALL write_line('Error Y boundary: Reflection coefficient must be 1 with outer boundary wrapping on',0,.TRUE.)
      STOP
    end if
  end if
  
  if (wrap_z) then
    if ((R_zmin.lt.0.9999d0).OR.(R_zmin.gt.1.0001d0).OR.(R_zmax.lt.0.9999d0).OR.(R_zmax.gt.1.0001d0)) then
      CALL write_line('Error Z boundary: Reflection coefficient must be 1 with outer boundary wrapping on',0,.TRUE.)
      STOP
    end if
  end if


  CALL write_line('FINISHED: check_solver_input_data',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE check_solver_input_data
