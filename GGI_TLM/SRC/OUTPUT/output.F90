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
! SUBROUTINE set_outputs_in_mesh
! SUBROUTINE initialise_outputs
! SUBROUTINE cell_output
! SUBROUTINE face_output
! SUBROUTINE cable_output
! SUBROUTINE write_frequency_domain_outputs
! SUBROUTINE finish_output
! SUBROUTINE set_output_time_information
! SUBROUTINE get_output_flag
!
! NAME
!     set_outputs_in_mesh
!
! DESCRIPTION
!     loop over all the required outputs and flag all the output cells/ faces 
!     in the arrays local_cell_output(i,j,k) or local_surface_output(i,j,k,face)
!     as required
!     
! COMMENTS
!     Output_surfaces could be made more efficient in parallel - every process gets the whole output plane...
!
! HISTORY
!
!     started 14/08/2012 CJS
!     output_surfaces included 14/09/2012 CJS
!     Parallel 23/11/2012 CJS
!
!
SUBROUTINE set_outputs_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE file_information

IMPLICIT NONE

! local variables

! START

  
  CALL write_line('CALLED: set_outputs_in_mesh',0,output_to_screen_flag)
  
  CALL set_output_points_in_mesh()
  
  CALL set_output_surfaces_in_mesh()
  
  CALL set_output_volumes_in_mesh()
  
  CALL set_output_volume_averages_in_mesh()
  
  CALL set_mode_outputs_in_mesh()
  
  CALL set_frequency_output_surfaces_in_mesh()
  
  CALL set_frequency_output_volumes_in_mesh()
  
  CALL set_frequency_domain_power_surfaces_in_mesh()
  
  CALL set_far_field_surfaces_in_mesh()
  
  CALL set_RCS_surfaces_in_mesh()
  
  CALL set_SAR_volumes_in_mesh()
  
  CALL write_line('FINISHED: set_outputs_in_mesh',0,output_to_screen_flag)

  RETURN

END SUBROUTINE set_outputs_in_mesh
!
! NAME
!     initialise_outputs
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!
!
SUBROUTINE initialise_outputs

USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE Cables
USE file_information

IMPLICIT NONE

! local variables

! START

  
  CALL write_line('CALLED: initialise_outputs',0,output_to_screen_flag)
  
  CALL initialise_output_points()
  
  CALL initialise_output_surfaces()
  
  CALL initialise_output_volumes()
  
  CALL initialise_output_volume_averages()
  
  CALL initialise_mode_outputs()
  
  CALL initialise_frequency_output_surfaces()
  
  CALL initialise_frequency_output_volumes()
  
  CALL initialise_frequency_domain_power_surfaces()
  
  CALL initialise_far_field_surfaces()
  
  CALL initialise_RCS_surfaces()
  
  CALL initialise_SAR_volumes()
  
  CALL initialise_cable_outputs()
  
  CALL write_line('FINISHED: initialise_outputs',0,output_to_screen_flag)

  RETURN

END SUBROUTINE initialise_outputs
!
! SUBROUTINE cell_output
!
! NAME
!     cell_output
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!     parallel 23/08/2012 CJS
!
!
SUBROUTINE cell_output

USE TLM_general
USE mesh
USE TLM_output
USE file_information
USE output_formats

IMPLICIT NONE

! local variables

! START
  
  CALL write_line('CALLED: cell_output',0,timestepping_output_to_screen_flag)

  CALL cell_output_point()

  CALL cell_output_volumes()

  CALL cell_output_volume_averages()
  
  CALL cell_output_frequency_output_volumes()

  CALL cell_output_SAR_volumes()
  
  CALL write_line('FINISHED: cell_output',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE cell_output
!
! SUBROUTINE face_output
!
! NAME
!     face_output
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!
!
SUBROUTINE face_output

USE TLM_general
USE mesh
USE TLM_output
USE file_information
USE output_formats
USE cell_parameters
USE constants

IMPLICIT NONE

! local variables

! START
  
  CALL write_line('CALLED: face_output',0,timestepping_output_to_screen_flag)

  CALL face_output_points()

  CALL face_output_surfaces()

  CALL mode_output()
  
  CALL face_output_frequency_output_surfaces()
  
  CALL face_output_frequency_domain_power_surfaces()
  
  CALL face_output_far_field()
  
  CALL face_output_RCS()

  CALL write_line('FINISHED: face_output',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE face_output
!
! NAME
!     write_frequency_domain_outputs
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 5/12/2012 CJS
!
!
SUBROUTINE write_frequency_domain_outputs

USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE Cables
USE file_information

IMPLICIT NONE

! local variables

! START

  
  CALL write_line('CALLED: write_frequency_domain_outputs',0,output_to_screen_flag)
    
  CALL write_frequency_output_surfaces()
    
  CALL write_frequency_output_volumes()
    
  CALL write_frequency_domain_power_surfaces()
  
  CALL write_far_field_surfaces()
  
  CALL write_RCS_surfaces()
  
  CALL write_SAR_volumes()
  
  CALL write_line('FINISHED: write_frequency_domain_outputs',0,output_to_screen_flag)

  RETURN

END SUBROUTINE write_frequency_domain_outputs

!
! NAME
!     finish_outputs
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 15/08/2012 CJS
!
!
SUBROUTINE finish_outputs

USE TLM_general
USE mesh
USE TLM_output
USE Cables
USE file_information

IMPLICIT NONE

! local variables

! START
  
  CALL write_line('CALLED: finish_outputs',0,output_to_screen_flag)
  
  if (n_output_points.gt.0) then
  
    CLOSE(unit=field_output_unit)
    
  end if
  
  if (n_output_surfaces.gt.0) then
  
    CLOSE(unit=surface_field_output_unit)
    
  end if
  
  if (n_output_volumes.gt.0) then
  
    CLOSE(unit=volume_field_output_unit)
    
  end if
  
  if (n_output_volume_averages.gt.0) then
  
    CLOSE(unit=volume_average_field_output_unit)
    
  end if
  
  if (n_output_modes.gt.0) then
  
    CLOSE(unit=mode_output_unit)
    
  end if
  
  if (n_far_field_surfaces.gt.0) then
  
    CLOSE(unit=far_field_output_unit)
    
  end if
  
  if (n_RCS_surfaces.gt.0) then
  
    CLOSE(unit=rcs_output_unit)
    
  end if
  
  if (n_SAR_volumes.gt.0) then
  
    CLOSE(unit=SAR_output_unit)
    
  end if
  
  if (n_frequency_output_surfaces.gt.0) then
  
    CLOSE(unit=frequency_output_surface_unit)
    
  end if
  
  if (n_frequency_output_volumes.gt.0) then
  
    CLOSE(unit=frequency_output_volume_unit)
    
  end if

! CABLE OUTPUTS
  if (n_cable_outputs.gt.0) then
  
    CLOSE(unit=cable_current_output_unit)
  
  end if  ! n_cable_outputs.gt.0
 
  CALL write_line('FINISHED: finish_outputs',0,output_to_screen_flag)

  RETURN

END SUBROUTINE finish_outputs
!
! NAME
!     set_output_time_information
!
! DESCRIPTION
!     set the output time information if not already specified i.e.
!     first timestep for output
!     last timestep for output
!     interval between outputs
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 13/09/2012 CJS
!
!
SUBROUTINE set_output_time_information( specified_timestep_information,	   &
                                        first_timestep,	   &
                                        last_timestep,	   &
                                        timestep_interval,   &
                                        specified_time_information,  &
                                        first_time,  &					 
                                        last_time,   &				 
                                        time_interval,	&
					number_of_output_timesteps )

USE TLM_general

IMPLICIT NONE

  logical		:: specified_timestep_information
  integer		:: first_timestep
  integer		:: last_timestep
  integer		:: timestep_interval
  
  logical		:: specified_time_information
  real*8		:: first_time
  real*8		:: last_time
  real*8		:: time_interval
  
  integer		:: number_of_output_timesteps

! local variables

  logical		:: output_flag
  
! START  
    
    if (specified_time_information) then   
! we need to convert time information into timestep information

      first_timestep=NINT(first_time/dt)
      last_timestep=NINT(last_time/dt)
      timestep_interval=NINT(time_interval/dt)
      
    end if

! some checks here to make sure the numbers are sensible
    if (first_timestep.lt.1) first_timestep=1
    if (last_timestep.lt.1) last_timestep=1
    if (timestep_interval.lt.1) timestep_interval=1

! simulate the timestep loop to calculate the number of output timesteps    
    number_of_output_timesteps=0
    do timestep=1,n_timesteps
    
      CALL get_output_flag(output_flag,first_timestep,last_timestep,timestep_interval )
			 
      if (output_flag) number_of_output_timesteps=number_of_output_timesteps+1

    end do
    
  RETURN
  
END SUBROUTINE set_output_time_information
!
! NAME
!     get_output_flag
!
! DESCRIPTION
!     given the current timestep (in TLM_general) plus: first timestep for output,
!     last timestep for output and interval between outputs, return a flag to indicate
!     whether we should write output at this timestep
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 13/09/2012 CJS
!
!
SUBROUTINE get_output_flag( output_flag,	&
			    first_timestep,	   &
                            last_timestep,	   &
                            timestep_interval )

USE TLM_general

IMPLICIT NONE

  logical		:: output_flag
  integer		:: first_timestep
  integer		:: last_timestep
  integer		:: timestep_interval

! local variables

  integer		:: timesteps_from_first

! START  

  output_flag=.FALSE.
  
  if ( (timestep.lt.first_timestep).OR.(timestep.gt.last_timestep) ) RETURN

  timesteps_from_first=timestep-first_timestep
  
  if (mod(timesteps_from_first,timestep_interval).EQ.0) then
    output_flag=.TRUE.
    RETURN
  end if
    
  RETURN
  
END SUBROUTINE get_output_flag
