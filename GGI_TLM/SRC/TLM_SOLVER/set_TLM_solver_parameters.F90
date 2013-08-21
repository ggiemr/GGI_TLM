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
! SUBROUTINE reset_TLM_solver_parameters
! SUBROUTINE set_TLM_solver_parameters
!
! NAME
!     reset_TLM_solver_parameters
!
! DESCRIPTION
!     reset all the parameters associated with running the TLM solver 
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 10/08/2012 CJS
!
!
SUBROUTINE reset_TLM_solver_parameters

USE TLM_general
USE TLM_excitation
USE TLM_output

IMPLICIT NONE

! local variables

! START
  
  CALL write_line('CALLED: reset_TLM_solver_parameters',0,output_to_screen_flag)

  simulation_time=0d0
  
  n_excitation_functions=0
  
  n_excitation_points=0
  
  n_output_points=0
  
  number_of_warnings=0
  
  CALL write_line('FINISHED: reset_TLM_solver_parameters',0,output_to_screen_flag)
  
  RETURN
  
END SUBROUTINE reset_TLM_solver_parameters
!
! NAME
!     set_TLM_solver_parameters
!
! DESCRIPTION
!     set all the parameters associated with running the TLM solver 
!     1. Calculate timestep
!     2. Set the number of timesteps
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!
!
SUBROUTINE set_TLM_solver_parameters

USE TLM_general
USE constants

IMPLICIT NONE

! local variables

! START
  
  CALL write_line('CALLED: set_TLM_solver_parameters',0,output_to_screen_flag)
  
  dt=dl/(2d0*c0)
  
  n_timesteps=NINT(simulation_time/dt)
  
  CALL write_line('FINISHED: set_TLM_solver_parameters',0,output_to_screen_flag)
  
  RETURN
  
END SUBROUTINE set_TLM_solver_parameters
