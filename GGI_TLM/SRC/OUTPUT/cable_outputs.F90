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
!SUBROUTINE initialise_cable_outputs
!SUBROUTINE cable_output
!
!
! NAME
!     initialise_cable_outputs
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
SUBROUTINE initialise_cable_outputs

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

  CALL write_line('CALLED: initialise_cable_outputs',0,output_to_screen_flag)

  if (n_cable_outputs.gt.0) then
  
    if (rank.eq.0) then
! rank 0 process only:

      OPEN(unit=cable_current_output_unit,file=trim(problem_name)//cable_current_output_extn)
    
      CALL write_time_domain_header_data(cable_current_output_unit,n_cable_outputs,n_timesteps)
 
    end if ! rank 0 process 
  
  end if  ! n_cable_outputs.gt.0
  
  CALL write_line('FINISHED: initialise_cable_outputs',0,output_to_screen_flag)

  RETURN

END SUBROUTINE initialise_cable_outputs
!
! SUBROUTINE cable_output
!
! NAME
!     cable_output
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
SUBROUTINE cable_output

USE TLM_general
USE TLM_output
USE file_information
USE output_formats
USE Cables
USE cell_parameters

IMPLICIT NONE

! local variables

  integer	:: output_number
  integer	:: output_cable
  integer	:: output_segment
  integer	:: output_point
  
  integer	:: n_conductors,conductor,bundle_segment_conductor
  
  integer	:: output_count
  
  real*8	:: output_time_face
  real*8	:: output_time_centre
  
  real*8	:: value
  
  integer 	:: i,p_rank
  integer 	:: talk_to,nreal
  integer	:: cz

! START
  
  CALL write_line('CALLED: cable_output',0,timestepping_output_to_screen_flag)

  output_time_centre=(timestep-1)*dt
  output_time_face=output_time_centre+dt/2d0
  
  output_count=0

  do output_number=1,n_cable_outputs

    output_cable=cable_output_list(output_number)%cable_number  
    output_segment=cable_output_list(output_number)%bundle_segment_number  
    output_point=cable_output_list(output_number)%output_point%point
    
    cz=cable_output_list(output_number)%output_point%cell%k
    
    n_conductors=cable_output_list(output_number)%n_conductors
    
    do conductor=1,n_conductors
    
      output_count=output_count+1
      
      bundle_segment_conductor=cable_output_list(output_number)%conductor_list(conductor)

! Calculate the output value in the process which owns the output cell 
      if (rank.eq.cell_rank(cz)) then
      
        if (output_point.EQ.centre) then
! use the cell centre current on the segment   
          value=bundle_segment_list(output_segment)%Iw_centre(bundle_segment_conductor)    
        else
! use the face current on the segment   
          value=bundle_segment_list(output_segment)%Iw_face(bundle_segment_conductor)    
        end if
      end if
      
#if defined(MPI)

! Make sure that the output current value gets to the rank 0 process
      if (rank.ne.0) then
  
        p_rank=cell_rank(cz)
	  
        if (p_rank.eq.rank) then
! send output data to rank 0 processor
          talk_to=0
          nreal=1
          CALL MPI_SEND(value,nreal,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)	
!          write(*,8000)'rank',rank,' sends to rank ',talk_to,' value',value,nreal
        end if ! output cell belongs to this process 
      
      else ! this is the rank 0 process
      
        p_rank=cell_rank(cz)
        talk_to=p_rank
        nreal=1
      
        if (p_rank.ne.0) then
	
! get the value from another process
          CALL MPI_RECV(value,nreal,MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
!          write(*,8000)'rank',rank,' gets from rank ',talk_to,'  value',value,nreal

        end if ! value must be read from another process
      
      end if ! rank 0 process

#endif

      if (rank.eq.0) then
      
        if (output_point.EQ.centre) then
! use the cell centre current on the segment   
          if (abs(value).lt.1d-99) value=0d0
          write(cable_current_output_unit,cable_current_output_format)output_time_centre,output_count,	&
	                                                              value,output_cable,conductor   
        else
! use the face current on the segment    
          if (abs(value).lt.1d-99) value=0d0
          write(cable_current_output_unit,cable_current_output_format)output_time_face,output_count,	&
	                                                              value,output_cable,conductor 
        end if ! centre or face current	
      
      end if ! rank 0 process
     
    end do ! next cable conductor

  end do ! next output number
  
  CALL write_line('FINISHED: cable_output',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE cable_output
