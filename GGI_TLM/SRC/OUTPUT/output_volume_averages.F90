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
! SUBROUTINE set_output_volume_averages_in_mesh
! SUBROUTINE initialise_output_volume_averages
! SUBROUTINE cell_output_volume_averages
!
! NAME
!     set_output_volume_averages_in_mesh
!
! DESCRIPTION
!     loop over all the required outputs and flag all the output cells 
!     in the array local_cell_output(i,j,k) 
!     
! COMMENTS
!
! HISTORY
!
!     started 11/02/2013 CJS
!
!
SUBROUTINE set_output_volume_averages_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE file_information

IMPLICIT NONE

! local variables

  integer	:: output_volume
  integer	:: volume_number
  integer	:: number_of_cells
  integer	:: cx,cy,cz
  integer	:: output_cell

! START
  
  CALL write_line('CALLED: set_output_volume_averages_in_mesh',0,output_to_screen_flag)
    
  if (n_output_volume_averages.gt.0) then

! OUTPUT volumeS
    do output_volume=1,n_output_volume_averages
    
      volume_number=output_volume_averages(output_volume)%volume_number
      number_of_cells=problem_volumes(volume_number)%number_of_cells
      output_volume_averages(output_volume)%number_of_cells=number_of_cells
      
! allocate cell list for the output volume      
      ALLOCATE( output_volume_averages(output_volume)%cell_list(1:number_of_cells) )
      
      do output_cell=1,number_of_cells
      
! copy the cells from the geometry structure to the local
! output_volume structure
	
        cx=problem_volumes(volume_number)%cell_list(output_cell)%cell%i
        cy=problem_volumes(volume_number)%cell_list(output_cell)%cell%j
        cz=problem_volumes(volume_number)%cell_list(output_cell)%cell%k
 
        if (rank.eq.cell_rank(cz)) then
! output point belongs to this processor
	  
          local_cell_output(cx,cy,cz)=1 
  
        end if ! output point belongs to this processor
	
        output_volume_averages(output_volume)%cell_list(output_cell)=	&
	       problem_volumes(volume_number)%cell_list(output_cell)%cell
                 
      end do !next cell in this volume	  
  
    end do ! next output volume

  end if ! n_output_volume_averages.GT.0

  
END SUBROUTINE set_output_volume_averages_in_mesh
!
! NAME
!     initialise_output_volume_averages
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 11/02/2013 CJS
!
!
SUBROUTINE initialise_output_volume_averages

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

  integer	:: output_volume
  integer	:: cx,cy,cz
  integer	:: number_of_cells
  integer	:: output_cell
  
! START

  
  CALL write_line('CALLED: initialise_output_volume_averages',0,output_to_screen_flag)
  
  if (n_output_volume_averages.eq.0) RETURN
      
! LOOP OVER OUTPUT VOLUMES SETTING MESH DATA
  do output_volume=1,n_output_volume_averages
  
    number_of_cells=output_volume_averages(output_volume)%number_of_cells
    
    ALLOCATE( output_volume_averages(output_volume)%cell_output_field_number_list(1:number_of_cells) )
    
    output_volume_averages(output_volume)%cell_output_field_number_list(1:number_of_cells)=0
    
    do output_cell=1,number_of_cells
  
      cx  =output_volume_averages(output_volume)%cell_list(output_cell)%i
      cy  =output_volume_averages(output_volume)%cell_list(output_cell)%j
      cz  =output_volume_averages(output_volume)%cell_list(output_cell)%k

      if (rank.eq.cell_rank(cz)) then
      
  	output_volume_averages(output_volume)%cell_output_field_number_list(output_cell)= &
        		    local_cell_output(cx,cy,cz) 
           
      end if ! cell belongs to this process
        	
    end do !next cell in this volume	
  
! set output timestep information
    CALL set_output_time_information(output_volume_averages(output_volume)%specified_timestep_information,  &
  				     output_volume_averages(output_volume)%first_timestep,   &
  				     output_volume_averages(output_volume)%last_timestep,    &
  				     output_volume_averages(output_volume)%timestep_interval, 	  &
  				     output_volume_averages(output_volume)%specified_time_information,	  &
  				     output_volume_averages(output_volume)%first_time,	  &		
  				     output_volume_averages(output_volume)%last_time, 	  &
  				     output_volume_averages(output_volume)%time_interval,   &
  				     output_volume_averages(output_volume)%number_of_output_timesteps )
  
  end do ! next output volume

  if (rank.eq.0) then
! rank 0 process only: write header for volume average field outputs
  
    OPEN(unit=volume_average_field_output_unit,file=trim(problem_name)//volume_average_field_output_extn)
  
    CALL write_time_domain_header_data(volume_average_field_output_unit,n_output_volume_averages,n_timesteps)
     
  end if

  RETURN

END SUBROUTINE initialise_output_volume_averages
!
! SUBROUTINE cell_output_volume_averages
!
! NAME
!     cell_output_volume_averages
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 11/03/2013 CJS
!
!
SUBROUTINE cell_output_volume_averages

USE TLM_general
USE mesh
USE TLM_output
USE file_information
USE output_formats
USE cell_parameters
USE constants

IMPLICIT NONE

! local variables

  integer 	:: output_volume
  integer 	:: output_cell
  integer 	:: output_field_number
  integer 	:: field_component
    
  logical	:: output_flag
  
  integer	:: n_data
  integer	:: loop
  integer	:: cz
 
  integer			:: local_n_cells
  real*8			:: local_value

! START
  
  CALL write_line('CALLED: cell_output_volume_averages',0,timestepping_output_to_screen_flag)

  do output_volume=1,n_output_volume_averages
  
    CALL get_output_flag(output_flag,	&
                         output_volume_averages(output_volume)%first_timestep,	&
                         output_volume_averages(output_volume)%last_timestep,	&
                         output_volume_averages(output_volume)%timestep_interval )
			 
    output_volume_averages(output_volume)%value=0d0  
			 
    if (output_flag) then
    

! loop over the output volume cells in this process and set the value
      do output_cell=1,output_volume_averages(output_volume)%number_of_cells
         
        cz  =output_volume_averages(output_volume)%cell_list(output_cell)%k

        if (rank.eq.cell_rank(cz)) then
	 
          output_field_number=output_volume_averages(output_volume)%cell_output_field_number_list(output_cell)
          field_component=output_volume_averages(output_volume)%field_component  

          output_volume_averages(output_volume)%value=output_volume_averages(output_volume)%value+	&
	                  cell_output_field(output_field_number,field_component)
        end if
	
      end do ! next cell
	 
#if defined(MPI)

! Sum the number of cell in each process to rank 0
      n_data=1
      call MPI_REDUCE(output_volume_averages(output_volume)%number_of_cells, local_n_cells, &
                      n_data,MPI_INTEGER, MPI_SUM, 0,MPI_COMM_WORLD,ierror)

! Sum the field values in each process to rank 0
      n_data=1
      call MPI_REDUCE(output_volume_averages(output_volume)%value, local_value, &
                      n_data,MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD,ierror)

#elif defined(SEQ)

      local_n_cells=output_volume_averages(output_volume)%number_of_cells
      local_value=output_volume_averages(output_volume)%value

#endif

! The rank 0 process writes the data to file

      if (rank.eq.0) then
	  
	if (local_n_cells.NE.0) then
	  local_value=local_value/dble(local_n_cells)
	else
	  local_value=0d0
	end if
	    
        if ( abs(local_value).lt.1D-30 )local_value=0d0
        write(volume_average_field_output_unit,time_domain_output_format)time,output_volume,local_value

      end if ! rank.eq.0
	  
    end if ! output_flag=TRUE
  
  end do ! next output_volume
  
  CALL write_line('FINISHED: cell_output_volume_averages',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE cell_output_volume_averages

