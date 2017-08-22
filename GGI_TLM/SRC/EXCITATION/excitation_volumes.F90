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
! SUBROUTINE set_excitation_volumes_in_mesh
!
! NAME
!     set_excitation_volumes_in_mesh
!
! DESCRIPTION
!     loop over all the required excitations and flag all the excitation cells
!     in the array local_cell_excitation(i,j,k) 
!     as required
!     
! COMMENTS
!     Excitation_volumes could be made more efficient in parallel - every process gets the whole excitation plane...
!  
!
! HISTORY
!
!     started 10/05/2017 CJS based on excitation_volumes.F90
!
!
SUBROUTINE set_excitation_volumes_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_excitation
USE file_information

IMPLICIT NONE

! local variables

  integer	:: excitation_point
  integer	:: excitation_volume
  integer	:: volume_number
  integer	:: number_of_cells
  integer	:: excitation_cell
  integer	:: n_cells
  integer	:: cx,cy,cz,face
  type(cell_point)	:: excitation_cell1
  type(cell_point)	:: excitation_cell2
    
  logical inside

! START

  
  CALL write_line('CALLED: set_excitation_volumes_in_mesh',0,output_to_screen_flag)
  
  if (n_excitation_volumes.gt.0) then

! excitation volumeS
    do excitation_volume=1,n_excitation_volumes
    
      volume_number=excitation_volumes(excitation_volume)%volume_number
      number_of_cells=problem_volumes(volume_number)%number_of_cells
      excitation_volumes(excitation_volume)%number_of_cells=number_of_cells
      
! allocate cell list for the excitation volume      
      ALLOCATE( excitation_volumes(excitation_volume)%cell_list(1:number_of_cells) )
      
      do excitation_cell=1,number_of_cells
      
! copy the cells from the geometry structure to the local
! excitation_volume structure
	
        cx=problem_volumes(volume_number)%cell_list(excitation_cell)%cell%i
        cy=problem_volumes(volume_number)%cell_list(excitation_cell)%cell%j
        cz=problem_volumes(volume_number)%cell_list(excitation_cell)%cell%k
        face=problem_volumes(volume_number)%cell_list(excitation_cell)%point
 
        if (face.eq.centre) then
! cell centre excitation
          
! Set the excitation cell in the local_cell_excitation array

          if (rank.eq.cell_rank(cz)) then
! excitation point belongs to this processor
  
            local_cell_excitation(cx,cy,cz)=1
    
!           write(*,*)'Setting cell centre excitation point',excitation_point
!           write(*,*)'Coordinates:',cx,cy,cz
  
          end if ! output point belongs to this processor
  
        else
 
          write(*,*)'ERROR, volume excitation should be at the cell centre'
          STOP 
   
        end if ! centre or face excitation
           	     
        excitation_volumes(excitation_volume)%cell_list(excitation_cell)=problem_volumes(volume_number)%cell_list(excitation_cell)
			 
      end do !next cell cell in this volume	  
  
    end do ! next excitation volume

  end if ! n_excitation_volumes.GT.0
  
  CALL write_line('FINISHED: set_excitation_volumes_in_mesh',0,output_to_screen_flag)

  RETURN

END SUBROUTINE set_excitation_volumes_in_mesh
!
! NAME
!     initialise_excitation_volumes
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 10/05/2017 CJS
!
!
SUBROUTINE initialise_excitation_volumes

USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE mesh
USE TLM_excitation
USE file_information
USE constants

IMPLICIT NONE

! local variables

  integer	:: excitation_volume
  integer	:: cx,cy,cz,face
  integer	:: number_of_cells
  integer	:: excitation_cell
    
! START

  
  CALL write_line('CALLED: initialise_excitation_volumes',0,output_to_screen_flag)

! EXCITATION volumeS

  if (n_excitation_volumes.GT.0) then

    if (rank.eq.0) then
      write(info_file_unit,*)'Number of Excitation volumes=',n_excitation_volumes
    end if

    do excitation_volume=1,n_excitation_volumes
    
      number_of_cells=excitation_volumes(excitation_volume)%number_of_cells
      
      if (rank.eq.0) then
        write(info_file_unit,*)'Excitation volume number',excitation_volume,' Number of cell_cells=',number_of_cells
      end if
      
      ALLOCATE( excitation_volumes(excitation_volume)%cell_excitation_field_number_list(1:number_of_cells) )
      
      do excitation_cell=1,number_of_cells
   
        cx  =excitation_volumes(excitation_volume)%cell_list(excitation_cell)%cell%i
        cy  =excitation_volumes(excitation_volume)%cell_list(excitation_cell)%cell%j
        cz  =excitation_volumes(excitation_volume)%cell_list(excitation_cell)%cell%k
        face=excitation_volumes(excitation_volume)%cell_list(excitation_cell)%point 
      
        if (face.eq.centre) then
! cell centre excitation
  
          excitation_volumes(excitation_volume)%cell_excitation_field_number_list(excitation_cell)=local_cell_excitation(cx,cy,cz)
  
        else
  
          write(*,*)'ERROR, volume excitation should be at the cell centre'
          STOP 
 
        end if ! centre or face excitation
		  
      end do !next cell in this volume	  

    end do ! next excitation volume

  end if ! n_excitation_volumes.GT.0 
  
  CALL write_line('FINISHED: initialise_excitation_volumes',0,output_to_screen_flag)

  RETURN

END SUBROUTINE initialise_excitation_volumes
!
! SUBROUTINE volume_excitation
!
! NAME
!     volume_excitation
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 10/05/2017 CJS
!
SUBROUTINE volume_excitation

USE TLM_general
USE mesh
USE TLM_excitation
USE constants

IMPLICIT NONE

! local variables

  integer 	:: excitation_volume,excitation_cell
  integer 	:: cz,face
  integer 	:: number_of_cells
  integer 	:: function_number
  integer 	:: excitation_array_point
  integer 	:: field_component
  real*8	:: offset,offset_min
  real*8 	:: value  

! START
  
  CALL write_line('CALLED: volume_excitation',0,timestepping_output_to_screen_flag)

  if (n_excitation_volumes.gt.0) then

! excitation volumeS
    do excitation_volume=1,n_excitation_volumes
    
      function_number=excitation_volumes(excitation_volume)%excitation_function_number
      field_component=excitation_volumes(excitation_volume)%field_component
      value=excitation_functions(function_number)%value(timestep)
      
      number_of_cells=excitation_volumes(excitation_volume)%number_of_cells
            
      do excitation_cell=1,number_of_cells
            
        cz    =excitation_volumes(excitation_volume)%cell_list(excitation_cell)%cell%k
        face  =excitation_volumes(excitation_volume)%cell_list(excitation_cell)%point
      
        if (rank.eq.cell_rank(cz)) then
      
          function_number=excitation_volumes(excitation_volume)%excitation_function_number
          field_component=excitation_volumes(excitation_volume)%field_component

	  excitation_array_point=excitation_volumes(excitation_volume)%cell_excitation_field_number_list(excitation_cell)
          value=excitation_functions(function_number)%value(timestep)
	  
          cell_excitation_field(excitation_array_point,field_component)=	&
	       cell_excitation_field(excitation_array_point,field_component)+value
	       
! Set hard or soft excitation type	       
          if (field_component.le.6) then
            cell_excitation_type(excitation_array_point,field_component)=excitation_volumes(excitation_volume)%source_type
	  end if
	  
	end if ! excitation point in this processor's mesh
		  
      end do !next cell cell in this volume	  

    end do ! next excitation volume

  end if ! n_excitation_volumes.GT.0
  
  CALL write_line('FINISHED: volume_excitation',0,timestepping_output_to_screen_flag)
  
  RETURN

END SUBROUTINE volume_excitation
