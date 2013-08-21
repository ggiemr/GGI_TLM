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
! SUBROUTINE set_SAR_volumes_in_mesh
! SUBROUTINE initialise_SAR_volumes
! SUBROUTINE cell_output_SAR_volumes
! SUBROUTINE write_SAR_volumes
!
! NAME
!     set_SAR_volumes_in_mesh
!
! DESCRIPTION
!
!     
! COMMENTS
!
! HISTORY
!
!     started 9/1/2013 CJS
!
!
SUBROUTINE set_SAR_volumes_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE file_information

IMPLICIT NONE

! local variables

  integer		:: output_volume
  integer		:: material_number
  integer		:: cx,cy,cz
  integer		:: cell_count

! START
  
  CALL write_line('CALLED: set_SAR_volumes_in_mesh',0,output_to_screen_flag)
    
  if (n_SAR_volumes.gt.0) then

    do output_volume=1,n_SAR_volumes
    
      material_number=SAR_volume_list(output_volume)%material_number
      
      SAR_volume_list(output_volume)%number_of_cells=0

! loop over this processor's mesh looking for the material number      
      do cz=nz1,nz2
        do cy=1,ny
          do cx=1,nx
	  
! check the material of this cell	  
	    if ( local_cell_material(cx,cy,cz).EQ.material_number) then
	      SAR_volume_list(output_volume)%number_of_cells=SAR_volume_list(output_volume)%number_of_cells+1
              local_cell_output(cx,cy,cz)=1
	    end if
	    
          end do ! next cx
	end do ! next cy
      end do ! next cz
      
      if (SAR_volume_list(output_volume)%number_of_cells.GT.0) then
        ALLOCATE( SAR_volume_list(output_volume)%cell_list(1:SAR_volume_list(output_volume)%number_of_cells) )
      end if

! loop over this processor's mesh allocating the cell_list   
      cell_count=0
      
      do cz=nz1,nz2
        do cy=1,ny
          do cx=1,nx
	  
! check the material of this cell	  
	    if ( local_cell_material(cx,cy,cz).EQ.material_number) then
	      cell_count=cell_count+1
	      SAR_volume_list(output_volume)%cell_list(cell_count)%i=cx
	      SAR_volume_list(output_volume)%cell_list(cell_count)%j=cy
	      SAR_volume_list(output_volume)%cell_list(cell_count)%k=cz
	    end if
	    
          end do ! next cx
	end do ! next cy
      end do ! next cz
        
    end do ! next SAR volume

  end if ! n_SAR_volumes.GT.0 
  
  CALL write_line('FINISHED: set_SAR_volumes_in_mesh',0,output_to_screen_flag)
  
END SUBROUTINE set_SAR_volumes_in_mesh
!
! NAME
!     initialise_SAR_volumes
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 9/1/2013 CJS
!
!
SUBROUTINE initialise_SAR_volumes

USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE TLM_volume_materials
USE file_information

IMPLICIT NONE

! local variables

  integer	:: output_volume
  integer	:: material_number
  integer	:: material_order
  integer	:: cx,cy,cz
  integer	:: cell_count,number_of_cells

! START
  
  if (n_SAR_volumes.gt.0) then

    do output_volume=1,n_SAR_volumes
    
      material_number=SAR_volume_list(output_volume)%material_number

! check that we have a simple material characterised by 0 order dispersion and electric conductivity only      
      if (volume_material_list(material_number)%type.EQ.	&
             volume_material_type_PEC) then
        write(*,*)'Cant calculate SAR in PEC'
        STOP
      else if (volume_material_list(material_number)%type.EQ.	&
               volume_material_type_PMC) then
        write(*,*)'Cant calculate SAR in PEC'
        STOP
      else

! check that we don't have a dispersive dielecric material      
        material_order=volume_material_list(material_number)%eps_S%a%order
        if (material_order.NE.0) then
          write(*,*)'Cant calculate SAR in dispersive dielectric materials'
          STOP
	end if
        material_order=volume_material_list(material_number)%eps_S%b%order
        if (material_order.NE.0) then
          write(*,*)'Cant calculate SAR in dispersive dielectric materials'
          STOP
	end if

! check that we don't have a dispersive magnetic material      
        material_order=volume_material_list(material_number)%mu_S%a%order
        if (material_order.NE.0) then
          write(*,*)'Cant calculate SAR in dispersive magnetic materials'
          STOP
	end if
        material_order=volume_material_list(material_number)%mu_S%b%order
        if (material_order.NE.0) then
          write(*,*)'Cant calculate SAR in dispersive magnetic materials'
          STOP
	end if
	
        if (volume_material_list(material_number)%sigma_m.NE.0d0) then
          write(*,*)'Cant calculate SAR in materials with magnetic conductivity'
          STOP
	end if
      
      end if
            
      SAR_volume_list(output_volume)%conductivity=volume_material_list(material_number)%sigma_e
      SAR_volume_list(output_volume)%cell_volume=dl*dl*dl
      SAR_volume_list(output_volume)%mass=SAR_volume_list(output_volume)%cell_volume*	&
                                          SAR_volume_list(output_volume)%density*	&
					  SAR_volume_list(output_volume)%number_of_cells 
      SAR_volume_list(output_volume)%SAR=0d0
       
      number_of_cells=SAR_volume_list(output_volume)%number_of_cells

      if (number_of_cells.GT.0) then
        ALLOCATE( SAR_volume_list(output_volume)%cell_output_field_number_list(1:number_of_cells) )
        ALLOCATE( SAR_volume_list(output_volume)%Ex(1:number_of_cells) )
        ALLOCATE( SAR_volume_list(output_volume)%Ey(1:number_of_cells) )
        ALLOCATE( SAR_volume_list(output_volume)%Ez(1:number_of_cells) )
	SAR_volume_list(output_volume)%Ex(1:number_of_cells)=(0d0,0d0)
	SAR_volume_list(output_volume)%Ey(1:number_of_cells)=(0d0,0d0)
	SAR_volume_list(output_volume)%Ez(1:number_of_cells)=(0d0,0d0)
      end if

! loop over this processor's cell_list   
      do cell_count=1,number_of_cells
      
	cx=SAR_volume_list(output_volume)%cell_list(cell_count)%i
	cy=SAR_volume_list(output_volume)%cell_list(cell_count)%j
	cz=SAR_volume_list(output_volume)%cell_list(cell_count)%k
	
	SAR_volume_list(output_volume)%cell_output_field_number_list(cell_count)=local_cell_output(cx,cy,cz)

      end do ! next cell in the cell list
        
    end do ! next SAR volume

    if (rank.eq.0) then
! rank 0 process only: write header for SAR outputs
  
      OPEN(unit=SAR_output_unit,file=trim(problem_name)//SAR_output_extn)

    end if ! rank 0 process

  end if ! n_SAR_volumes.GT.0   

  RETURN

END SUBROUTINE initialise_SAR_volumes

!
! NAME
!     cell_output_SAR_volumes
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 9/1/2013 CJS
!
!
SUBROUTINE cell_output_SAR_volumes

USE TLM_general
USE mesh
USE TLM_output
USE file_information
USE output_formats
USE cell_parameters
USE constants

IMPLICIT NONE

! local variables

  integer	:: output_volume
  integer	:: material_number
  integer	:: cell_count,number_of_cells
  
  integer	:: output_field_number
  
  real*8	:: Ex_value,Ey_value,Ez_value
  
  real*8	:: frequency
  complex*16	:: ejwt

! START
  
  CALL write_line('CALLED: cell_output_SAR_volumes',0,timestepping_output_to_screen_flag)
  
  if (n_SAR_volumes.gt.0) then

    do output_volume=1,n_SAR_volumes
    
      frequency=SAR_volume_list(output_volume)%frequency
      ejwt=exp(-j*2d0*pi*frequency*time)
       
      number_of_cells=SAR_volume_list(output_volume)%number_of_cells

! loop over this processor's cell_list   
      do cell_count=1,number_of_cells
	
	output_field_number=SAR_volume_list(output_volume)%cell_output_field_number_list(cell_count)

        Ex_value=cell_output_field(output_field_number,Ex)  
	Ey_value=cell_output_field(output_field_number,Ey)  
	Ez_value=cell_output_field(output_field_number,Ez) 
	 
	SAR_volume_list(output_volume)%Ex(cell_count)=SAR_volume_list(output_volume)%Ex(cell_count)+ejwt*Ex_value*dt
	SAR_volume_list(output_volume)%Ey(cell_count)=SAR_volume_list(output_volume)%Ey(cell_count)+ejwt*Ey_value*dt
	SAR_volume_list(output_volume)%Ez(cell_count)=SAR_volume_list(output_volume)%Ez(cell_count)+ejwt*Ez_value*dt

      end do ! next cell in the cell list
        
    end do ! next SAR volume

  end if ! n_SAR_volumes.GT.0   
  
  CALL write_line('FINISHED: cell_output_SAR_volumes',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE cell_output_SAR_volumes

! Name write_SAR_volumes
!     
!
! Description
!     
!
! Comments:
!      
!
! History
!
!     started 9/1/2013 CJS 
!

SUBROUTINE write_SAR_volumes

USE file_information
USE TLM_output
USE mesh
Use TLM_general
Use constants

IMPLICIT NONE

! variables passed to subroutine

! local_variables

  integer :: output_volume
  
  integer :: cx,cy,cz
  
  complex*16 :: Ex_value
  complex*16 :: Ey_value
  complex*16 :: Ez_value
 
  integer,allocatable	:: n_cells_rank(:)
  
  integer	:: talk_to
  integer	:: n_integer
  integer	:: loop
  integer 	:: cell_count
  integer 	:: cell_count_save
 
  integer	:: local_n_cells
  integer 	:: output_cell
  type(SAR_volume_type)	:: local_SAR_volume
  
  integer,allocatable	:: integer_array(:)
  complex*16,allocatable:: complex_array(:)

! START

  write(*,*)'CALLED: write_SAR_volumes'

! write SAR volume

  if (n_SAR_volumes.ne.0) then

    do output_volume=1,n_SAR_volumes
    
! Stage 1: send the number of SAR_volume cells in each process to rank 0

      if (rank.eq.0) then
  
        ALLOCATE ( n_cells_rank(0:np-1) )
    
        n_cells_rank(0)=SAR_volume_list(output_volume)%number_of_cells
      
      end if
      
#if defined(MPI)

      if (rank.ne.0) then 
    
        talk_to=0
        n_integer=1
        CALL MPI_SEND(SAR_volume_list(output_volume)%number_of_cells,n_integer,	&
	              MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    

      else if (rank.eq.0) then 
    
        do talk_to=1,np-1
          n_integer=1
          CALL MPI_RECV(n_cells_rank(talk_to),n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
        end do
    
      end if

#endif

! Stage 2: Count the total number of SAR output cells
      if (rank.eq.0) then
      
        local_n_cells=0
        do loop=0,np-1
          local_n_cells=local_n_cells+n_cells_rank(loop)
        end do
	
! Stage 3: Allocate memory for the whole output dataset in the rank 0 process

        local_SAR_volume%number_of_cells=local_n_cells
        ALLOCATE( local_SAR_volume%cell_list(1:local_n_cells) )
        ALLOCATE( local_SAR_volume%Ex(1:local_n_cells) )
        ALLOCATE( local_SAR_volume%Ey(1:local_n_cells) )
        ALLOCATE( local_SAR_volume%Ez(1:local_n_cells) )

      end if

! STAGE 4: communicate all of the data to the rank 0 process

      if (rank.eq.0) then 
  
! set values already in rank 0 process    
	
        local_SAR_volume%conductivity=SAR_volume_list(output_volume)%conductivity
	local_SAR_volume%cell_volume=SAR_volume_list(output_volume)%cell_volume
	local_SAR_volume%density=SAR_volume_list(output_volume)%density
	local_SAR_volume%frequency=SAR_volume_list(output_volume)%frequency
	  
        cell_count=0
        do loop=1,n_cells_rank(0)
	  
          cell_count=cell_count+1
	  
	  local_SAR_volume%cell_list(cell_count)%i=	&
	    SAR_volume_list(output_volume)%cell_list(cell_count)%i
	  local_SAR_volume%cell_list(cell_count)%j=	&
	    SAR_volume_list(output_volume)%cell_list(cell_count)%j
	  local_SAR_volume%cell_list(cell_count)%k=	&
	    SAR_volume_list(output_volume)%cell_list(cell_count)%k
	  local_SAR_volume%Ex(cell_count)=	&
	    SAR_volume_list(output_volume)%Ex(cell_count)
	  local_SAR_volume%Ey(cell_count)=	&
	    SAR_volume_list(output_volume)%Ey(cell_count)
	  local_SAR_volume%Ez(cell_count)=	&
	    SAR_volume_list(output_volume)%Ez(cell_count)
        end do ! next rank 0 process

! calculate the total mass of all material cells	
        local_SAR_volume%mass=local_SAR_volume%cell_volume*	 &
                              local_SAR_volume%density* &
			      local_SAR_volume%number_of_cells 
      end if
      
#if defined(MPI)      

      if (rank.ne.0) then 

! send cx,cy,cz,face,value to rank 0 process    
        talk_to=0
        n_integer=SAR_volume_list(output_volume)%number_of_cells
	
	if (n_integer.ne.0) then
	
          ALLOCATE( integer_array(1:n_integer) )
          ALLOCATE( complex_array(1:n_integer) )
! send cx      
          integer_array(1:n_integer)=SAR_volume_list(output_volume)%cell_list(1:n_integer)%i
          CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send cy   
          integer_array(1:n_integer)=SAR_volume_list(output_volume)%cell_list(1:n_integer)%j
          CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send cz    
          integer_array(1:n_integer)=SAR_volume_list(output_volume)%cell_list(1:n_integer)%k
          CALL MPI_SEND(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send Ex value  
          complex_array(1:n_integer)=SAR_volume_list(output_volume)%Ex(1:n_integer)
          CALL MPI_SEND(complex_array,n_integer,MPI_DOUBLE_COMPLEX,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send Ey value  
          complex_array(1:n_integer)=SAR_volume_list(output_volume)%Ey(1:n_integer)
          CALL MPI_SEND(complex_array,n_integer,MPI_DOUBLE_COMPLEX,talk_to,0,MPI_COMM_WORLD,status,ierror)	    
! send Ez value  
          complex_array(1:n_integer)=SAR_volume_list(output_volume)%Ez(1:n_integer)
          CALL MPI_SEND(complex_array,n_integer,MPI_DOUBLE_COMPLEX,talk_to,0,MPI_COMM_WORLD,status,ierror)	    

          DEALLOCATE( integer_array )
          DEALLOCATE( complex_array )
	  
        end if ! number of SAR cells ne.0
	
      else if (rank.eq.0) then 
  	
        do talk_to=1,np-1
      
          cell_count_save=cell_count
	  
          n_integer=n_cells_rank(talk_to)
	  
	  if (n_integer.ne.0) then
	  
            ALLOCATE( integer_array(1:n_integer) )
            ALLOCATE( complex_array(1:n_integer) )
! get cx
            CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
            cell_count=cell_count_save
            do loop=1,n_integer
              cell_count=cell_count+1
	      local_SAR_volume%cell_list(cell_count)%i=integer_array(loop)
            end do
! get cy
            CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
            cell_count=cell_count_save
            do loop=1,n_integer
              cell_count=cell_count+1
	      local_SAR_volume%cell_list(cell_count)%j=integer_array(loop)
            end do
! get cz
            CALL MPI_RECV(integer_array,n_integer,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)	
            cell_count=cell_count_save
            do loop=1,n_integer
              cell_count=cell_count+1
	      local_SAR_volume%cell_list(cell_count)%k=integer_array(loop)
            end do
! get Ex
            CALL MPI_RECV(complex_array,n_integer,MPI_DOUBLE_COMPLEX,talk_to,0,MPI_COMM_WORLD,status,ierror)	
            cell_count=cell_count_save
            do loop=1,n_integer
              cell_count=cell_count+1
	      local_SAR_volume%Ex(cell_count)=complex_array(loop)
            end do
! get Ey
            CALL MPI_RECV(complex_array,n_integer,MPI_DOUBLE_COMPLEX,talk_to,0,MPI_COMM_WORLD,status,ierror)	
            cell_count=cell_count_save
            do loop=1,n_integer
              cell_count=cell_count+1
	      local_SAR_volume%Ey(cell_count)=complex_array(loop)
            end do
! get Ez
            CALL MPI_RECV(complex_array,n_integer,MPI_DOUBLE_COMPLEX,talk_to,0,MPI_COMM_WORLD,status,ierror)	
            cell_count=cell_count_save
            do loop=1,n_integer
              cell_count=cell_count+1
	      local_SAR_volume%Ez(cell_count)=complex_array(loop)
            end do

            DEALLOCATE( integer_array )
            DEALLOCATE( complex_array )
	
	  end if ! number of SAR cells ne.0
	
        end do ! next process to talk to

! check that we have the correct amount of data  

        if (cell_count.ne.local_n_cells) then
          write(*,*)'Error in write_SAR_volumes'
          write(*,*)'Data counting error',cell_count
          write(*,*)'Cell count',cell_count
          write(*,*)'local_n_cells',local_n_cells
          STOP
        end if
      
      end if ! rank 0 process

#endif

! STAGE 5: The rank 0 process writes the data to file

      if (rank.eq.0) then
      
        local_SAR_volume%SAR=0d0
    
        do output_cell=1,local_n_cells
         	
	  Ex_value=local_SAR_volume%Ex(output_cell)
	  Ey_value=local_SAR_volume%Ey(output_cell)
	  Ez_value=local_SAR_volume%Ez(output_cell)
	  
          local_SAR_volume%SAR=local_SAR_volume%SAR+0.5D0*( abs(Ex_value)**2+	&
			                                    abs(Ey_value)**2+	&
			                                    abs(Ez_value)**2 )
        
        end do !next cell   
	
!	write(*,*)'n_SAR_cells=',local_SAR_volume%number_of_cells
!	write(*,*)'conductivity=',local_SAR_volume%conductivity
!	write(*,*)'frequency=',local_SAR_volume%frequency
!	write(*,*)'density=',local_SAR_volume%density
!	write(*,*)'cell_volume=',local_SAR_volume%cell_volume
!	write(*,*)'mass=',local_SAR_volume%mass
     
        local_SAR_volume%SAR=local_SAR_volume%SAR*local_SAR_volume%conductivity*local_SAR_volume%cell_volume/	&
	                     local_SAR_volume%mass

        write(SAR_output_unit,*)local_SAR_volume%frequency,output_volume,local_SAR_volume%SAR
  
        DEALLOCATE ( n_cells_rank )
        DEALLOCATE ( local_SAR_volume%cell_list )
        DEALLOCATE ( local_SAR_volume%Ex )
        DEALLOCATE ( local_SAR_volume%Ey )
        DEALLOCATE ( local_SAR_volume%Ez )
   
      end if ! rank.eq.0
  
    end do ! next SAR_volume
  
  end if !(n_SAR_volumes.ne.0) 

  write(*,*)'FINISHED: write_SAR_volumes'
  return
  
END SUBROUTINE write_SAR_volumes

