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
!SUBROUTINE initialise_mode_stir_surfaces
!SUBROUTINE mode_stir_surfaces
!FUNCTION random_int(imin,imax)
!
! NAME
!     initialise_mode_stir_surfaces
!
! DESCRIPTION
!     
!     
! COMMENTS
!     We really need to only include faces which have PEC properties - a surface may be partly overwritten
!     by another surface with different properties.
!
! HISTORY
!
!     started 14/01/2013 CJS
!
!
SUBROUTINE initialise_mode_stir_surfaces

USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE mesh
USE mode_stir
USE TLM_surface_materials
USE file_information

IMPLICIT NONE

! local variables

  integer	:: mode_stir_surface
  integer	:: geometric_surface
  integer	:: surface_number
  integer	:: material_number
  integer	:: n_faces,n_mode_stir_faces
  integer	:: n_ports,port
  integer	:: face_loop
  integer	:: orientation
  
  integer	:: material_type
  
  integer	:: cx,cy,cz,face
  
  type(cell_point)	:: face1
  type(cell_point)	:: face2
  
  integer	:: swap,port1,port2
  
  integer	:: rank0_n_ports
  
  integer	:: n_data,talk_to
  
  integer	:: total_number_of_mode_stir_faces
  integer	:: total_number_of_mode_stir_ports
  
  integer	:: voltage_sign
  
! function types

  integer	::  random_int
  
! START

  
  CALL write_line('CALLED: initialise_mode_stir_surfaces',0,output_to_screen_flag)
   
  if (n_mode_stir_surfaces.ne.0) then

    do mode_stir_surface=1,n_mode_stir_surfaces

! count the number of TLM cell faces in this mode stirred surface          
      n_mode_stir_faces=0
 
      do geometric_surface=1,mode_stir_surface_list(mode_stir_surface)%number_of_surfaces
      
        surface_number= mode_stir_surface_list(mode_stir_surface)%surface_list(geometric_surface)
	
	if (surface_number.EQ.0) then
! outer boundary surface
! count the number of faces contributing to this mode stirred boundary from this process
	  
	  if (rank.eq.0) then ! include zmin face 
	    n_mode_stir_faces=n_mode_stir_faces+nx*ny
	  end if
	  
	  if (rank.eq.np-1) then ! include zmax face 
	    n_mode_stir_faces=n_mode_stir_faces+nx*ny
	  end if
	  
          n_mode_stir_faces=n_mode_stir_faces+2*(nx*(nz2-nz1+1))+2*(ny*(nz2-nz1+1)) ! ymin, ymax, xmin, xmax faces
	  
        else
! normal geometric surface

! Check that the surface has been given PEC properties - this is a requirement for the 
! mode stirred boundary to work correctly.

          material_number=problem_surfaces(surface_number)%surface_material_number
	  material_type=surface_material_list(material_number)%type
	  if (material_type.NE.surface_material_type_PEC) then
	    write(*,*)'Error in initialise_mode_stir_surfaces'
	    write(*,*)'Mode stirred boundary geometry surfaces must be given PEC properties'
	    write(*,*)'Mode_stir_surface number',mode_stir_surface
	    write(*,*)'Geometry surface number',surface_number
	    STOP
	  end if

!loop over faces checking individual face material
	  
	  n_faces=0
	  do face_loop=1,problem_surfaces(surface_number)%number_of_faces
	  
	    CALL get_min_face(problem_surfaces(surface_number)%face_list(face_loop),face2)
	    cx=face2%cell%i
	    cy=face2%cell%j
	    cz=face2%cell%k
	    face=face2%point
	    
	    material_number=abs( local_surface_material(cx  ,cy  ,cz  ,face) )
	    material_type=surface_material_list(material_number)%type
	    if (material_type.EQ.surface_material_type_PEC) then
	      n_faces=n_faces+1
	    end if
	    
	  end do
	  
	  n_mode_stir_faces=n_mode_stir_faces+n_faces

        end if ! surface number.ne.0

      end do ! next geometric surface
      
      n_ports=n_mode_stir_faces*2
      
      mode_stir_surface_list(mode_stir_surface)%number_of_faces=n_mode_stir_faces
      mode_stir_surface_list(mode_stir_surface)%number_of_ports=n_ports
      
#if defined(MPI)

      if (np.gt.1) then
! send the total number of faces and number of port voltages to the rank 0 process
! The rank 0 process must be capable of storing all the mode stir voltage pulses for the surface
        n_data=1
      
        CALL MPI_ALLREDUCE( n_mode_stir_faces,mode_stir_surface_list(mode_stir_surface)%total_number_of_mode_stir_faces	&
                           ,n_data,MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)
        CALL MPI_ALLREDUCE( n_ports,mode_stir_surface_list(mode_stir_surface)%total_number_of_mode_stir_ports	&
                           ,n_data,MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierror)
			   
      end if

#endif
      
      if (np.eq.1) then
      
        mode_stir_surface_list(mode_stir_surface)%total_number_of_mode_stir_faces=n_mode_stir_faces
        mode_stir_surface_list(mode_stir_surface)%total_number_of_mode_stir_ports=n_ports
	
      end if

! allocate memory for the face list and the port voltage list
	
      ALLOCATE( mode_stir_surface_list(mode_stir_surface)%face_list(1:n_mode_stir_faces) )
      ALLOCATE( mode_stir_surface_list(mode_stir_surface)%port1(1:n_mode_stir_faces) )
      ALLOCATE( mode_stir_surface_list(mode_stir_surface)%port2(1:n_mode_stir_faces) )

      if (rank.ne.0) then
      
        ALLOCATE( mode_stir_surface_list(mode_stir_surface)%port_voltage_list(1:n_ports) )
        mode_stir_surface_list(mode_stir_surface)%port_voltage_list(1:n_ports)=0d0
	
      else
      
        rank0_n_ports=mode_stir_surface_list(mode_stir_surface)%total_number_of_mode_stir_ports
        ALLOCATE( mode_stir_surface_list(mode_stir_surface)%port_voltage_list(1:rank0_n_ports) )
        mode_stir_surface_list(mode_stir_surface)%port_voltage_list(1:rank0_n_ports)=0d0
      
      end if

! Fill the face list         
      n_mode_stir_faces=0
    
      do geometric_surface=1,mode_stir_surface_list(mode_stir_surface)%number_of_surfaces
      
        surface_number= mode_stir_surface_list(mode_stir_surface)%surface_list(geometric_surface)
	
	if (surface_number.EQ.0) then
! outer boundary surface
! count the number of faces contributing to this mode stirred boundary from this process
	  
	  if (rank.eq.0) then ! include zmin face 
	  
	    cz=1
	    face=face_zmin
	    do cx=1,nx
	      do cy=1,ny
	        n_mode_stir_faces=n_mode_stir_faces+1
		mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%i=cx
		mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%j=cy
		mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%k=cz
		mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%point=face
		mode_stir_surface_list(mode_stir_surface)%port1(n_mode_stir_faces)=Vx_zmin
		mode_stir_surface_list(mode_stir_surface)%port2(n_mode_stir_faces)=Vy_zmin
	      end do
	    end do
	    
	  end if
	  
	  if (rank.eq.np-1) then ! include zmax face 
	  
	    cz=nz
	    face=face_zmax
	    do cx=1,nx
	      do cy=1,ny
	        n_mode_stir_faces=n_mode_stir_faces+1
		mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%i=cx
		mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%j=cy
		mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%k=cz
		mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%point=face
		mode_stir_surface_list(mode_stir_surface)%port1(n_mode_stir_faces)=Vx_zmax
		mode_stir_surface_list(mode_stir_surface)%port2(n_mode_stir_faces)=Vy_zmax
	      end do
	    end do

	  end if

! xmin face	  
	  cx=1
	  face=face_xmin
	  do cz=nz1,nz2
	    do cy=1,ny
	      n_mode_stir_faces=n_mode_stir_faces+1
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%i=cx
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%j=cy
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%k=cz
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%point=face
	      mode_stir_surface_list(mode_stir_surface)%port1(n_mode_stir_faces)=Vy_xmin
	      mode_stir_surface_list(mode_stir_surface)%port2(n_mode_stir_faces)=Vz_xmin
	    end do
	  end do

! xmax face	  
	  cx=nx
	  face=face_xmax
	  do cz=nz1,nz2
	    do cy=1,ny
	      n_mode_stir_faces=n_mode_stir_faces+1
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%i=cx
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%j=cy
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%k=cz
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%point=face
	      mode_stir_surface_list(mode_stir_surface)%port1(n_mode_stir_faces)=Vy_xmax
	      mode_stir_surface_list(mode_stir_surface)%port2(n_mode_stir_faces)=Vz_xmax
	    end do
	  end do

! ymin face	  
	  cy=1
	  face=face_ymin
	  do cz=nz1,nz2
	    do cx=1,nx
	      n_mode_stir_faces=n_mode_stir_faces+1
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%i=cx
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%j=cy
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%k=cz
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%point=face
	      mode_stir_surface_list(mode_stir_surface)%port1(n_mode_stir_faces)=Vx_ymin
	      mode_stir_surface_list(mode_stir_surface)%port2(n_mode_stir_faces)=Vz_ymin
	    end do
	  end do

! ymax face	  
	  cy=ny
	  face=face_ymax
	  do cz=nz1,nz2
	    do cx=1,nx
	      n_mode_stir_faces=n_mode_stir_faces+1
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%i=cx
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%j=cy
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%cell%k=cz
	      mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%point=face
	      mode_stir_surface_list(mode_stir_surface)%port1(n_mode_stir_faces)=Vx_ymax
	      mode_stir_surface_list(mode_stir_surface)%port2(n_mode_stir_faces)=Vz_ymax
	    end do
	  end do
	  
        else
! normal geometric surface

          orientation=mode_stir_surface_list(mode_stir_surface)%surface_orientation_list(geometric_surface)
	  
	  do face_loop=1,problem_surfaces(surface_number)%number_of_faces
	  
	  
	    CALL get_min_face(problem_surfaces(surface_number)%face_list(face_loop),face2)
	    cx=face2%cell%i
	    cy=face2%cell%j
	    cz=face2%cell%k
	    face=face2%point
	    
	    material_number=abs( local_surface_material(cx  ,cy  ,cz  ,face) )
	    material_type=surface_material_list(material_number)%type
	    if (material_type.EQ.surface_material_type_PEC) then
	    
	      n_mode_stir_faces=n_mode_stir_faces+1
	    
	      face1=problem_surfaces(surface_number)%face_list(face_loop)
	    
	      if (orientation.eq.1) then
	    
	        mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)=face1
	      
	      else
 ! get the opposite face	    
	        CALL get_other_side_of_face(face1,face2)
	        mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)=face2
	      
	      end if
	    
	      face=mode_stir_surface_list(mode_stir_surface)%face_list(n_mode_stir_faces)%point
	    
	      if      (face.eq.face_xmin) then
	        mode_stir_surface_list(mode_stir_surface)%port1(n_mode_stir_faces)=Vy_xmin
	        mode_stir_surface_list(mode_stir_surface)%port2(n_mode_stir_faces)=Vz_xmin
	      else if (face.eq.face_xmax) then
	        mode_stir_surface_list(mode_stir_surface)%port1(n_mode_stir_faces)=Vy_xmax
	        mode_stir_surface_list(mode_stir_surface)%port2(n_mode_stir_faces)=Vz_xmax
	      else if (face.eq.face_ymin) then
	        mode_stir_surface_list(mode_stir_surface)%port1(n_mode_stir_faces)=Vx_ymin
	        mode_stir_surface_list(mode_stir_surface)%port2(n_mode_stir_faces)=Vz_ymin
	      else if (face.eq.face_ymax) then
	        mode_stir_surface_list(mode_stir_surface)%port1(n_mode_stir_faces)=Vx_ymax
	        mode_stir_surface_list(mode_stir_surface)%port2(n_mode_stir_faces)=Vz_ymax
	      else if (face.eq.face_zmin) then
	        mode_stir_surface_list(mode_stir_surface)%port1(n_mode_stir_faces)=Vx_zmin
	        mode_stir_surface_list(mode_stir_surface)%port2(n_mode_stir_faces)=Vy_zmin
	      else if (face.eq.face_zmax) then
	        mode_stir_surface_list(mode_stir_surface)%port1(n_mode_stir_faces)=Vx_zmax
	        mode_stir_surface_list(mode_stir_surface)%port2(n_mode_stir_faces)=Vy_zmax
	      end if  
	    
	    end if ! include this face in the list
	    
	  end do ! next face in this geometric surface

        end if ! surface number.ne.0

      end do ! next geometric surface

! check
      if (n_mode_stir_faces.NE.mode_stir_surface_list(mode_stir_surface)%number_of_faces) then
        write(*,*)'Error in initialise_mode_stir_surfaces'
	write(*,*)'face counting discrepancy'
	write(*,*)n_mode_stir_faces,mode_stir_surface_list(mode_stir_surface)%number_of_faces
	STOP
      end if
            
! tell the rank 0 process the number of mode_stir ports in each process

      if (rank.eq.0) then
        ALLOCATE( mode_stir_surface_list(mode_stir_surface)%number_of_ports_rank(0:np-1) )
	mode_stir_surface_list(mode_stir_surface)%number_of_ports_rank(0)=	&
	     mode_stir_surface_list(mode_stir_surface)%number_of_ports      
      end if
      
#if defined(MPI)

      n_data=1
      if (rank.ne.0) then      
        talk_to=0
	n_ports=mode_stir_surface_list(mode_stir_surface)%number_of_ports
        CALL MPI_SEND(n_ports,n_data,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)
      else
! rank.eq.0
	do talk_to=1,np-1
	  CALL MPI_RECV(n_ports,n_data,MPI_INTEGER,talk_to,0,MPI_COMM_WORLD,status,ierror)
	  mode_stir_surface_list(mode_stir_surface)%number_of_ports_rank(talk_to)=n_ports
	end do
      end if ! rank

#endif

! set up the port pairing list in the rank 0 process
  
      if (rank.eq.0) then
  
! create port pairing list
        total_number_of_mode_stir_ports=mode_stir_surface_list(mode_stir_surface)%total_number_of_mode_stir_ports
	
        ALLOCATE( mode_stir_surface_list(mode_stir_surface)%mode_stir_voltage_pairing(1:total_number_of_mode_stir_ports) )
    
        do port=1,total_number_of_mode_stir_ports
          mode_stir_surface_list(mode_stir_surface)%mode_stir_voltage_pairing(port)=port
        end do
    
! randomise the pairing array
        do port=1,total_number_of_mode_stir_ports
    
          port1=port
          port2=random_int(1,total_number_of_mode_stir_ports)
          swap=mode_stir_surface_list(mode_stir_surface)%mode_stir_voltage_pairing(port1)
      
          mode_stir_surface_list(mode_stir_surface)%mode_stir_voltage_pairing(port1)=	&
               mode_stir_surface_list(mode_stir_surface)%mode_stir_voltage_pairing(port2)
	   
          mode_stir_surface_list(mode_stir_surface)%mode_stir_voltage_pairing(port2)=swap
      
        end do
	    
! randomise the signs of the connection ports
    
        total_number_of_mode_stir_faces=mode_stir_surface_list(mode_stir_surface)%total_number_of_mode_stir_faces
	
        ALLOCATE( mode_stir_surface_list(mode_stir_surface)%sign(1:total_number_of_mode_stir_faces) )

        do face_loop=1,total_number_of_mode_stir_faces
    
	  voltage_sign=random_int(1,2)
	  if (voltage_sign.eq.2) voltage_sign=-1
	  mode_stir_surface_list(mode_stir_surface)%sign(face_loop)=voltage_sign
       
        end do
	
	write(info_file_unit,*)'Mode stir surface number:',mode_stir_surface
	write(info_file_unit,*)'Number of faces=',total_number_of_mode_stir_faces
	write(info_file_unit,*)'Number of ports=',total_number_of_mode_stir_ports
	
      end if ! rank.eq.0
	
    end do ! next mode stir surface
          
  end if ! n_mode_stir_surfaces.ne.0
  
  CALL write_line('FINISHED: initialise_mode_stir_surfaces',0,output_to_screen_flag)

  RETURN

END SUBROUTINE initialise_mode_stir_surfaces
!
! NAME
!     mode_stir_surfaces
!
! DESCRIPTION
!     Implement the mode stirring boundary condition
!     
! COMMENTS
!     Called immediately before the connect process
!
! HISTORY
!
!     started 14/1/2013 CJS
!     
!
!
SUBROUTINE mode_stir_surfaces

USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE mesh
USE mode_stir

IMPLICIT NONE

! local variables

  integer	:: mode_stir_surface
  integer	:: face_loop
  
  integer	:: cx,cy,cz,face
  
  integer	:: port_count
  integer	:: port_pair
  integer	:: first_port,last_port
    
  integer	:: port1,port2,sign

  integer	:: n_data
  integer	:: talk_to
  
  real*8	:: Vswap
  real*8	:: Rm
  real*8	:: amplitude
  real*8	:: r_random
  
! START
  
  CALL write_line('CALLED: mode_stir_surfaces',0,timestepping_output_to_screen_flag)
   
  if (n_mode_stir_surfaces.ne.0) then

    do mode_stir_surface=1,n_mode_stir_surfaces

! Copy the mode_stir surface port voltages into the port_voltage_list
      port_count=0

      do face_loop=1,mode_stir_surface_list(mode_stir_surface)%number_of_faces
      
        cx=mode_stir_surface_list(mode_stir_surface)%face_list(face_loop)%cell%i
	cy=mode_stir_surface_list(mode_stir_surface)%face_list(face_loop)%cell%j
	cz=mode_stir_surface_list(mode_stir_surface)%face_list(face_loop)%cell%k
	face=mode_stir_surface_list(mode_stir_surface)%face_list(face_loop)%point
	port1=mode_stir_surface_list(mode_stir_surface)%port1(face_loop)
	port2=mode_stir_surface_list(mode_stir_surface)%port2(face_loop)
        port_count=port_count+1
	mode_stir_surface_list(mode_stir_surface)%port_voltage_list(port_count)=V(port1,cx,cy,cz)
        port_count=port_count+1
	mode_stir_surface_list(mode_stir_surface)%port_voltage_list(port_count)=V(port2,cx,cy,cz)
	
      end do ! next face

! Send all the port_voltage_lists to the rank 0 process

#if defined(MPI)      

      if (rank.ne.0) then      
      
        talk_to=0
	n_data=mode_stir_surface_list(mode_stir_surface)%number_of_ports
	if (n_data.ne.0) then
          CALL MPI_SEND(mode_stir_surface_list(mode_stir_surface)%port_voltage_list(1:n_data),n_data,	&
	                MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
	end if
	
      else
! rank.eq.0
        port_count=mode_stir_surface_list(mode_stir_surface)%number_of_ports_rank(0)
        do talk_to=1,np-1
          n_data=mode_stir_surface_list(mode_stir_surface)%number_of_ports_rank(talk_to)
	  if (n_data.ne.0) then
	    first_port=port_count+1
	    last_port=first_port+n_data-1
            CALL MPI_RECV(mode_stir_surface_list(mode_stir_surface)%port_voltage_list(first_port:last_port),n_data,&
	                  MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
	    port_count=last_port
	  end if
        end do
	
      end if ! rank

#endif
      
      if (rank.eq.0) then
      
! Do the port swapping process

        Rm=mode_stir_surface_list(mode_stir_surface)%R_ms

        do port_pair=1,mode_stir_surface_list(mode_stir_surface)%total_number_of_mode_stir_faces
	  port_count=port_pair*2-1
          port1=mode_stir_surface_list(mode_stir_surface)%mode_stir_voltage_pairing(port_count)
	  port2=mode_stir_surface_list(mode_stir_surface)%mode_stir_voltage_pairing(port_count+1)
	  sign=mode_stir_surface_list(mode_stir_surface)%sign(port_pair)
	  
	  Vswap=sign*Rm*mode_stir_surface_list(mode_stir_surface)%port_voltage_list(port1)
	  mode_stir_surface_list(mode_stir_surface)%port_voltage_list(port1)=	&
	        sign*Rm*mode_stir_surface_list(mode_stir_surface)%port_voltage_list(port2)
	  mode_stir_surface_list(mode_stir_surface)%port_voltage_list(port2)=Vswap
	  
        end do
	  
	if (timestep.eq.1) then
! add the impulse excitation

          amplitude=mode_stir_surface_list(mode_stir_surface)%impulse_amplitude
	  
          do port_pair=1,mode_stir_surface_list(mode_stir_surface)%total_number_of_mode_stir_faces
	    port_count=port_pair*2-1
            port1=mode_stir_surface_list(mode_stir_surface)%mode_stir_voltage_pairing(port_count)
	    port2=mode_stir_surface_list(mode_stir_surface)%mode_stir_voltage_pairing(port_count+1)
            
	    CALL random_number(r_random)
	    mode_stir_surface_list(mode_stir_surface)%port_voltage_list(port1)=	&
	      mode_stir_surface_list(mode_stir_surface)%port_voltage_list(port1)+(2d0*r_random-1d0)*amplitude
	    
	    CALL random_number(r_random)
	    mode_stir_surface_list(mode_stir_surface)%port_voltage_list(port2)= &
	      mode_stir_surface_list(mode_stir_surface)%port_voltage_list(port2)+(2d0*r_random-1d0)*amplitude
	  
          end do

	end if
	
      end if

! Send the rank 0 mode stirred port_voltage_lists back to their own processes

#if defined(MPI)      

      if (rank.eq.0) then   
         
        port_count=mode_stir_surface_list(mode_stir_surface)%number_of_ports_rank(0)
        do talk_to=1,np-1
          n_data=mode_stir_surface_list(mode_stir_surface)%number_of_ports_rank(talk_to)
	  if (n_data.ne.0) then
	    first_port=port_count+1
	    last_port=first_port+n_data-1
            CALL MPI_SEND(mode_stir_surface_list(mode_stir_surface)%port_voltage_list(first_port:last_port),n_data,&
	                  MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
	    port_count=last_port
	  end if
        end do
	
      else ! rank.ne.0
      
        talk_to=0
	n_data=mode_stir_surface_list(mode_stir_surface)%number_of_ports
	if (n_data.ne.0) then
          CALL MPI_RECV(mode_stir_surface_list(mode_stir_surface)%port_voltage_list(1:n_data),n_data,	&
	                MPI_DOUBLE_PRECISION,talk_to,0,MPI_COMM_WORLD,status,ierror)
	end if
	
      end if ! rank

#endif
      
! Copy the mode_stirred port_voltage_list data back to the surface port voltages
      port_count=0

      do face_loop=1,mode_stir_surface_list(mode_stir_surface)%number_of_faces
      
        cx=mode_stir_surface_list(mode_stir_surface)%face_list(face_loop)%cell%i
	cy=mode_stir_surface_list(mode_stir_surface)%face_list(face_loop)%cell%j
	cz=mode_stir_surface_list(mode_stir_surface)%face_list(face_loop)%cell%k
	face=mode_stir_surface_list(mode_stir_surface)%face_list(face_loop)%point
	port1=mode_stir_surface_list(mode_stir_surface)%port1(face_loop)
	port2=mode_stir_surface_list(mode_stir_surface)%port2(face_loop)
        port_count=port_count+1
	V(port1,cx,cy,cz)=mode_stir_surface_list(mode_stir_surface)%port_voltage_list(port_count)
        port_count=port_count+1
	V(port2,cx,cy,cz)=mode_stir_surface_list(mode_stir_surface)%port_voltage_list(port_count)

      end do ! next face

    end do ! next mode stir surface
        
  end if ! n_mode_stir_surfaces.ne.0
  
  
  CALL write_line('FINISHED: mode_stir_surfaces',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE mode_stir_surfaces
!
! NAME
!     return a random integer in the given range
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 15/01/2013 CJS
!
!
FUNCTION random_int(imin,imax)

  integer random_int
  integer imin,imax
  
  real*8 r_random,range
  
! START

  CALL random_number(r_random)
  
  range=imax-imin+1
  
  random_int=int(r_random*range)+imin
  
  RETURN
  
END FUNCTION random_int
