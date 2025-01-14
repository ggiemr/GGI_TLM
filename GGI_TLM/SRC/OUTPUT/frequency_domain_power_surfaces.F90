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
! SUBROUTINE set_frequency_domain_power_surfaces_in_mesh
! SUBROUTINE initialise_frequency_domain_power_surface
! SUBROUTINE face_output_frequency_domain_power_surfaces
!
! NAME
!     set_frequency_domain_power_surfaces_in_mesh
!
! DESCRIPTION
!
!     
! COMMENTS
!
! HISTORY
!
!     started 8/2/2013 CJS
!
!
SUBROUTINE set_frequency_domain_power_surfaces_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE file_information

IMPLICIT NONE

! local variables

  integer		:: output_surface
  integer		:: surface_number
  integer		:: number_of_faces
  integer		:: output_face
  integer		:: cx,cy,cz,face
  type(cell_point)	:: output_face1
  type(cell_point)	:: output_face2

! START
  
  CALL write_line('CALLED: set_frequency_domain_power_surfaces_in_mesh',0,output_to_screen_flag)
    
  if (n_frequency_domain_power_surfaces.gt.0) then

! FREQUENCY OUTPUT SURFACES

    if (rank.eq.0) then
      write(info_file_unit,*)'Number of Frequency Domain Power Surfaces=',n_frequency_domain_power_surfaces
    end if
    
    do output_surface=1,n_frequency_domain_power_surfaces
    
      surface_number=frequency_domain_power_surface(output_surface)%surface_number
      number_of_faces=problem_surfaces(surface_number)%number_of_faces
      frequency_domain_power_surface(output_surface)%number_of_faces=number_of_faces
      
      if (rank.eq.0) then
        write(info_file_unit,*)'Frequency Domain Power Surface number',output_surface,' Number of cell_faces=',number_of_faces
      end if
      
! allocate face list for the output surface      
      ALLOCATE( frequency_domain_power_surface(output_surface)%face_list(1:number_of_faces) )
      
      do output_face=1,number_of_faces
      
! copy the faces from the geometry structure to the local
! frequency_domain_power_surface structure
	
        cx=problem_surfaces(surface_number)%face_list(output_face)%cell%i
        cy=problem_surfaces(surface_number)%face_list(output_face)%cell%j
        cz=problem_surfaces(surface_number)%face_list(output_face)%cell%k
        face=problem_surfaces(surface_number)%face_list(output_face)%point
 
! check the side of the surface on which we want output and change if required,
	
        if      (face.eq.face_xmin) then      
	  if (.NOT.frequency_domain_power_surface(output_surface)%output_on_outward_normal) then ! output on other side of face
	    cx=cx-1
	    face=face_xmax
	  end if
         else if (face.eq.face_xmax) then	 
	  if (.NOT.frequency_domain_power_surface(output_surface)%output_on_outward_normal) then ! output on other side of face
	    cx=cx+1
	    face=face_xmin
	  end if
        else if (face.eq.face_ymin) then	 
	  if (.NOT.frequency_domain_power_surface(output_surface)%output_on_outward_normal) then ! output on other side of face
	    cy=cy-1
	    face=face_ymax
	  end if
        else if (face.eq.face_ymax) then	 
	  if (.NOT.frequency_domain_power_surface(output_surface)%output_on_outward_normal) then ! output on other side of face
	    cy=cy+1
	    face=face_ymin
	  end if
        else if (face.eq.face_zmin) then	 
	  if (.NOT.frequency_domain_power_surface(output_surface)%output_on_outward_normal) then ! output on other side of face
	    cz=cz-1
	    face=face_zmax
	  end if
        else if (face.eq.face_zmax) then	 
	  if (.NOT.frequency_domain_power_surface(output_surface)%output_on_outward_normal) then ! output on other side of face
	    cz=cz+1
	    face=face_zmin
	  end if
        end if
  	     
! Set the output face in the local_surface_output array
! We must set the output point number on the min face 

        output_face1%cell%i=cx
        output_face1%cell%j=cy
        output_face1%cell%k=cz
        output_face1%point=face
	
        CALL get_min_face(output_face1,output_face2)
	
        cx=output_face2%cell%i
        cy=output_face2%cell%j
        cz=output_face2%cell%k
	face=output_face2%point

        if (rank.eq.cell_face_rank(cz,face)) then
! output point belongs to this processor
	  
          local_surface_output(cx  ,cy  ,cz  ,face)=1 
  
        end if ! output point belongs to this processor

! preserve the face which includes the normal information in the frequency_domain_power_surface%face_list	
        frequency_domain_power_surface(output_surface)%face_list(output_face)=output_face1
                 
      end do !next cell face in this surface	  
  
    end do ! next frequency_output surface

  end if ! n_frequency_domain_power_surface.GT.0 
  
  CALL write_line('FINISHED: set_frequency_domain_power_surfaces_in_mesh',0,output_to_screen_flag)
  
END SUBROUTINE set_frequency_domain_power_surfaces_in_mesh
!
! NAME
!     initialise_frequency_domain_power_surfaces
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/2/2013 CJS
!
!
SUBROUTINE initialise_frequency_domain_power_surfaces

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

  integer	:: output_surface
  integer	:: cx,cy,cz,face
  integer	:: number_of_faces
  integer	:: number_of_frequencies
  integer	:: output_face
  integer	:: n_frames

  integer 	:: face_number
  integer 	:: number_of_points
  integer 	:: point_number
  type(xyz)	:: point1,point2,point3,point4

! START
  
  if (n_frequency_domain_power_surfaces.gt.0) then

! OUTPUT SURFACES
    do output_surface=1,n_frequency_domain_power_surfaces
    
      number_of_faces=frequency_domain_power_surface(output_surface)%number_of_faces
      
      write(*,*)'Frequency output surface number',output_surface,' number of cell_faces=',number_of_faces
      
      number_of_frequencies=frequency_domain_power_surface(output_surface)%n_frequencies
      
      ALLOCATE( frequency_domain_power_surface(output_surface)%face_output_field_number_list(1:number_of_faces) )
      ALLOCATE( frequency_domain_power_surface(output_surface)%E1(1:number_of_faces,1:number_of_frequencies) )
      ALLOCATE( frequency_domain_power_surface(output_surface)%E2(1:number_of_faces,1:number_of_frequencies) )
      ALLOCATE( frequency_domain_power_surface(output_surface)%H1(1:number_of_faces,1:number_of_frequencies) )
      ALLOCATE( frequency_domain_power_surface(output_surface)%H2(1:number_of_faces,1:number_of_frequencies) )
      
      frequency_domain_power_surface(output_surface)%E1(1:number_of_faces,1:number_of_frequencies)=(0d0,0d0)
      frequency_domain_power_surface(output_surface)%E2(1:number_of_faces,1:number_of_frequencies)=(0d0,0d0)
      frequency_domain_power_surface(output_surface)%H1(1:number_of_faces,1:number_of_frequencies)=(0d0,0d0)
      frequency_domain_power_surface(output_surface)%H2(1:number_of_faces,1:number_of_frequencies)=(0d0,0d0)
      
      ALLOCATE( frequency_domain_power_surface(output_surface)%Power(1:number_of_frequencies) )
      frequency_domain_power_surface(output_surface)%Power(1:number_of_frequencies)=(0d0,0d0)
      
      do output_face=1,number_of_faces
   
        cx  =frequency_domain_power_surface(output_surface)%face_list(output_face)%cell%i
        cy  =frequency_domain_power_surface(output_surface)%face_list(output_face)%cell%j
        cz  =frequency_domain_power_surface(output_surface)%face_list(output_face)%cell%k
        face=frequency_domain_power_surface(output_surface)%face_list(output_face)%point 
	
	if (rank.eq.cell_face_rank(cz,face)) then
   	     
          if      (face.eq.face_xmin) then      
            frequency_domain_power_surface(output_surface)%face_output_field_number_list(output_face)=	&
	    			local_surface_output(cx  ,cy  ,cz  ,face_xmin) 
          else if (face.eq.face_xmax) then	 
            frequency_domain_power_surface(output_surface)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx+1,cy  ,cz  ,face_xmin)
          else if (face.eq.face_ymin) then	 
            frequency_domain_power_surface(output_surface)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx  ,cy  ,cz  ,face_ymin)
          else if (face.eq.face_ymax) then	 
            frequency_domain_power_surface(output_surface)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx  ,cy+1,cz  ,face_ymin)
          else if (face.eq.face_zmin) then	 
            frequency_domain_power_surface(output_surface)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx  ,cy  ,cz  ,face_zmin)
          else if (face.eq.face_zmax) then	 
            frequency_domain_power_surface(output_surface)%face_output_field_number_list(output_face)=	&
	  			local_surface_output(cx  ,cy  ,cz+1,face_zmin)
          end if
	
	end if ! cell belongs to this process
		  
      end do !next cell face in this surface	  
  
    end do ! next output surface

  end if ! n_frequency_domain_power_surfaces.GT.0
  

  RETURN

END SUBROUTINE initialise_frequency_domain_power_surfaces

!
! SUBROUTINE face_output_frequency_domain_power_surfaces
!
! NAME
!     face_output_frequency_domain_power_surfaces
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/2/2013 CJS
!
!
SUBROUTINE face_output_frequency_domain_power_surfaces

USE TLM_general
USE mesh
USE TLM_output
USE file_information
USE output_formats
USE cell_parameters
USE constants

IMPLICIT NONE

! local variables

  integer 	:: output_surface
  integer 	:: output_face
  integer 	:: output_field_number
  integer 	:: number_of_faces
  integer 	:: face
  integer 	:: cz
  integer	:: side
  
  real*8	:: E1,E2,H1,H2
  
  integer	:: frequency_loop
  real*8	:: frequency
  complex*16	:: ejwt

! START
  
  CALL write_line('CALLED: face_output_frequency_domain_power_surfaces',0,timestepping_output_to_screen_flag)

! loop over surfaces  
  do output_surface=1,n_frequency_domain_power_surfaces
      
    number_of_faces=frequency_domain_power_surface(output_surface)%number_of_faces
    
! loop over frequency 
    
    do frequency_loop=1,frequency_domain_power_surface(output_surface)%n_frequencies
    
      frequency=                   frequency_domain_power_surface(output_surface)%fmin+	&
                (frequency_loop-1)*frequency_domain_power_surface(output_surface)%fstep
	      
      ejwt=exp(-j*2d0*pi*frequency*time)       ! *dt  ! remove dt normalisation here
    
      do output_face=1,number_of_faces
         
        output_field_number=frequency_domain_power_surface(output_surface)%face_output_field_number_list(output_face)
        face=frequency_domain_power_surface(output_surface)%face_list(output_face)%point
        cz  =frequency_domain_power_surface(output_surface)%face_list(output_face)%cell%k
      
        if (rank.eq.cell_face_rank(cz,face)) then
! the output point is in this process so set the value   	     
          if	(face.eq.face_xmin) then
	    side=1
            E1=face_output_field(output_field_number,side,Ey)
            E2=face_output_field(output_field_number,side,Ez)
            H1=face_output_field(output_field_number,side,Hy)
            H2=face_output_field(output_field_number,side,Hz)
          else if (face.eq.face_xmax) then
	    side=2
            E1=face_output_field(output_field_number,side,Ey)
            E2=face_output_field(output_field_number,side,Ez)
            H1=face_output_field(output_field_number,side,Hy)
            H2=face_output_field(output_field_number,side,Hz)
          else if (face.eq.face_ymin) then
	    side=1	
            E1=face_output_field(output_field_number,side,Ez)
            E2=face_output_field(output_field_number,side,Ex)
            H1=face_output_field(output_field_number,side,Hz)
            H2=face_output_field(output_field_number,side,Hx)
          else if (face.eq.face_ymax) then
	    side=2	
            E1=face_output_field(output_field_number,side,Ez)
            E2=face_output_field(output_field_number,side,Ex)
            H1=face_output_field(output_field_number,side,Hz)
            H2=face_output_field(output_field_number,side,Hx)
          else if (face.eq.face_zmin) then
	    side=1	
            E1=face_output_field(output_field_number,side,Ex)
            E2=face_output_field(output_field_number,side,Ey)
            H1=face_output_field(output_field_number,side,Hx)
            H2=face_output_field(output_field_number,side,Hy)
          else if (face.eq.face_zmax) then
	    side=2	
            E1=face_output_field(output_field_number,side,Ex)
            E2=face_output_field(output_field_number,side,Ey)
            H1=face_output_field(output_field_number,side,Hx)
            H2=face_output_field(output_field_number,side,Hy)
          end if
 	
	  frequency_domain_power_surface(output_surface)%E1(output_face,frequency_loop)=	&
	      frequency_domain_power_surface(output_surface)%E1(output_face,frequency_loop)+E1*ejwt   
 	
	  frequency_domain_power_surface(output_surface)%E2(output_face,frequency_loop)=	&
	      frequency_domain_power_surface(output_surface)%E2(output_face,frequency_loop)+E2*ejwt   
 	
	  frequency_domain_power_surface(output_surface)%H1(output_face,frequency_loop)=	&
	      frequency_domain_power_surface(output_surface)%H1(output_face,frequency_loop)+H1*ejwt   
 	
	  frequency_domain_power_surface(output_surface)%H2(output_face,frequency_loop)=	&
	      frequency_domain_power_surface(output_surface)%H2(output_face,frequency_loop)+H2*ejwt   

        end if ! output face belongs to this process
          
      end do !next cell face in this surface	  
  
    end do ! next frequency
  
  end do ! next frequency_domain_power_surface

  
  CALL write_line('FINISHED: face_output_frequency_domain_power_surfaces',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE face_output_frequency_domain_power_surfaces

! Name write_frequency_domain_power_surfaces
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
!     started 8/2/2013 CJS adapted from Fieldsolve
!

SUBROUTINE write_frequency_domain_power_surfaces

USE file_information
USE output_formats
USE TLM_output
USE mesh
Use TLM_general
Use constants

IMPLICIT NONE

! variables passed to subroutine

! local_variables

  integer :: output_surface
  
  integer :: output_face
  
  integer :: cx,cy,cz,face
  
  integer :: sign
  
  integer	:: number_of_faces
  integer	:: number_of_frequencies,frequency_loop
  
  real*8	:: frequency
  
  complex*16 :: E1,E2,H1,H2,Pwr
 
  integer,allocatable	:: n_surfaces_rank(:)
  
  integer	:: talk_to  
  integer	:: n_complex
  complex*16,allocatable:: local_power(:)

! START

! write frequency_domain_power_surface_data

  if (n_frequency_domain_power_surfaces.ne.0) then

    write(*,*)'CALLED: write_frequency_domain_power_surfaces'

! open the output file
    if (rank.eq.0) then

      OPEN(unit=frequency_domain_power_surface_unit,file=trim(problem_name)//frequency_domain_power_surface_extn)
      write(frequency_domain_power_surface_unit,'(A)')'# Frequency Domain Data'
      write(frequency_domain_power_surface_unit,'(A,I10)')'# Number of data points:',n_frequency_domain_power_surfaces
      write(frequency_domain_power_surface_unit,'(A,I10)')'# Number of frequencies:',	&
            frequency_domain_power_surface(1)%n_frequencies

    end if
    
! loop over surfaces and calculate the power out of each surface  
    do output_surface=1,n_frequency_domain_power_surfaces
      
      number_of_faces=frequency_domain_power_surface(output_surface)%number_of_faces
    
! loop over frequency 
    
      do frequency_loop=1,frequency_domain_power_surface(output_surface)%n_frequencies
      
        Pwr=(0d0,0d0)
	
        do output_face=1,number_of_faces
         
          face=frequency_domain_power_surface(output_surface)%face_list(output_face)%point
          cz  =frequency_domain_power_surface(output_surface)%face_list(output_face)%cell%k
      
          if (rank.eq.cell_face_rank(cz,face)) then
! the output point is in this process so set the value   	     
            if	(face.eq.face_xmin) then
	      sign=1
            else if (face.eq.face_xmax) then
	      sign=-1
            else if (face.eq.face_ymin) then
	      sign=1
            else if (face.eq.face_ymax) then
	      sign=-1
            else if (face.eq.face_zmin) then
	      sign=1
            else if (face.eq.face_zmax) then
	      sign=-1
            end if
 	
	    E1=frequency_domain_power_surface(output_surface)%E1(output_face,frequency_loop)   	
	    E2=frequency_domain_power_surface(output_surface)%E2(output_face,frequency_loop)  
	    H1=frequency_domain_power_surface(output_surface)%H1(output_face,frequency_loop)   	
	    H2=frequency_domain_power_surface(output_surface)%H2(output_face,frequency_loop)  
	  
	    Pwr=Pwr+sign*(E1*conjg(H2)-E2*conjg(H1))*dl*dl/2d0

          end if ! output face belongs to this process
          
        end do !next cell face in this surface
      	  
        frequency_domain_power_surface(output_surface)%Power(frequency_loop)=Pwr
        
      end do ! next frequency
  
    end do ! next frequency_domain_power_surface

! Sum the power across all of the processes

    do output_surface=1,n_frequency_domain_power_surfaces

      number_of_frequencies=frequency_domain_power_surface(output_surface)%n_frequencies
      
      n_complex=number_of_frequencies
      ALLOCATE( local_Power(1:n_complex) )

#if defined(MPI)

      CALL MPI_REDUCE(frequency_domain_power_surface(output_surface)%Power, local_Power, &
                      n_complex, MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD,ierror)

#elif defined(SEQ)

      local_Power=frequency_domain_power_surface(output_surface)%Power

#endif

      if (rank.eq.0) then
! write frequency domain power to file    

        do frequency_loop=1,number_of_frequencies
	
          frequency=               frequency_domain_power_surface(output_surface)%fmin+	&
                (frequency_loop-1)*frequency_domain_power_surface(output_surface)%fstep
		
          write(frequency_domain_power_surface_unit,frequency_domain_output_format)				&
	        frequency,output_surface,dble(local_Power(frequency_loop)),dimag(local_Power(frequency_loop)),	&
		dabs(dble(local_Power(frequency_loop))),							&
		atan2( imag(local_Power(frequency_loop)),dble(local_Power(frequency_loop)) ),			&
		10d0*log10( dabs(dble(local_Power(frequency_loop))) )
		
	end do  

      end if ! rank.eq.0
      
      if (allocated( local_Power ) ) deallocate( local_Power )
      
    end do ! next output surface
  
    CLOSE(unit=frequency_domain_power_surface_unit)

    write(*,*)'FINISHED: write_frequency_domain_power_surfaces'
 
  end if !  n_frequency_domain_power_surfaces.ne.0
  
  RETURN
  
END SUBROUTINE write_frequency_domain_power_surfaces

