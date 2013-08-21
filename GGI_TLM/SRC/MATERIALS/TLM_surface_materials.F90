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
! SUBROUTINE set_surface_material_mesh
! SUBROUTINE calculate_surface_material_filter_coefficients
! SUBROUTINE allocate_surface_material_filter_data
! SUBROUTINE surface_material_update
!
! NAME
!     set_surface_material_mesh
!
! DESCRIPTION
!     
!     Count the number of faces of each material and construct the 
!     appropriate lists for the main solver
!     
! COMMENTS
!     Allocate a local mesh array and fill with material data. This ensures that
!     a face doesn't get given more than one material property.
!     In overlapping surfaces the last material assignment takes precidence
!
! HISTORY
!
!     started 14/08/2012 CJS
!
!     11/09/2012 - define local_surface_material on min faces taking
!                  proper account of non-symmetric materials and surface normals
!     frequency warping included 12/02/2013 CJS
!
!
SUBROUTINE set_surface_material_mesh

USE TLM_general
USE geometry
USE TLM_surface_materials
USE mesh
USE cell_parameters
USE file_information

IMPLICIT NONE

! local variables

  integer	:: material_number
  integer	:: surface_number
  integer	:: face,number_of_faces
  integer	:: total_number_of_material_faces
  integer	:: total_number_of_material_faces_all_procsesses
  integer	:: n_data
  integer	:: cx,cy,cz,cell_face
  
  integer	:: i
  integer	:: face_number
  integer	:: sign

! START
  
  CALL write_line('CALLED: set_surface_material_mesh',0,output_to_screen_flag)

! INITIALISE SURFACE MATERIALS  
    
  CALL write_line_integer('n_surface_materials=',n_surface_materials,0,output_to_screen_flag)

  if (rank.eq.0) then
  
    write(info_file_unit,*)'____________________________________________________'
    write(info_file_unit,*)''
    write(info_file_unit,*)'Number of surface materials=',n_surface_materials
    write(info_file_unit,*)''
    write(info_file_unit,*)'Material_number  Total_number_of_material_faces'
  
  end if

  do material_number=1,n_surface_materials
  
    total_number_of_material_faces=0
    
    CALL write_line_integer('Surface material number=',material_number,0,output_to_screen_flag)
    CALL write_line_integer('Number of geometric surfaces=',	&
                             surface_material_list(material_number)%n_surfaces,0,output_to_screen_flag)

! loop over the geometric surfaces with this material type and initialise the TLM faces
    do i=1,surface_material_list(material_number)%n_surfaces
    
      surface_number=surface_material_list(material_number)%surface_list(i)
      sign=surface_material_list(material_number)%surface_orientation_list(i)

      CALL write_line_integer('Geometric surface material number=',surface_number,0,output_to_screen_flag)
      
      problem_surfaces(surface_number)%surface_material_number=material_number
      
      number_of_faces=problem_surfaces(surface_number)%number_of_faces

      CALL write_line_integer('Number of faces=',number_of_faces,0,output_to_screen_flag)
      
      total_number_of_material_faces=total_number_of_material_faces+number_of_faces
      
      do face=1,number_of_faces
    
	cx=problem_surfaces(surface_number)%face_list(face)%cell%i
	cy=problem_surfaces(surface_number)%face_list(face)%cell%j
	cz=problem_surfaces(surface_number)%face_list(face)%cell%k
	cell_face=problem_surfaces(surface_number)%face_list(face)%point
	
!	write(*,*)cx,cy,cz,cell_face

! give the cell face the material number and the opposite face -material number
! this allows us to keep track of the orientation of the surface and also means that
! there is no ambiguity if a face and its opposite are given different material numbers - the
! surface property is set by the last specification
! Note we only need the min surfaces to be set, the opposite side is negative.
   	     
	if ( (cell_face.eq.face_xmin).AND.(rank.eq.cell_face_rank(cz,cell_face)) ) then      
          local_surface_material(cx  ,cy  ,cz  ,face_xmin)= sign*material_number	 
	    
	else if ( (cell_face.eq.face_xmax).AND.(rank.eq.cell_face_rank(cz,cell_face)) ) then      
          local_surface_material(cx+1,cy  ,cz  ,face_xmin)=-sign*material_number	
	     
	else if ( (cell_face.eq.face_ymin).AND.(rank.eq.cell_face_rank(cz,cell_face)) ) then      
          local_surface_material(cx  ,cy  ,cz  ,face_ymin)= sign*material_number	
	     
	else if ( (cell_face.eq.face_ymax).AND.(rank.eq.cell_face_rank(cz,cell_face)) ) then      
          local_surface_material(cx  ,cy+1,cz  ,face_ymin)=-sign*material_number	
	     
	else if ( (cell_face.eq.face_zmin).AND.(rank.eq.cell_face_rank(cz,cell_face)) ) then      
          local_surface_material(cx  ,cy  ,cz  ,face_zmin)= sign*material_number	 
	    
	else if ( (cell_face.eq.face_zmax).AND.(rank.eq.cell_face_rank(cz+1,face_zmin)) ) then      
          local_surface_material(cx  ,cy  ,cz+1,face_zmin)=-sign*material_number	   
	  	   
        end if
	 	
      end do ! next face
    
    end do ! next surface number

#if defined(MPI)

    n_data=1
    call MPI_REDUCE(total_number_of_material_faces,total_number_of_material_faces_all_procsesses , &
                    n_data,MPI_INTEGER, MPI_SUM, 0,MPI_COMM_WORLD,ierror)
		    
#elif defined(SEQ)

    total_number_of_material_faces_all_procsesses=total_number_of_material_faces
    
#endif			

    if (rank.eq.0) then
      CALL write_line_integer('Total number of material faces=',	&
                              total_number_of_material_faces_all_procsesses,0,output_to_screen_flag)
      write(info_file_unit,'(2I14)')material_number,total_number_of_material_faces_all_procsesses
    end if

  end do ! next surface material number
  
  CALL write_line('FINISHED: set_surface_material_mesh',0,output_to_screen_flag)

  RETURN

END SUBROUTINE set_surface_material_mesh
!
! Name calculate_surface_material_filter_coefficients
!     
!
! Description:
!      Once the timestep is found we can determine the surface material filter coefficients
!      
!
! Comments:
!      
!      
!
! History
!
!     started 04/09/12 CJS
!

SUBROUTINE calculate_surface_material_filter_coefficients()

USE TLM_general
USE TLM_surface_materials
USE filter_types
USE filter_operators
USE filter_functions
USE constants

IMPLICIT NONE

! variables passed to subroutine

! local_variables

  integer material_number
  
  type(Zfilter) :: Zfilter_temp

! function_types

! START

  write(*,*)'CALLED: calculate_surface_material_filter_coefficients'

  do material_number=1,n_surface_materials
   
! Z transform of Z parameter filter	
    
    if (surface_material_list(material_number)%type.EQ.surface_material_type_DISPERSIVE) then

!      Zfilter_temp=s_to_z(surface_material_list(material_number)%Z11_S,dt)     
      Zfilter_temp=s_to_z_warp(surface_material_list(material_number)%Z11_S,dt,bicubic_warp_flag,frequency_scale)     
      call Z_fast_slow_docomposition( Zfilter_temp, 		&
                     surface_material_list(material_number)%Z11_f, 	&
                     surface_material_list(material_number)%Z11_Z    )
    
!      Zfilter_temp=s_to_z(surface_material_list(material_number)%Z12_S,dt)     
      Zfilter_temp=s_to_z_warp(surface_material_list(material_number)%Z12_S,dt,bicubic_warp_flag,frequency_scale)     
      call Z_fast_slow_docomposition( Zfilter_temp, 		&
                     surface_material_list(material_number)%Z12_f, 	&
                     surface_material_list(material_number)%Z12_Z    )
    
!      Zfilter_temp=s_to_z(surface_material_list(material_number)%Z21_S,dt)     
      Zfilter_temp=s_to_z_warp(surface_material_list(material_number)%Z21_S,dt,bicubic_warp_flag,frequency_scale)     
      call Z_fast_slow_docomposition( Zfilter_temp, 		&
                     surface_material_list(material_number)%Z21_f, 	&
                     surface_material_list(material_number)%Z21_Z    )
     
!      Zfilter_temp=s_to_z(surface_material_list(material_number)%Z22_S,dt)     
      Zfilter_temp=s_to_z_warp(surface_material_list(material_number)%Z22_S,dt,bicubic_warp_flag,frequency_scale)     
      call Z_fast_slow_docomposition( Zfilter_temp, 		&
                     surface_material_list(material_number)%Z22_f, 	&
                     surface_material_list(material_number)%Z22_Z    )
		     
    end if ! material type is dispersive
    
  end do ! next material to set

  write(*,*)'FINISHED: calculate_surface_material_filter_coefficients'
  return
  
END SUBROUTINE calculate_surface_material_filter_coefficients
!
! NAME
!     allocate_surface_material_filter_data
!
! DESCRIPTION
!     allocate the memory required for surface material filter data
!     
! COMMENTS
!
!
! HISTORY
!
!     started 4/09/2012 CJS
!
!
SUBROUTINE allocate_surface_material_filter_data

USE TLM_general
USE TLM_surface_materials
USE filter_types
USE filter_operators
USE filter_functions
USE mesh

IMPLICIT NONE

! local variables
  
  integer	:: material_number
  integer	:: face
  integer	:: surface_filter_number
  integer       :: z11aorder,z11border
  integer       :: z12aorder,z12border
  integer       :: z21aorder,z21border
  integer       :: z22aorder,z22border

! START
  
  CALL write_line('CALLED: allocate_surface_material_filter_data',0,output_to_screen_flag)
  
  n_surface_material_faces=0
  surface_material_storage=0

! loop over special faces
  DO face=1,n_special_faces
  
    if (face_update_code_to_material_data(face,1).ne.0) then
! this is a material face add to the memory allocation list if required

      material_number=abs(face_update_code_to_material_data(face,1))
      
      if ( (surface_material_list(material_number)%type.eq.surface_material_type_DISPERSIVE) ) then
	  	      
        n_surface_material_faces=n_surface_material_faces+1

        face_update_code_to_material_data(face,2)=surface_material_storage+1
        surface_material_storage=surface_material_storage+2   ! note 2 filters per face	    

      end if ! surface material at this cell face

    end if ! is this a material face? 

  end do ! next special face

! allocate surface material filter data storage for each filter required
  if (surface_material_storage.ne.0) then
    allocate (surface_material_Z11_filter_data(1:surface_material_storage)) 
    allocate (surface_material_Z12_filter_data(1:surface_material_storage)) 
    allocate (surface_material_Z21_filter_data(1:surface_material_storage)) 
    allocate (surface_material_Z22_filter_data(1:surface_material_storage)) 
  end if
 
! loop over filters allocating memory for the particular material order
 

! loop over special faces
  DO face=1,n_special_faces
  
    if (face_update_code_to_material_data(face,1).ne.0) then
! this is a material face add to the memory allocation list if required

      material_number=abs(face_update_code_to_material_data(face,1))
      
      if ( (surface_material_list(material_number)%type.eq.surface_material_type_DISPERSIVE) ) then
	  	      
	z11aorder=surface_material_list(material_number)%Z11_Z%a%order
	z11border=surface_material_list(material_number)%Z11_Z%b%order

	z12aorder=surface_material_list(material_number)%Z12_Z%a%order
	z12border=surface_material_list(material_number)%Z12_Z%b%order

	z21aorder=surface_material_list(material_number)%Z21_Z%a%order
	z21border=surface_material_list(material_number)%Z21_Z%b%order

	z22aorder=surface_material_list(material_number)%Z22_Z%a%order
	z22border=surface_material_list(material_number)%Z22_Z%b%order

! first polarisation 
        surface_filter_number=face_update_code_to_material_data(face,2)
	  
        surface_material_Z11_filter_data(surface_filter_number)=allocate_Zfilter_response(z11aorder,z11border)	
        surface_material_Z12_filter_data(surface_filter_number)=allocate_Zfilter_response(z12aorder,z12border)	
        surface_material_Z21_filter_data(surface_filter_number)=allocate_Zfilter_response(z21aorder,z21border)	
        surface_material_Z22_filter_data(surface_filter_number)=allocate_Zfilter_response(z22aorder,z22border)

! second polarisation 
        surface_filter_number=face_update_code_to_material_data(face,2)+1
	  
        surface_material_Z11_filter_data(surface_filter_number)=allocate_Zfilter_response(z11aorder,z11border)	
        surface_material_Z12_filter_data(surface_filter_number)=allocate_Zfilter_response(z12aorder,z12border)	
        surface_material_Z21_filter_data(surface_filter_number)=allocate_Zfilter_response(z21aorder,z21border)	
        surface_material_Z22_filter_data(surface_filter_number)=allocate_Zfilter_response(z22aorder,z22border)

      end if ! surface material at this cell face

    end if ! is this a material face? 

  end do ! next special face
 
 
  CALL write_line('FINISHED: allocate_surface_material_filter_data',0,output_to_screen_flag)


  RETURN

END SUBROUTINE allocate_surface_material_filter_data
!
! Name surface_material_update
!     
!
! Description
!    Surface material solution: Calculate the voltages on either side of
!    a surface material characterised by a frequency dependent impedance matrix.
!
!    Called with incident voltage pulses from either side, wave impedance either side, material
!    type number and surface filter number
!
!    Update of the surface filter is done here, timeshifting for the next timestep is also done   
!
! Comments:
!
! History
!
!     started 5/09/12 CJS
!

SUBROUTINE surface_material_update(Vli,Zl,Vri,Zr,Vl,Vr,material_number,reverse_material,surface_filter_number)

USE TLM_general
USE TLM_excitation
USE TLM_output
USE TLM_surface_materials
USE filter_types
USE filter_operators
USE filter_functions
USE constants


IMPLICIT NONE

! variables passed to subroutine

  real*8 	:: Vli,Zl,Vri,Zr,Vl,Vr
  integer	:: material_number
  logical	:: reverse_material
  integer	:: surface_filter_number

! local_variables

  real*8 z11f,z12f,z21f,z22f
  real*8 z11_i1s,z12_i2s,z21_i1s,z22_i2s
  
  real*8 U11,U12,U21,U22,i1,i2
  real*8 rhs1,rhs2
  real*8 det
  
  real*8 V1i,V2i,V1,V2,Z1,Z2
 
! function_types

! START
 
   Z11f=surface_material_list(material_number)%Z11_f
   Z12f=surface_material_list(material_number)%Z12_f
   Z21f=surface_material_list(material_number)%Z21_f
   Z22f=surface_material_list(material_number)%Z22_f

! Z11 slow response
   call evaluate_Zfilter(surface_material_list(material_number)%Z11_Z,		&
	                 surface_material_Z11_filter_data(surface_filter_number),	&
		         0d0)
   Z11_i1s=surface_material_Z11_filter_data(surface_filter_number)%f

! Z12 slow response
   call evaluate_Zfilter(surface_material_list(material_number)%Z12_Z,		&
	                 surface_material_Z12_filter_data(surface_filter_number),	&
		         0d0)
   Z12_i2s=surface_material_Z12_filter_data(surface_filter_number)%f

! Z21 slow response
   call evaluate_Zfilter(surface_material_list(material_number)%Z21_Z,		&
	                 surface_material_Z21_filter_data(surface_filter_number),	&
		         0d0)
   Z21_i1s=surface_material_Z21_filter_data(surface_filter_number)%f

! Z22 slow response
   call evaluate_Zfilter(surface_material_list(material_number)%Z22_Z,		&
	                 surface_material_Z22_filter_data(surface_filter_number),	&
		         0d0)
   Z22_i2s=surface_material_Z22_filter_data(surface_filter_number)%f

   if (.NOT.reverse_material) then
     
     Z1=Zl
     Z2=Zr
       
   else
     
     Z1=Zr
     Z2=Zl
     
   end if
 
   U11=Z1+Z11f
   U12=   Z12f
   U21=   Z21f
   U22=Z2+Z22f
   
! solve for port currents   
   
   det=U11*U22-U21*U12
   
   if (det.ne.0d0) then  ! go ahead and solve
   
     if (.NOT.reverse_material) then
     
       V1i=Vli
       V2i=Vri
       
     else
     
       V1i=Vri
       V2i=Vli
     
     end if
       
! calculate RHS   
     rhs1=2d0*v1i-Z11_i1s-Z12_i2s
     rhs2=2d0*v2i-Z21_i1s-Z22_i2s
     
! calculate port currents     
     i1=( U22*rhs1-U12*rhs2)/det
     i2=(-U21*rhs1+U11*rhs2)/det
   
! calculate port voltages    
     V1=2d0*v1i-i1*Z1
     V2=2d0*v2i-i2*Z2

! update filter data for next timestep     
     call evaluate_Zfilter(surface_material_list(material_number)%Z11_Z,		&
	                   surface_material_Z11_filter_data(surface_filter_number),	&
	  	           i1)
     call evaluate_Zfilter(surface_material_list(material_number)%Z12_Z,		&
	                   surface_material_Z12_filter_data(surface_filter_number),	&
		           i2)
     call evaluate_Zfilter(surface_material_list(material_number)%Z21_Z,		&
	                   surface_material_Z21_filter_data(surface_filter_number),	&
		           i1)
     call evaluate_Zfilter(surface_material_list(material_number)%Z22_Z,		&
	                   surface_material_Z22_filter_data(surface_filter_number),	&
		           i2)
   
     if (.NOT.reverse_material) then
     
       Vl=V1
       Vr=V2
       
     else
     
       Vl=V2
       Vr=V1
     
     end if
     
   else ! det.eq.0d0
   
     write(*,*)'ERROR: surface_material_update'
     write(*,*)'determinant=0'
     write(*,*)'Z11f=',Z11f
     write(*,*)'Z12f=',Z12f
     write(*,*)'Z21f=',Z21f
     write(*,*)'Z22f=',Z22f
     write(*,*)'Zl=',Zl
     write(*,*)'Zr=',Zr
     stop
     
   end if
   
! timeshift filter ready for the next timestep
   call timeshift_Zfilter(surface_material_Z11_filter_data(surface_filter_number))
   call timeshift_Zfilter(surface_material_Z12_filter_data(surface_filter_number))
   call timeshift_Zfilter(surface_material_Z21_filter_data(surface_filter_number))
   call timeshift_Zfilter(surface_material_Z22_filter_data(surface_filter_number))

!  write(*,*)'FINISHED: surface_material_update'
  return
  
END SUBROUTINE surface_material_update
