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
! SUBROUTINE set_huygens_surface_in_mesh
!
! NAME
!     set_huygens_surface_in_mesh
!
! DESCRIPTION
!     loop over all the required excitations and flag all the excitation cells/ faces 
!     in the arrays local_cell_excitation(i,j,k) or local_surface_excitation(i,j,k,face)
!     as required
!     
! COMMENTS
!  
!
! HISTORY
!
!     started 17/09/2012 CJS
!     Parallel 23/11/2012 CJS
!     Huygens surface 4/12/2012 CJS
!
!
SUBROUTINE set_huygens_surface_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_excitation
USE file_information

IMPLICIT NONE

! local variables

  integer	:: excitation_point
  integer	:: excitation_surface
  integer	:: surface_number
  integer	:: number_of_faces
  integer	:: excitation_face
  integer	:: huygens_face
  integer	:: n_faces
  integer	:: cx,cy,cz,face
  type(cell_point)	:: excitation_face1
  type(cell_point)	:: excitation_face2
    
  logical inside

! START
  
  CALL write_line('CALLED: set_huygens_surface_in_mesh',0,output_to_screen_flag)  
  
! SET HUYGENS SURFACE FACES  
  
  if (n_huygens_surfaces.ne.0) then
  
    if (.NOT.huygens_surface%outer_surface_flag) then ! not outer boundary

      surface_number=huygens_surface%surface_number

      n_faces=huygens_surface%n_surface_patches
      
      do huygens_face=1,n_faces
      
        cx=huygens_surface%cx(huygens_face)
	cy=huygens_surface%cy(huygens_face)
	cz=huygens_surface%cz(huygens_face)
	face=huygens_surface%face(huygens_face)
	
! Set the excitation face in the local_surface_excitation array
! We must set the excitation point number on the min face 

        excitation_face1%cell%i=cx
        excitation_face1%cell%j=cy
        excitation_face1%cell%k=cz
        excitation_face1%point=face
	
        CALL get_min_face(excitation_face1,excitation_face2)
	
        cx=excitation_face2%cell%i
        cy=excitation_face2%cell%j
        cz=excitation_face2%cell%k
	face=excitation_face2%point

! note, the face is always on the min surface of the cell
	
        if (rank.eq.cell_face_rank(cz,face)) then
! huygens surface face belongs to this processor	  
          local_surface_excitation(cx  ,cy  ,cz  ,face)=1 	
	end if
     
      end do ! next huygens surface face
  
    end if ! .NOT.huygens_surface%outer_surface_flag
  
  end if ! n_huygens_surfaces.ne.0
  
  CALL write_line('FINISHED: set_huygens_surface_in_mesh',0,output_to_screen_flag)

  RETURN

END SUBROUTINE set_huygens_surface_in_mesh
!
! NAME
!     initialise_huygens_surface
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
!     Huygens surface 4/12/2012 CJS
!
!
SUBROUTINE initialise_huygens_surface

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

  integer	:: cx,cy,cz,face
  integer	:: number_of_faces
  integer	:: excitation_face
  integer	:: huygens_face

  integer 	:: face_number
  integer 	:: number_of_points
  integer 	:: point_number
  type(xyz)	:: point1,point2,point3,point4
    
! START

  
  CALL write_line('CALLED: initialise_huygens_surface',0,output_to_screen_flag)

! HUYGENS SURFACES  
   
  if (n_huygens_surfaces.ne.0) then
  
    if (.NOT.huygens_surface%outer_surface_flag) then ! not outer boundary
  
      number_of_faces=huygens_surface%n_surface_patches
      
      if (rank.eq.0) then
        write(info_file_unit,*)'Huygens surface number',1,' Number of cell_faces=',number_of_faces
      end if
            
      do huygens_face=1,number_of_faces
      
        cx  =huygens_surface%cx(huygens_face)
        cy  =huygens_surface%cy(huygens_face)
        cz  =huygens_surface%cz(huygens_face)
        face=huygens_surface%face(huygens_face)
      
        if (rank.eq.cell_face_rank(cz,face)) then

          if      (face.eq.face_xmin) then      
            huygens_surface%face_excitation_field_number_list(huygens_face)=	&
	  			local_surface_excitation(cx  ,cy  ,cz  ,face_xmin) 
          else if (face.eq.face_xmax) then	 
            huygens_surface%face_excitation_field_number_list(huygens_face)=	&
	  			local_surface_excitation(cx+1,cy  ,cz  ,face_xmin)
          else if (face.eq.face_ymin) then	 
            huygens_surface%face_excitation_field_number_list(huygens_face)=	&
	  			local_surface_excitation(cx  ,cy  ,cz  ,face_ymin)
          else if (face.eq.face_ymax) then	 
            huygens_surface%face_excitation_field_number_list(huygens_face)=	&
	  			local_surface_excitation(cx  ,cy+1,cz  ,face_ymin)
          else if (face.eq.face_zmin) then	 
            huygens_surface%face_excitation_field_number_list(huygens_face)=	&
	  			local_surface_excitation(cx  ,cy  ,cz  ,face_zmin)
          else if (face.eq.face_zmax) then	 
            huygens_surface%face_excitation_field_number_list(huygens_face)=	&
	  			local_surface_excitation(cx  ,cy  ,cz+1,face_zmin)
          end if
	  
	  if (huygens_surface%face_excitation_field_number_list(huygens_face).eq.0) then
	    write(*,*)'Error in initialise_huygens_surface'
	    write(*,*)'No excitation surface found'
	    write(*,*)cx,cy,cz,face
	    STOP
	  end if
	  
	end if ! excitation point in this processor's mesh
		  
      end do !next cell face in the huygens surface	  
  
    end if ! .NOT.huygens_surface%outer_surface_flag
  
  end if ! n_huygens_surfaces.ne.0
  
  
  CALL write_line('FINISHED: initialise_huygens_surface',0,output_to_screen_flag)

  RETURN

END SUBROUTINE initialise_huygens_surface
!
! NAME
!     set_huygens_surface_data
!
! DESCRIPTION
!     
!     
! COMMENTS
!     Based on the Fieldsolve Huygens surface code.
!     I don't think we need to copy cell face coordinates and normals from the surface
!     structure... Even if we do, should it be the min face? Only use offset and normal direction in the update.
!
!     Thinking about applying Huygens surfaces separately on the two sides of a surface - not yet done
!     and still needs to be worked out.
!
!
! HISTORY
!
!     started 4/12/2012 CJS
!
!
SUBROUTINE set_huygens_surface_data

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

  integer	:: cx,cy,cz,face
  integer	:: number_of_faces
  
  real*8 r,theta,phi
  real*8 dxdt,dxdp
  real*8 dydt,dydp
  real*8 dzdt,dzdp
  real*8 vx,vy,vz
  real*8 kx,ky,kz
  real*8 norm
  real*8 Ei(3),Hi(3)
  
  integer	:: surface_number
  integer	:: cell_face
  integer	:: excitation_face
  integer	:: n_faces
  integer	:: face_count
  
  type(cell_point)	:: huygens_face
  type(xyz)		:: face_centre_coordinate
  type(xyz)		:: Kvector
  
  real*8 offset_min,offset
  real*8 hxmin,hxmax,hymin,hymax,hzmin,hzmax
  real*8 x,y,z
  
! START
  
  CALL write_line('CALLED: set_huygens_surface_data',0,output_to_screen_flag)
  
  if (n_huygens_surfaces.ne.0) then
  
    if (rank.eq.0) write(info_file_unit,*)'____________________________________________________'
    if (rank.eq.0) write(info_file_unit,*)''
    if (rank.eq.0) write(info_file_unit,*)'Huygens surface'
    if (rank.eq.0) write(info_file_unit,*)''

! calculate incident vector
    r=1d0
    theta=huygens_surface%Ktheta
    phi=huygens_surface%Kphi
    
    CALL rthetaphi_to_xyz_point(r,theta,phi,Kvector)
    
    kx=Kvector%x
    ky=Kvector%y
    kz=Kvector%z
    
    huygens_surface%Ki(1)=kx
    huygens_surface%Ki(2)=ky
    huygens_surface%Ki(3)=kz
    
    if (rank.eq.0) write(info_file_unit,8000)'K vector=',huygens_surface%Ki(1),huygens_surface%Ki(2),huygens_surface%Ki(3)
8000 format(A,3F12.6)

    dxdt= r*cos(theta)*cos(phi)
    dxdp=-r*sin(theta)*sin(phi)
    dydt= r*cos(theta)*sin(phi)
    dydp= r*sin(theta)*cos(phi)
    dzdt=-r*sin(theta)
    dzdp=0d0
    vx=huygens_surface%Etheta*dxdt+huygens_surface%Ephi*dxdp
    vy=huygens_surface%Etheta*dydt+huygens_surface%Ephi*dydp
    vz=huygens_surface%Etheta*dzdt+huygens_surface%Ephi*dzdp
    
    norm=sqrt(vx*vx+vy*vy+vz*vz)
    
    if (norm.eq.0d0) then
      write(*,*)'Error calculating Huygens surface data'
      write(*,8000)'E vector(theta)=',dxdt,dydt,dzdt
      write(*,8000)'E vector(phi)=',dxdp,dydp,dzdp
      write(*,8000)'E vector=',huygens_surface%Etheta,huygens_surface%Ephi
      STOP
    end if

    Ei(1)=vx/norm
    Ei(2)=vy/norm
    Ei(3)=vz/norm    
!note K cross E = H

    call vector_product(kx,ky,kz,Ei(1),Ei(2),Ei(3),Hi(1),Hi(2),Hi(3))
    Hi(:)=Hi(:)/Z0
    
    huygens_surface%Ei(1)=Ei(1)
    huygens_surface%Ei(2)=Ei(2)
    huygens_surface%Ei(3)=Ei(3)
    huygens_surface%Hi(1)=Hi(1)
    huygens_surface%Hi(2)=Hi(2)
    huygens_surface%Hi(3)=Hi(3)
    
    if (rank.eq.0) write(info_file_unit,8000)'E vector=',huygens_surface%Ei(1),huygens_surface%Ei(2),huygens_surface%Ei(3)
    if (rank.eq.0) write(info_file_unit,8000)'H vector=',huygens_surface%Hi(1),huygens_surface%Hi(2),huygens_surface%Hi(3)

! work out the number of points in the huygens surface , face_material holds the surface number
! so flag the surfaces whose surface number = huyges_surface_number

    if (huygens_surface%outer_surface_flag) then ! outer boundary Huygens surface
    
      n_faces=0
      
! xmin and xmax faces
      n_faces=n_faces+2*(nz2-nz1+1)*ny
  
! ymin and ymax faces
      n_faces=n_faces+2*(nz2-nz1+1)*nx
    
! zmin face
      if (rank.eq.0) then
        n_faces=n_faces+nx*ny
      end if
    
! zmax face
      if (rank.eq.np-1) then
        n_faces=n_faces+nx*ny
      end if
      
    else
    
      surface_number=huygens_surface%surface_number
      n_faces=problem_surfaces(surface_number)%number_of_faces
    
    end if
    
! allocate memeory for huygens surface data

    huygens_surface%n_surface_patches=n_faces
  
    if (rank.eq.0) write(info_file_unit,*)'Number of surface patches=',huygens_surface%n_surface_patches
  
    ALLOCATE( huygens_surface%offset(1:huygens_surface%n_surface_patches) )
    ALLOCATE( huygens_surface%cx(1:huygens_surface%n_surface_patches) )
    ALLOCATE( huygens_surface%cy(1:huygens_surface%n_surface_patches) )
    ALLOCATE( huygens_surface%cz(1:huygens_surface%n_surface_patches) )
    ALLOCATE( huygens_surface%face(1:huygens_surface%n_surface_patches) )
    ALLOCATE( huygens_surface%nx(1:huygens_surface%n_surface_patches) )
    ALLOCATE( huygens_surface%ny(1:huygens_surface%n_surface_patches) )
    ALLOCATE( huygens_surface%nz(1:huygens_surface%n_surface_patches) )        
    ALLOCATE( huygens_surface%face_excitation_field_number_list(1:huygens_surface%n_surface_patches) )

! Loop over huygens surface faces setting cell faces and normal directions
! Also work out the excitation function offset   
      
    if (.NOT.huygens_surface%outer_surface_flag) then ! not outer boundary

      surface_number=huygens_surface%surface_number

      n_faces=huygens_surface%n_surface_patches
      
      do excitation_face=1,n_faces

! get the correct cell and face for the excitation    
        cx=problem_surfaces(surface_number)%face_list(excitation_face)%cell%i
        cy=problem_surfaces(surface_number)%face_list(excitation_face)%cell%j
        cz=problem_surfaces(surface_number)%face_list(excitation_face)%cell%k
        face=problem_surfaces(surface_number)%face_list(excitation_face)%point
 
! check the side of the surface on which we want the excitation and reverse if required,
        if      (face.eq.face_xmin) then      
	  if (.NOT.huygens_surface%excitation_on_outward_normal) then ! excitation on other side of face
	    cx=cx-1
	    face=face_xmax
	  end if
         else if (face.eq.face_xmax) then	 
	  if (.NOT.huygens_surface%excitation_on_outward_normal) then ! excitation on other side of face
	    cx=cx+1
	    face=face_xmin
	  end if
        else if (face.eq.face_ymin) then	 
	  if (.NOT.huygens_surface%excitation_on_outward_normal) then ! excitation on other side of face
	    cy=cy-1
	    face=face_ymax
	  end if
        else if (face.eq.face_ymax) then	 
	  if (.NOT.huygens_surface%excitation_on_outward_normal) then ! excitation on other side of face
	    cy=cy+1
	    face=face_ymin
	  end if
        else if (face.eq.face_zmin) then	 
	  if (.NOT.huygens_surface%excitation_on_outward_normal) then ! excitation on other side of face
	    cz=cz-1
	    face=face_zmax
	  end if
        else if (face.eq.face_zmax) then	 
	  if (.NOT.huygens_surface%excitation_on_outward_normal) then ! excitation on other side of face
	    cz=cz+1
	    face=face_zmin
	  end if
        end if
	
        huygens_face%cell%i=cx
        huygens_face%cell%j=cy
        huygens_face%cell%k=cz
        huygens_face%point=face

! extent of the huygens surface mesh in x,y and z	
	huygens_surface%xmin=problem_surfaces(surface_number)%mesh_xmin
	huygens_surface%xmax=problem_surfaces(surface_number)%mesh_xmax
	huygens_surface%ymin=problem_surfaces(surface_number)%mesh_ymin
	huygens_surface%ymax=problem_surfaces(surface_number)%mesh_ymax
	huygens_surface%zmin=problem_surfaces(surface_number)%mesh_zmin
	huygens_surface%zmax=problem_surfaces(surface_number)%mesh_zmax
	
     	huygens_surface%nx(excitation_face)=0
     	huygens_surface%ny(excitation_face)=0
     	huygens_surface%nz(excitation_face)=0
	
! Work out the normal direction
	
	if (face.eq.face_xmin) then
	  huygens_surface%nx(excitation_face)=-1
	else if (face.eq.face_xmax) then
	  huygens_surface%nx(excitation_face)= 1
	else if (face.eq.face_ymin) then
	  huygens_surface%ny(excitation_face)=-1
	else if (face.eq.face_ymax) then
	  huygens_surface%ny(excitation_face)= 1
	else if (face.eq.face_zmin) then
	  huygens_surface%nz(excitation_face)=-1
	else if (face.eq.face_zmax) then
	  huygens_surface%nz(excitation_face)= 1
	else
	  write(*,*)'Error setting up Huygens surface'
	  write(*,*)'Face not set'
	  STOP
	end if
		
        huygens_surface%cx(excitation_face)=huygens_face%cell%i
	huygens_surface%cy(excitation_face)=huygens_face%cell%j
	huygens_surface%cz(excitation_face)=huygens_face%cell%k
	huygens_surface%face(excitation_face)=huygens_face%point
	
! work out the excitation function offset

	CALL get_cell_point_coordinate(huygens_face,face_centre_coordinate)
    	offset=kx*face_centre_coordinate%x+	&
	       ky*face_centre_coordinate%y+	&
	       kz*face_centre_coordinate%z
	
    	huygens_surface%offset(excitation_face)=offset

      end do ! next face belonging to the Huygens surface
      
    else  ! outer boundary

! extent of huygens surface mesh in x,y and z    
      huygens_surface%xmin=mesh_xmin
      huygens_surface%xmax=mesh_xmax
      huygens_surface%ymin=mesh_ymin
      huygens_surface%ymax=mesh_ymax
      huygens_surface%zmin=mesh_zmin
      huygens_surface%zmax=mesh_zmax
    
      face_count=0
! xmin and xmax faces
      do cz=nz1,nz2
        do cy=1,ny
    
          cx=1   
	         
	  face_count=face_count+1
          huygens_surface%cx(face_count)=cx
	  huygens_surface%cy(face_count)=cy
	  huygens_surface%cz(face_count)=cz
	  huygens_surface%face(face_count)=face_xmin	
     	  huygens_surface%nx(face_count)=-1
     	  huygens_surface%ny(face_count)=0
     	  huygens_surface%nz(face_count)=0
    
          cx=nx    
	         
	  face_count=face_count+1
          huygens_surface%cx(face_count)=cx
	  huygens_surface%cy(face_count)=cy
	  huygens_surface%cz(face_count)=cz
	  huygens_surface%face(face_count)=face_xmax	
     	  huygens_surface%nx(face_count)=1
     	  huygens_surface%ny(face_count)=0
     	  huygens_surface%nz(face_count)=0
		
        end do    ! next y cell
      end do      ! next z cell
  
! ymin and ymax faces
      do cz=nz1,nz2
        do cx=1,nx
    
          cy=1
	         
	  face_count=face_count+1
          huygens_surface%cx(face_count)=cx
	  huygens_surface%cy(face_count)=cy
	  huygens_surface%cz(face_count)=cz
	  huygens_surface%face(face_count)=face_ymin	
     	  huygens_surface%nx(face_count)=0
     	  huygens_surface%ny(face_count)=-1
     	  huygens_surface%nz(face_count)=0
      
          cy=ny
	         
	  face_count=face_count+1
          huygens_surface%cx(face_count)=cx
	  huygens_surface%cy(face_count)=cy
	  huygens_surface%cz(face_count)=cz
	  huygens_surface%face(face_count)=face_ymax	
     	  huygens_surface%nx(face_count)=0
     	  huygens_surface%ny(face_count)=1
     	  huygens_surface%nz(face_count)=0
     	  
        end do  ! next x cell
      end do      ! next z cell
    
! zmin face
      cz=1
  
      if (rank.eq.0) then
  
        do cy=1,ny
          do cx=1,nx
      	         
	    face_count=face_count+1
            huygens_surface%cx(face_count)=cx
	    huygens_surface%cy(face_count)=cy
	    huygens_surface%cz(face_count)=cz
	    huygens_surface%face(face_count)=face_zmin	
     	    huygens_surface%nx(face_count)=0
     	    huygens_surface%ny(face_count)=0
     	    huygens_surface%nz(face_count)=-1
 	  
          end do  ! next x cell
        end do    ! next y cell
    
      end if   ! rank.eq.0
    
! zmax face
      cz=nz
  
      if (rank.eq.np-1) then
    
        do cy=1,ny
          do cx=1,nx
      	         
	    face_count=face_count+1
            huygens_surface%cx(face_count)=cx
	    huygens_surface%cy(face_count)=cy
	    huygens_surface%cz(face_count)=cz
	    huygens_surface%face(face_count)=face_zmax	
     	    huygens_surface%nx(face_count)=0
     	    huygens_surface%ny(face_count)=0
     	    huygens_surface%nz(face_count)=1     
  	  
          end do  ! next x cell
        end do    ! next y cell
    
      end if ! rank.eq.np-1
      
! check
      if (face_count.ne.n_faces) then
        write(*,*)'Error in set_huygens_surface_data'
	write(*,*)'Face_count.ne.n_faces'
	write(*,*)'face_count=',face_count,' n_faces=',n_faces
	STOP
      end if
	
! work out the excitation function offset
      
      do face=1,n_faces
      
        huygens_face%cell%i=huygens_surface%cx(face)
        huygens_face%cell%j=huygens_surface%cy(face)
        huygens_face%cell%k=huygens_surface%cz(face)
	huygens_face%point =huygens_surface%face(face)
	
	CALL get_cell_point_coordinate(huygens_face,face_centre_coordinate)
    	offset=kx*face_centre_coordinate%x+	&
	       ky*face_centre_coordinate%y+	&
	       kz*face_centre_coordinate%z
	
    	huygens_surface%offset(face)=offset

      end do ! next face belonging to the Huygens surface

! in order to find offest_min we check the 12 corners of the mesh

    end if   ! outer boundary huygens surface

! the minimum offset is calculated from the 8 corners of the smallest cuboid containing the 
! huygens surface
    offset_min=1d30

! check corner 1    
    offset=kx*huygens_surface%xmin+	&
	   ky*huygens_surface%ymin+	&
	   kz*huygens_surface%zmin
    offset_min=min(offset,offset_min)

! check corner 2    
    offset=kx*huygens_surface%xmax+	&
	   ky*huygens_surface%ymin+	&
	   kz*huygens_surface%zmin
    offset_min=min(offset,offset_min)

! check corner 3    
    offset=kx*huygens_surface%xmin+	&
	   ky*huygens_surface%ymax+	&
	   kz*huygens_surface%zmin
    offset_min=min(offset,offset_min)

! check corner 4    
    offset=kx*huygens_surface%xmax+	&
	   ky*huygens_surface%ymax+	&
	   kz*huygens_surface%zmin
    offset_min=min(offset,offset_min)

! check corner 5    
    offset=kx*huygens_surface%xmin+	&
	   ky*huygens_surface%ymin+	&
	   kz*huygens_surface%zmax
    offset_min=min(offset,offset_min)

! check corner 6    
    offset=kx*huygens_surface%xmax+	&
	   ky*huygens_surface%ymin+	&
	   kz*huygens_surface%zmax
    offset_min=min(offset,offset_min)

! check corner 7   
    offset=kx*huygens_surface%xmin+	&
	   ky*huygens_surface%ymax+	&
	   kz*huygens_surface%zmax
    offset_min=min(offset,offset_min)

! check corner 8    
    offset=kx*huygens_surface%xmax+	&
	   ky*huygens_surface%ymax+	&
	   kz*huygens_surface%zmax
    offset_min=min(offset,offset_min)
      
    huygens_surface%offset_min=offset_min-dl/2d0
    
    if (rank.eq.0) write(info_file_unit,*)'Offset for source term=',huygens_surface%offset_min    
    if (rank.eq.0) write(info_file_unit,*)''

  end if ! n_huygens surfaces.ne.0
  
  CALL write_line('FINISHED: set_huygens_surface_data',0,output_to_screen_flag)

  RETURN

END SUBROUTINE set_huygens_surface_data
!
! SUBROUTINE Huygens_surface_excitation
!
! NAME
!     Huygens_surface_excitation
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started  4/12/2012 CJS
!
!
SUBROUTINE Huygens_surface_excitation

USE TLM_general
USE mesh
USE TLM_excitation
USE constants

IMPLICIT NONE

! local variables

  integer 	:: huygens_face
  integer 	:: face,side
  integer 	:: cz
  integer 	:: number_of_faces
  integer 	:: function_number
  integer 	:: excitation_array_point
  integer 	:: field_component
  real*8	:: offset,offset_min
  real*8 	:: value  
  
  real*8	:: Js(3),Ms(3)
  real*8	:: normx,normy,normz

! START
  
  CALL write_line('CALLED: Huygens_surface_excitation',0,timestepping_output_to_screen_flag)
  
  if (n_huygens_surfaces.ne.0) then
  
    if (.NOT.huygens_surface%outer_surface_flag) then ! not outer boundary
  
      number_of_faces=huygens_surface%n_surface_patches
            
      do huygens_face=1,number_of_faces
      
        cz  =huygens_surface%cz(huygens_face)
        face  =huygens_surface%face(huygens_face)
      
        if (rank.eq.cell_face_rank(cz,face)) then

! calculate the value of all 6 field components on this cell face

          function_number=huygens_surface%excitation_function_number
	  
	  offset    =huygens_surface%offset(huygens_face)
	  offset_min=huygens_surface%offset_min
	  
	  CALL get_interpolated_excitation_value(offset,offset_min,function_number,value,huygens_face)
	
          normx=huygens_surface%nx(huygens_face)
          normy=huygens_surface%ny(huygens_face)
          normz=huygens_surface%nz(huygens_face)
	  
          Js(1)= (normy*value*(huygens_surface%Hi(3))-normz*value*(huygens_surface%Hi(2)))
          Js(2)= (normz*value*(huygens_surface%Hi(1))-normx*value*(huygens_surface%Hi(3)))
          Js(3)= (normx*value*(huygens_surface%Hi(2))-normy*value*(huygens_surface%Hi(1)))
    
          Ms(1)=-(normy*value*(huygens_surface%Ei(3))-normz*value*(huygens_surface%Ei(2)))
          Ms(2)=-(normz*value*(huygens_surface%Ei(1))-normx*value*(huygens_surface%Ei(3)))
          Ms(3)=-(normx*value*(huygens_surface%Ei(2))-normy*value*(huygens_surface%Ei(1)))  
 	     
          if	(face.eq.face_xmin) then
	    side=1
          else if (face.eq.face_xmax) then
	    side=2
          else if (face.eq.face_ymin) then
	    side=1	
          else if (face.eq.face_ymax) then
	    side=2	
          else if (face.eq.face_zmin) then
	    side=1	
          else if (face.eq.face_zmax) then
	    side=2	
          end if

   
	  excitation_array_point=huygens_surface%face_excitation_field_number_list(huygens_face)

! We will need something like the following if we ever allow a Huygens surface to be applied
! to only one side of a boundary face...
!	face_excitation_field(excitation_array_point,side,Jx)=Js(1)
!	face_excitation_field(excitation_array_point,side,Jy)=Js(2)
!	face_excitation_field(excitation_array_point,side,Jz)=Js(3)
!	face_excitation_field(excitation_array_point,side,Mx)=Ms(1)
!	face_excitation_field(excitation_array_point,side,My)=Ms(2)
!	face_excitation_field(excitation_array_point,side,Mz)=Ms(3)
		
! Set the same excitation both sides of the face at the moment - need to think about this...
	face_excitation_field(excitation_array_point,1,Jx)=Js(1)
	face_excitation_field(excitation_array_point,2,Jx)=Js(1)
	face_excitation_field(excitation_array_point,1,Jy)=Js(2)
	face_excitation_field(excitation_array_point,2,Jy)=Js(2)
	face_excitation_field(excitation_array_point,1,Jz)=Js(3)
	face_excitation_field(excitation_array_point,2,Jz)=Js(3)
	face_excitation_field(excitation_array_point,1,Mx)=Ms(1)
	face_excitation_field(excitation_array_point,2,Mx)=Ms(1)
	face_excitation_field(excitation_array_point,1,My)=Ms(2)
	face_excitation_field(excitation_array_point,2,My)=Ms(2)
	face_excitation_field(excitation_array_point,1,Mz)=Ms(3)
	face_excitation_field(excitation_array_point,2,Mz)=Ms(3)
	  
	end if ! excitation point in this processor's mesh
		  
      end do !next cell face in this surface	  
  
    end if ! .NOT.huygens_surface%outer_surface_flag
  
  end if ! n_huygens_surfaces.ne.0
  
  CALL write_line('FINISHED: Huygens_surface_excitation',0,timestepping_output_to_screen_flag)
  
  RETURN

END SUBROUTINE Huygens_surface_excitation
