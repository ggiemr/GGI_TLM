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
! SUBROUTINE outer_boundary
!
! NAME
!     outer_boundary
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
SUBROUTINE outer_boundary

USE TLM_general
USE mesh
USE TLM_excitation
USE constants

IMPLICIT NONE

! local variables

! START


! local_variables

  integer cx,cy,cz
  
  real*8 Vx,Vy,Vz
  
! Huygens surface stuff

  integer	:: function_number
  integer	:: huygens_face
  
  real*8	:: offset,offset_min
  real*8 	:: value  
  
  real*8	:: Js(3),Ms(3)
  real*8	:: normx,normy,normz

! START
  
  CALL write_line('CALLED: outer_boundary',0,timestepping_output_to_screen_flag)
    
! xmin and xmax faces
  do cz=nz1,nz2
    do cy=1,ny
    
      cx=1          
      V(Vy_xmin,cx,cy,cz)=V(Vy_xmin,cx,cy,cz)*R_xmin
      V(Vz_xmin,cx,cy,cz)=V(Vz_xmin,cx,cy,cz)*R_xmin
    
      cx=nx    
      V(Vy_xmax,cx,cy,cz)=V(Vy_xmax,cx,cy,cz)*R_xmax
      V(Vz_xmax,cx,cy,cz)=V(Vz_xmax,cx,cy,cz)*R_xmax
		
    end do    ! next y cell
  end do      ! next z cell
  
! ymin and ymax faces
  do cz=nz1,nz2
    do cx=1,nx
    
      cy=1
      V(Vx_ymin,cx,cy,cz)=V(Vx_ymin,cx,cy,cz)*R_ymin
      V(Vz_ymin,cx,cy,cz)=V(Vz_ymin,cx,cy,cz)*R_ymin
    
      cy=ny
      V(Vx_ymax,cx,cy,cz)=V(Vx_ymax,cx,cy,cz)*R_ymax
      V(Vz_ymax,cx,cy,cz)=V(Vz_ymax,cx,cy,cz)*R_ymax
    	  
    end do  ! next x cell
  end do      ! next z cell
    
! zmin face
  cz=1
  
  if (rank.eq.0) then
  
    do cy=1,ny
      do cx=1,nx
      
        V(Vx_zmin,cx,cy,cz)=V(Vx_zmin,cx,cy,cz)*R_zmin
        V(Vy_zmin,cx,cy,cz)=V(Vy_zmin,cx,cy,cz)*R_zmin
  	  
      end do  ! next x cell
    end do    ! next y cell
    
  end if  
    
! zmax face
  cz=nz
  
  if (rank.eq.np-1) then
    
    do cy=1,ny
      do cx=1,nx
      
        V(Vx_zmax,cx,cy,cz)=V(Vx_zmax,cx,cy,cz)*R_zmax
        V(Vy_zmax,cx,cy,cz)=V(Vy_zmax,cx,cy,cz)*R_zmax
  	  
      end do  ! next x cell
    end do    ! next y cell
    
  end if  
  
! ADD HUYGENS SURFACE SOURCE TERMS IF REQUIRED
    
  if (n_huygens_surfaces.eq.0) RETURN
  
  if (huygens_surface%outer_surface_flag) then ! outer boundary huygens surface...

    function_number=huygens_surface%excitation_function_number
    offset_min=huygens_surface%offset_min

! note that these loops mimic those used in initialise_huygens_surface

   huygens_face=0
   
! xmin and xmax faces
    do cz=nz1,nz2
      do cy=1,ny
    
        cx=1         
	huygens_face=huygens_face+1 
	offset=huygens_surface%offset(huygens_face)
	
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
	
        V(Vy_xmin,cx,cy,cz)=V(Vy_xmin,cx,cy,cz)-Js(2)*dl*Z0/2d0-Ms(3)*dl/2d0
        V(Vz_xmin,cx,cy,cz)=V(Vz_xmin,cx,cy,cz)-Js(3)*dl*Z0/2d0+Ms(2)*dl/2d0
    
        cx=nx    
	huygens_face=huygens_face+1 
	offset=huygens_surface%offset(huygens_face)
	
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

        V(Vy_xmax,cx,cy,cz)=V(Vy_xmax,cx,cy,cz)-Js(2)*dl*Z0/2d0+Ms(3)*dl/2d0
        V(Vz_xmax,cx,cy,cz)=V(Vz_xmax,cx,cy,cz)-Js(3)*dl*Z0/2d0-Ms(2)*dl/2d0
		
      end do    ! next y cell
    end do      ! next z cell
  
! ymin and ymax faces
    do cz=nz1,nz2
      do cx=1,nx
    
        cy=1
	huygens_face=huygens_face+1 
	offset=huygens_surface%offset(huygens_face)
	
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

        V(Vx_ymin,cx,cy,cz)=V(Vx_ymin,cx,cy,cz)-Js(1)*dl*Z0/2d0+Ms(3)*dl/2d0
        V(Vz_ymin,cx,cy,cz)=V(Vz_ymin,cx,cy,cz)-Js(3)*dl*Z0/2d0-Ms(1)*dl/2d0
    
        cy=ny
	huygens_face=huygens_face+1 
	offset=huygens_surface%offset(huygens_face)
	
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

        V(Vx_ymax,cx,cy,cz)=V(Vx_ymax,cx,cy,cz)-Js(1)*dl*Z0/2d0-Ms(3)*dl/2d0
        V(Vz_ymax,cx,cy,cz)=V(Vz_ymax,cx,cy,cz)-Js(3)*dl*Z0/2d0+Ms(1)*dl/2d0
    	  
      end do  ! next x cell
    end do      ! next z cell
    
! zmin face
    cz=1
  
    if (rank.eq.0) then
  
      do cy=1,ny
        do cx=1,nx
      
	  huygens_face=huygens_face+1 
	  offset=huygens_surface%offset(huygens_face)
	
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

          V(Vx_zmin,cx,cy,cz)=V(Vx_zmin,cx,cy,cz)-Js(1)*dl*Z0/2d0-Ms(2)*dl/2d0
          V(Vy_zmin,cx,cy,cz)=V(Vy_zmin,cx,cy,cz)-Js(2)*dl*Z0/2d0+Ms(1)*dl/2d0
  	  
        end do  ! next x cell
      end do    ! next y cell
    
    end if  
    
! zmax face
    cz=nz
  
    if (rank.eq.np-1) then
    
      do cy=1,ny
        do cx=1,nx
      
	  huygens_face=huygens_face+1 
	  offset=huygens_surface%offset(huygens_face)
	
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

          V(Vx_zmax,cx,cy,cz)=V(Vx_zmax,cx,cy,cz)-Js(1)*dl*Z0/2d0+Ms(2)*dl/2d0
          V(Vy_zmax,cx,cy,cz)=V(Vy_zmax,cx,cy,cz)-Js(2)*dl*Z0/2d0-Ms(1)*dl/2d0
  	  
        end do  ! next x cell
      end do    ! next y cell
    
    end if  
    
! check that we have updated all the faces
    if (huygens_face.ne.huygens_surface%n_surface_patches) then
      write(*,*)'Error in initialise_huygens_surface'
      write(*,*)'Face_count.ne.n_faces'
      write(*,*)'face_count=',huygens_face,' n_faces=',huygens_surface%n_surface_patches
      STOP
    end if
  
  end if ! huygens_surface%outer_surface_flag

  CALL write_line('CALLED: outer_boundary',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE outer_boundary
