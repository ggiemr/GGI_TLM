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
! SUBROUTINE set_periodic_boundary_far_field_surfaces_in_mesh
! SUBROUTINE initialise_periodic_boundary_far_field_surfaces
! SUBROUTINE face_output_periodic_boundary_far_field
!
! NAME
!     set_periodic_boundary_far_field_surfaces_in_mesh
!
! DESCRIPTION
!
!     
! COMMENTS
!
! HISTORY
!
!     started 1/5/2014 CJS based on the HIRF-SE code Fieldsolve 
!
!
SUBROUTINE set_periodic_boundary_far_field_surfaces_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE file_information

IMPLICIT NONE

! local variables

  integer	:: surface_number
  integer	:: surface
  integer	:: number_of_faces
  integer	:: output_face,face_count
  integer	:: cx,cy,cz,face
  type(cell_point)	:: output_face1
  type(cell_point)	:: output_face2

! START
  
  CALL write_line('CALLED: set_periodic_boundary_far_field_surfaces_in_mesh',0,output_to_screen_flag)
  
  if (n_PB_far_field_surfaces.ne.0) then
  
    do surface=1,n_PB_far_field_surfaces
  
      surface_number=PB_far_field_surface(surface)%surface_number
    
      number_of_faces=problem_surfaces(surface_number)%number_of_faces/4
      PB_far_field_surface(surface)%number_of_faces=number_of_faces
      
! allocate face list for the output surface      
      ALLOCATE( PB_far_field_surface(surface)%face_list(1:number_of_faces) )

! note that the number of faces is filtered down by a factor of 4 here...      
      face_count=0
      do output_face=1,problem_surfaces(surface_number)%number_of_faces
      
! copy the faces from the geometry structure to the local
! frequency_output_surface structure
	
        cx=problem_surfaces(surface_number)%face_list(output_face)%cell%i
        cy=problem_surfaces(surface_number)%face_list(output_face)%cell%j
        cz=problem_surfaces(surface_number)%face_list(output_face)%cell%k
        face=problem_surfaces(surface_number)%face_list(output_face)%point

! check the side of the surface on which we want output and change if required,

        if      (face.eq.face_xmin) then      
          if (.NOT.PB_far_field_surface(surface)%output_on_outward_normal) then ! output on other side of face
            cx=cx-1
            face=face_xmax
          end if
         else if (face.eq.face_xmax) then        
          if (.NOT.PB_far_field_surface(surface)%output_on_outward_normal) then ! output on other side of face
            cx=cx+1
            face=face_xmin
          end if
        else if (face.eq.face_ymin) then         
          if (.NOT.PB_far_field_surface(surface)%output_on_outward_normal) then ! output on other side of face
            cy=cy-1
            face=face_ymax
          end if
        else if (face.eq.face_ymax) then         
          if (.NOT.PB_far_field_surface(surface)%output_on_outward_normal) then ! output on other side of face
            cy=cy+1
            face=face_ymin
          end if
        else if (face.eq.face_zmin) then         
          if (.NOT.PB_far_field_surface(surface)%output_on_outward_normal) then ! output on other side of face
            cz=cz-1
            face=face_zmax
          end if
        else if (face.eq.face_zmax) then         
          if (.NOT.PB_far_field_surface(surface)%output_on_outward_normal) then ! output on other side of face
            cz=cz+1
            face=face_zmin
          end if
        end if
	
	if ((cx.gt.nx/2).AND.(cy.gt.ny/2)) then
! we are in the appropriate sub-cell for output so add this cell to the output list  
  	     
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
        
	    local_surface_output(cx,cy,cz,face)=1 

          end if ! output point belongs to this processor

! preserve the face which includes the normal information in the PB_far_field_surface%face_list	
          face_count=face_count+1
          PB_far_field_surface(surface)%face_list(face_count)=output_face1
	       
        end if !we are in the appropriate sub-cell for output so add this face to the output list        
	       
      end do !next cell face in this surface	
	       
    end do ! periodic boundary far field surface	

  end if ! n_PB_far_field_surfaces.ne.0
  
  CALL write_line('FINISHED: set_periodic_boundary_far_field_surfaces_in_mesh',0,output_to_screen_flag)
  
END SUBROUTINE set_periodic_boundary_far_field_surfaces_in_mesh
!
! NAME
!     initialise_periodic_boundary_far_field_surfaces
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 1/5/2014 CJS based on the HIRF-SE code Fieldsolve 
!
!
SUBROUTINE initialise_periodic_boundary_far_field_surfaces

USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE TLM_excitation
USE constants
USE file_information

IMPLICIT NONE

! local variables

  integer	:: surface 
  integer	:: cx,cy,cz,face
  integer	:: number_of_faces
  integer	:: output_face
  
  type(xyz)	:: xyz_point
  type(xyz)	:: k
  
  real*8 	:: r,theta,phi,kx,ky,kz,phi_in
  real*8 	:: dxdt,dxdp,dydt,dydp,dzdt,dzdp,vx,vy,vz,norm
  
  real*8	:: local_dmin,local_dmax
  real*8	:: dmin,dmax,dmrcs
  integer	:: n_data
  
  real*8	:: far_field_time_period
  real*8	:: far_field_distance

! function variables
 
  real*8	:: xyz_dot

! START
  
  CALL write_line('CALLED: initialise_periodic_boundary_far_field_surfaces',0,timestepping_output_to_screen_flag)
  
  if (n_PB_far_field_surfaces.ne.0) then
  
    if (rank.eq.0) write(info_file_unit,*)'____________________________________________________'
    if (rank.eq.0) write(info_file_unit,*)''
    if (rank.eq.0) write(info_file_unit,*)'Periodic boundary far field surfaces'
    if (rank.eq.0) write(info_file_unit,*)''
  
    do surface=1,n_PB_far_field_surfaces
    
      if (rank.eq.0) write(info_file_unit,*)''
      if (rank.eq.0) write(info_file_unit,*)'Far field surface number',surface
      if (rank.eq.0) write(info_file_unit,*)''
    
! calculate incident field vector
      r    =1d0
      theta=huygens_surface%Ktheta
      phi  =huygens_surface%Kphi

! save initial phi so that polarisation ambiguity may be resolved for theta=0  
      phi_in=phi

      CALL rthetaphi_to_xyz_point(r,theta,phi,k)
    
      if (PB_far_field_surface(surface)%r_t_option.eq.'r') then
! define K vector for the reflected wave by reversing kz
    
        k%z=-k%z
    
      end if

      CALL xyz_point_to_rthetaphi(k,r,theta,phi)

! resolve polarisation ambiguity if theta=0  
      if (theta.eq.0d0) then
        phi=phi_in
      end if
         
      PB_far_field_surface(surface)%theta=theta
      PB_far_field_surface(surface)%phi=phi
        
      PB_far_field_surface(surface)%Kmrcs(1)=k%x
      PB_far_field_surface(surface)%Kmrcs(2)=k%y
      PB_far_field_surface(surface)%Kmrcs(3)=k%z

! calculate vectors perpendicular to K in the theta and phi directions
      if ((theta.eq.0d0).OR.(theta.eq.180d0)) then
        PB_far_field_surface(surface)%Vtheta(1)=1d0
        PB_far_field_surface(surface)%Vtheta(2)=0d0
        PB_far_field_surface(surface)%Vtheta(3)=0d0
        PB_far_field_surface(surface)%Vphi(1)=0d0
        PB_far_field_surface(surface)%Vphi(2)=1d0
        PB_far_field_surface(surface)%Vphi(3)=0d0
      else
        dxdt= r*cos(theta)*cos(phi)
        dxdp=-r*sin(theta)*sin(phi)
        dydt= r*cos(theta)*sin(phi)
        dydp= r*sin(theta)*cos(phi)
        dzdt=-r*sin(theta)
        dzdp=0d0
    
        vx=dxdt
        vy=dydt
        vz=dzdt
        norm=sqrt(vx*vx+vy*vy+vz*vz)
        if (norm.eq.0d0) then
          write(*,*)'Error calculating rcs Vtheta'
          write(*,8000)'vector(theta)=',dxdt,dydt,dzdt
          write(*,8000)'vector(phi)=',dxdp,dydp,dzdp
8000      format(A,3F12.6)
          STOP
        end if
        PB_far_field_surface(surface)%Vtheta(1)=vx/norm
        PB_far_field_surface(surface)%Vtheta(2)=vy/norm
        PB_far_field_surface(surface)%Vtheta(3)=vz/norm

        vx=dxdp
        vy=dydp
        vz=dzdp
        norm=sqrt(vx*vx+vy*vy+vz*vz)
        if (norm.eq.0d0) then
          write(*,*)'Error calculating rcs Vphi'
          write(*,8000)'vector(theta)=',dxdt,dydt,dzdt
          write(*,8000)'vector(phi)=',dxdp,dydp,dzdp
          STOP
        end if
        PB_far_field_surface(surface)%Vphi(1)=vx/norm
        PB_far_field_surface(surface)%Vphi(2)=vy/norm
        PB_far_field_surface(surface)%Vphi(3)=vz/norm

      end if
    
      number_of_faces=PB_far_field_surface(surface)%number_of_faces
     
      if (rank.eq.0) write(info_file_unit,*)'Number of PCB far field faces=',number_of_faces
    
      ALLOCATE( PB_far_field_surface(surface)%face_output_field_number_list(1:number_of_faces) )
    
       PB_far_field_surface(surface)%face_output_field_number_list(1:number_of_faces)=0
        
      local_dmin=1d30
      local_dmax=-1d30

      do output_face=1,number_of_faces
   
        cx  =PB_far_field_surface(surface)%face_list(output_face)%cell%i
        cy  =PB_far_field_surface(surface)%face_list(output_face)%cell%j
        cz  =PB_far_field_surface(surface)%face_list(output_face)%cell%k
        face=PB_far_field_surface(surface)%face_list(output_face)%point 

        if (rank.eq.cell_face_rank(cz,face)) then
           
    	  if	(face.eq.face_xmin) then      
    	    PB_far_field_surface(surface)%face_output_field_number_list(output_face)= &
        	  	      local_surface_output(cx  ,cy  ,cz  ,face_xmin) 
    	  else if (face.eq.face_xmax) then       
    	    PB_far_field_surface(surface)%face_output_field_number_list(output_face)= &
        	  	      local_surface_output(cx+1,cy  ,cz  ,face_xmin)
    	  else if (face.eq.face_ymin) then       
    	    PB_far_field_surface(surface)%face_output_field_number_list(output_face)= &
          		      local_surface_output(cx  ,cy  ,cz  ,face_ymin)
    	  else if (face.eq.face_ymax) then       
    	    PB_far_field_surface(surface)%face_output_field_number_list(output_face)= &
          		      local_surface_output(cx  ,cy+1,cz  ,face_ymin)
    	  else if (face.eq.face_zmin) then       
    	    PB_far_field_surface(surface)%face_output_field_number_list(output_face)= &
        		      local_surface_output(cx  ,cy  ,cz  ,face_zmin)
    	  else if (face.eq.face_zmax) then       
    	    PB_far_field_surface(surface)%face_output_field_number_list(output_face)= &
          		      local_surface_output(cx  ,cy  ,cz+1,face_zmin)
    	  end if

! work out the distance of this point along the K vector in the far field direction	    
	  CALL get_cell_point_coordinate(PB_far_field_surface(surface)%face_list(output_face),xyz_point)
	  dmrcs=xyz_dot(K,xyz_point)
	  local_dmin=min(dmrcs,local_dmin)
	  local_dmax=max(dmrcs,local_dmax)
	
        end if ! cell belongs to this process
        	
      end do !next cell face in this surface	

! work out dmin and dmax across all processes
  
#if defined(MPI)

      n_data=1
      call MPI_ALLREDUCE(local_dmin,dmin,n_data,MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD,ierror)
      call MPI_ALLREDUCE(local_dmax,dmax,n_data,MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD,ierror)

#elif defined(SEQ)

      dmin=local_dmin
      dmax=local_dmax
    
#endif

! work out the time period which we need to save for the far field in the monostatic direction

      if (rank.eq.0) write(info_file_unit,*)'(R.k)min=',dmin
      if (rank.eq.0) write(info_file_unit,*)'(R.k)max=',dmax

      far_field_time_period=(dmax-dmin)/c0+simulation_time
      far_field_distance=far_field_time_period*c0
      PB_far_field_surface(surface)%n_far_field_points=n_rcs_points_per_cell*far_field_distance/dl+2
      PB_far_field_surface(surface)%far_field_point_offset=n_rcs_points_per_cell*dmin/dl-1
    
      if (rank.eq.0) write(info_file_unit,*)'Far field time period=',far_field_time_period
      if (rank.eq.0) write(info_file_unit,*)'Far field distance	=',far_field_distance
      if (rank.eq.0) write(info_file_unit,*)'n_far_field_points=',PB_far_field_surface(surface)%n_far_field_points
      if (rank.eq.0) write(info_file_unit,*)'far_field_point_offset=',PB_far_field_surface(surface)%far_field_point_offset
      if (rank.eq.0) write(info_file_unit,*)''
    
      ALLOCATE( PB_far_field_surface(surface)%Etheta(1:PB_far_field_surface(surface)%n_far_field_points) )
      ALLOCATE( PB_far_field_surface(surface)%Htheta(1:PB_far_field_surface(surface)%n_far_field_points) )
      ALLOCATE( PB_far_field_surface(surface)%Ephi(1:PB_far_field_surface(surface)%n_far_field_points) )
      ALLOCATE( PB_far_field_surface(surface)%Hphi(1:PB_far_field_surface(surface)%n_far_field_points) )
    
      PB_far_field_surface(surface)%Etheta(1:PB_far_field_surface(surface)%n_far_field_points)=0d0
      PB_far_field_surface(surface)%Htheta(1:PB_far_field_surface(surface)%n_far_field_points)=0d0
      PB_far_field_surface(surface)%Ephi(1:PB_far_field_surface(surface)%n_far_field_points)  =0d0
      PB_far_field_surface(surface)%Hphi(1:PB_far_field_surface(surface)%n_far_field_points)  =0d0

      if (rank.eq.0) then
! rank 0 process only: open file for periodic boundary far field output
  
        OPEN(unit=PB_far_field_output_unit,file=trim(problem_name)//PB_far_field_output_extn)
  
        OPEN(unit=PB_far_field_time_output_unit,file=trim(problem_name)//PB_far_field_time_output_extn)

      end if ! rank 0 process
   
     end do ! next periodic boundary far field surface
   
   end if ! n_PB_far_field_surfaces.ne.0
      
   PB_far_field_output_warning=.FALSE.
  
  CALL write_line('FINISHED: initialise_periodic_boundary_far_field_surfaces',0,timestepping_output_to_screen_flag)
 
  RETURN

END SUBROUTINE initialise_periodic_boundary_far_field_surfaces

!
! SUBROUTINE face_output_periodic_boundary_far_field
!
! NAME
!     face_output_periodic_boundary_far_field
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 1/5/2014 CJS based on the HIRF-SE code Fieldsolve 
!
!
SUBROUTINE face_output_periodic_boundary_far_field

USE TLM_general
USE mesh
USE TLM_output
USE file_information
USE output_formats
USE cell_parameters
USE constants

IMPLICIT NONE

! local variables

  integer	:: surface
  integer 	:: output_surface
  integer 	:: output_face
  integer 	:: output_field_number
  
  integer 	:: number_of_faces
  integer 	:: face
  integer 	:: cz
  type(xyz)	:: xyz_point
  type(xyz)	:: K
  integer	:: side
  integer	:: normx,normy,normz
  real*8	:: field(6)
  
  real*8	:: Jsx,Jsy,Jsz,Msx,Msy,Msz
  real*8	:: theta,phi
  real*8	:: ct,cp,st,sp
  real*8	:: Ntheta,Nphi
  real*8	:: Ltheta,Lphi
  real*8	:: Etheta,Ephi
  real*8	:: Htheta,Hphi
    
  real*8 	:: dA
  
  integer	:: point
  real*8	:: dmrcs
  
  integer	:: i

! function variables

  real*8	:: xyz_dot

! START
  
  CALL write_line('CALLED: face_output_periodic_boundary_far_field',0,timestepping_output_to_screen_flag)

  if (n_PB_far_field_surfaces.eq.0) RETURN
  
  do surface=1,n_PB_far_field_surfaces
  
    dA=dl*dl
   
    theta=PB_far_field_surface(surface)%theta
    phi=PB_far_field_surface(surface)%phi
    
    ct=cos(theta)
    st=sin(theta)
    cp=cos(phi)
    sp=sin(phi)  
    
    K%x=PB_far_field_surface(surface)%Kmrcs(1)
    K%y=PB_far_field_surface(surface)%Kmrcs(2)
    K%z=PB_far_field_surface(surface)%Kmrcs(3)
    
    number_of_faces=PB_far_field_surface(surface)%number_of_faces
    
    do output_face=1,number_of_faces
    
      output_field_number=PB_far_field_surface(surface)%face_output_field_number_list(output_face)
               
      face=PB_far_field_surface(surface)%face_list(output_face)%point
      cz  =PB_far_field_surface(surface)%face_list(output_face)%cell%k
      
      normx=0
      normy=0
      normz=0
      
      if (rank.eq.cell_face_rank(cz,face)) then
! inside this processor's mesh
 	     
        if	(face.eq.face_xmin) then
	  side=1
	  normx=1
        else if (face.eq.face_xmax) then
	  side=2
	  normx=-1
        else if (face.eq.face_ymin) then
	  side=1	
	  normy=1
        else if (face.eq.face_ymax) then
	  side=2	
	  normy=-1
        else if (face.eq.face_zmin) then
	  side=1	
	  normz=1
        else if (face.eq.face_zmax) then
	  side=2	
	  normz=-1
        end if
	
        field(1:6)=face_output_field(output_field_number,side,1:6)  
            
        Jsx=   (normy*field(Hz)-normz*field(Hy))
        Jsy=   (normz*field(Hx)-normx*field(Hz))
        Jsz=   (normx*field(Hy)-normy*field(Hx))
      
        Msx= - (normy*field(Ez)-normz*field(Ey))
        Msy= - (normz*field(Ex)-normx*field(Ez))
        Msz= - (normx*field(Ey)-normy*field(Ex))
      
! work out the distance of this point along the K vector in the far field direction	
 
	CALL get_cell_point_coordinate(PB_far_field_surface(surface)%face_list(output_face),xyz_point)
	dmrcs=xyz_dot(K,xyz_point)

! add the offset distance due to the time	    
        dmrcs=-dmrcs+time*c0
	    
! calculate the point to update in the far field array
        point=INT(dble(n_rcs_points_per_cell)*dmrcs/dl)	   
	       
! add the offset distance due to the geometry of the surface	    
        point=point-PB_far_field_surface(surface)%far_field_point_offset
	    
        if ( (point.ge.1).AND.(point.le.PB_far_field_surface(surface)%n_far_field_points) ) then
! calculate far field and add contributions to this point
		
          Ntheta=( Jsx*ct*cp+Jsy*ct*sp-Jsz*st)*dA/(4d0*pi)
          Nphi  =(-Jsx*sp   +Jsy*cp	  )*dA/(4d0*pi)
	 
          Ltheta=( Msx*ct*cp+Msy*ct*sp-Msz*st)*dA/(4d0*pi)
          Lphi  =(-Msx*sp   +Msy*cp	  )*dA/(4d0*pi)
	
          Etheta=-(Lphi+Z0*Ntheta)
	  Ephi  = (Ltheta-Z0*Nphi)
	
	  Htheta= (Nphi-Ltheta/Z0)
	  Hphi  =-(Ntheta+Lphi/Z0)
	    
          PB_far_field_surface(surface)%Etheta(point)=PB_far_field_surface(surface)%Etheta(point)+Etheta
          PB_far_field_surface(surface)%Ephi(point)  =PB_far_field_surface(surface)%Ephi(point)  +Ephi
          PB_far_field_surface(surface)%Htheta(point)=PB_far_field_surface(surface)%Htheta(point)+Htheta
          PB_far_field_surface(surface)%Hphi(point)  =PB_far_field_surface(surface)%Hphi(point)  +Hphi
	
        else
          
	  if (.NOT.PB_far_field_output_warning) then
            PB_far_field_output_warning=.TRUE.
	    number_of_warnings=number_of_warnings+1
	    write(warning_file_unit,*)' '
	    write(warning_file_unit,*)'Periodic boundary far field output warning: Trying to put far field point outside array'
            write(warning_file_unit,*)'point=',point,' n_points=',PB_far_field_surface(surface)%n_far_field_points
	    write(warning_file_unit,*)' '
	  end if
	  
        end if
	 
      end if ! in this processor's mesh
      
    end do

  end do ! next periodic boundary far field surface    
  
  CALL write_line('FINISHED: face_output_periodic_boundary_far_field',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE face_output_periodic_boundary_far_field
!
! Name write_periodic_boundary_far_field_surfaces
!     
!
! Description
!     calculate far E field and output to file
!
! Comments:
!      
!
! History
!
!     started 1/5/2014 CJS based on the HIRF-SE code Fieldsolve 
!

SUBROUTINE write_periodic_boundary_far_field_surfaces

USE file_information
USE TLM_output
USE TLM_excitation
USE mesh
Use TLM_general
Use constants

IMPLICIT NONE

! variables passed to subroutine

! local_variables


  real*8, allocatable	:: local_Etheta(:)
  real*8, allocatable	:: local_Ephi(:)
  real*8, allocatable 	:: local_Htheta(:)
  real*8, allocatable	:: local_Hphi(:)
  
  real*8 		:: REtheta,REphi,RHtheta,RHphi
  complex*16 		:: Cpower
  real*8     		:: Rpower,Ipower

  integer		:: n_data
  integer		:: dmin,dmax,d_local
  
  integer		:: timestep_loop
  real*8		:: t_local,dt_local
  real*8		:: dA
  
  integer		:: n_surface_patches
  real*8		:: Area

  integer		:: n_frequencies,frequency_loop
  real*8		:: frequency
  
  integer		:: function_number
  
  real*8		:: normx=0.0
  real*8		:: normy=0.0
  real*8		:: normz=1.0
  
  real*8		:: theta,phi
  real*8		:: ct,cp,st,sp
	  
  real*8	 ::  Exi
  real*8	 ::  Eyi
  real*8	 ::  Ezi
   
  real*8	 ::  Hxi
  real*8	 ::  Hyi
  real*8	 ::  Hzi
     
  real*8	 ::  Jsxi
  real*8	 ::  Jsyi
  real*8	 ::  Jszi
  
  real*8	 ::   Msxi
  real*8	 ::   Msyi
  real*8	 ::   Mszi
  
  real*8	 ::   Nthetai
  real*8	 ::   Nphii  

  real*8	 ::   Lthetai
  real*8	 ::   Lphii  

  complex*16	 ::  Ethetai
  complex*16	 ::  Ephii  

  complex*16	 ::  Hthetai
  complex*16	 ::  Hphii  
  
  complex*16		:: Fourier_const
  complex*16 		:: Etheta,Ephi,Htheta,Hphi,Etot,Htot
  
  integer		:: i

! function_types
  
! START

  write(*,*)'CALLED: write_periodic_boundary_far_field_surfaces'


  if (n_PB_far_field_surfaces.eq.0) RETURN
  
  do i=1,n_PB_far_field_surfaces
      
    n_data=PB_far_field_surface(i)%n_far_field_points
    
    allocate( local_Etheta(1:n_data) )
    allocate( local_Ephi(1:n_data) )
    allocate( local_Htheta(1:n_data) )
    allocate( local_Hphi(1:n_data) )
  
#if defined(MPI)

    call MPI_REDUCE(PB_far_field_surface(i)%Etheta, local_Etheta, &
                    n_data,MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(PB_far_field_surface(i)%Ephi, local_Ephi, &
                    n_data,MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(PB_far_field_surface(i)%Htheta, local_Htheta, &
                    n_data,MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(PB_far_field_surface(i)%Hphi, local_Hphi, &
                    n_data,MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD,ierror)
		    
#elif defined(SEQ)

    local_Etheta=PB_far_field_surface(i)%Etheta
    local_Ephi  =PB_far_field_surface(i)%Ephi
    local_Htheta=PB_far_field_surface(i)%Htheta
    local_Hphi  =PB_far_field_surface(i)%Hphi

#endif

    if (rank.eq.0) then

! write down range response to file    
      dmin=1
      dmax=PB_far_field_surface(i)%n_far_field_points
      do d_local=dmin,dmax 
    
        t_local=(d_local-2)*dl/(n_rcs_points_per_cell*c0)-dt/2d0    

        REtheta=local_Etheta(d_local)
        REphi=local_Ephi(d_local)
        RHtheta=local_Htheta(d_local)
        RHphi=local_Hphi(d_local)
      
        if (abs(REtheta).lt.1e-30) REtheta=0d0
        if (abs(REphi).lt.1e-30)   REphi=0d0
        if (abs(RHtheta).lt.1e-30) RHtheta=0d0
        if (abs(RHphi).lt.1e-30)   RHphi=0d0
      
        write(pb_far_field_time_output_unit,8050)t_local,i,REtheta,REphi,RHtheta,RHphi
     
8050    format(E14.6,I5,4E14.6)
      end do ! next timestep
   
      dA=dl*dl
      dt_local=dl/(n_rcs_points_per_cell*c0)    
      n_surface_patches=PB_far_field_surface(i)%number_of_faces
      Area=dA*n_surface_patches
   
      theta=huygens_surface%Ktheta
      phi  =huygens_surface%Kphi
    
      ct=cos(theta)
      st=sin(theta)
      cp=cos(phi)
      sp=sin(phi)  

      n_frequencies=int( (PB_far_field_surface(i)%fmax-PB_far_field_surface(i)%fmin)/PB_far_field_surface(i)%fstep )+1

      do frequency_loop=1,n_frequencies

        frequency=PB_far_field_surface(i)%fmin+(frequency_loop-1)*PB_far_field_surface(i)%fstep
           
! calculate the incident field at the frequency of interest for RCS calculation

        function_number=huygens_surface%excitation_function_number
	
        Ethetai=(0d0,0d0)
	Ephii  =(0d0,0d0)
	
	Hthetai=(0d0,0d0)
	Hphii  =(0d0,0d0)
 
        do timestep_loop=1,n_timesteps

          t_local=(timestep_loop-1)*dt

! evaluate the function at the current time
          Fourier_const=dt*exp(-j*2d0*pi*frequency*t_local)
	  
	  normx=0.0
	  normy=0.0
	  normz=1.0
	  
          Exi=huygens_surface%Ei(1)*excitation_functions(function_number)%value(timestep_loop)
          Eyi=huygens_surface%Ei(2)*excitation_functions(function_number)%value(timestep_loop)
          Ezi=huygens_surface%Ei(3)*excitation_functions(function_number)%value(timestep_loop)
	  
          Hxi=huygens_surface%Hi(1)*excitation_functions(function_number)%value(timestep_loop)
          Hyi=huygens_surface%Hi(2)*excitation_functions(function_number)%value(timestep_loop)
          Hzi=huygens_surface%Hi(3)*excitation_functions(function_number)%value(timestep_loop)
            
          Jsxi=   (normy*Hzi-normz*Hyi)
          Jsyi=   (normz*Hxi-normx*Hzi)
          Jszi=   (normx*Hyi-normy*Hxi)
      
          Msxi= - (normy*Ezi-normz*Eyi)
          Msyi= - (normz*Exi-normx*Ezi)
          Mszi= - (normx*Eyi-normy*Exi)
      		
          Nthetai=( Jsxi*ct*cp+Jsyi*ct*sp-Jszi*st)*Area/(4d0*pi)
          Nphii  =(-Jsxi*sp   +Jsyi*cp	  )*Area/(4d0*pi)
	 
          Lthetai=( Msxi*ct*cp+Msyi*ct*sp-Mszi*st)*Area/(4d0*pi)
          Lphii  =(-Msxi*sp   +Msyi*cp	  )*Area/(4d0*pi)
	
          Ethetai=Ethetai-(Lphii+Z0*Nthetai)*Fourier_const
	  Ephii  =Ephii  +(Lthetai-Z0*Nphii)*Fourier_const
	
	  Hthetai=Hthetai+(Nphii-Lthetai/Z0)*Fourier_const
	  Hphii  =Hphii  -(Nthetai+Lphii/Z0)*Fourier_const
	  
        end do ! next timestep
     
! Fourier transform the far field

        Etheta=(0d0,0d0)
        Ephi=(0d0,0d0)
        Htheta=(0d0,0d0)
        Hphi=(0d0,0d0)
	
        dmin=1
        dmax=PB_far_field_surface(i)%n_far_field_points 
    
        do d_local=dmin,dmax 
	
          t_local=d_local*dl/(n_rcs_points_per_cell*c0)
! evaluate the function at the current time

          Fourier_const=dt_local*exp(-j*2d0*pi*frequency*t_local)
	  
	  Etheta=Etheta+local_Etheta(d_local)*Fourier_const
	  Ephi  =Ephi  +local_Ephi(d_local)  *Fourier_const
	  Htheta=Htheta+local_Htheta(d_local)*Fourier_const
	  Hphi  =Hphi  +local_Hphi(d_local)  *Fourier_const
	  
        end do ! next timestep
                  
! Poynting vector over transmission surface.      
          Cpower=Etheta*conjg(Hphi)-Ephi*conjg(Htheta)
          Rpower=abs(Cpower)
	  
! Poynting vector for incident wave over xy surface  
!          Ipower=mod_Einc**2*abs(huygens_surface%Ei(1)*huygens_surface%Hi(2)-huygens_surface%Ei(2)*huygens_surface%Hi(1))
          Ipower=abs(Ethetai*conjg(Hphii)-Ephii*conjg(Hthetai))
	                  
	  if ((Rpower.gt.0d0).AND.(Ipower.ne.0d0)) then
	    Rpower=Rpower/Ipower
	    Ipower=-10.0*log10(Rpower)
	  else
	    Rpower=0d0
	    Ipower=0d0
	  end if
         
          write(pb_far_field_output_unit,8060)i,frequency,abs(Etheta),abs(Ephi),sqrt(Real(Ipower)),	&
                                        sqrt(Rpower),Rpower,Ipower
8060      format(i5,7E12.4)
	
      end do ! next frequency
    
    end if   !(rank.eq.0)
    
    deallocate( local_Etheta )
    deallocate( local_Ephi )
    deallocate( local_Htheta )
    deallocate( local_Hphi )

  end do ! next periodic boundary far field surface    
  
  write(*,*)'FINISHED: write_periodic_boundary_far_field_surfaces'
  return
  
END SUBROUTINE write_periodic_boundary_far_field_surfaces

