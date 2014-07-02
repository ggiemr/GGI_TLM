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
! SUBROUTINE set_far_field_surfaces_in_mesh
! SUBROUTINE initialise_far_field_surfaces
! SUBROUTINE face_output_far_field
!
! NAME
!     set_far_field_surfaces_in_mesh
!
! DESCRIPTION
!
!     
! COMMENTS
!
! HISTORY
!
!     started 5/12/2012 CJS
!
!
SUBROUTINE set_far_field_surfaces_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE file_information

IMPLICIT NONE

! local variables

  integer	:: ff_surface
  integer	:: surface_number
  integer	:: number_of_faces
  integer	:: output_face
  integer	:: cx,cy,cz,face
  type(cell_point)	:: output_face1
  type(cell_point)	:: output_face2

! START
  
  CALL write_line('CALLED: set_far_field_surfaces_in_mesh',0,output_to_screen_flag)
  
  do ff_surface=1,n_far_field_surfaces
  
    surface_number=far_field_surface(ff_surface)%surface_number
    number_of_faces=problem_surfaces(surface_number)%number_of_faces
    far_field_surface(ff_surface)%number_of_faces=number_of_faces
      
! allocate face list for the output surface      
    ALLOCATE( far_field_surface(ff_surface)%face_list(1:number_of_faces) )
      
    do output_face=1,number_of_faces
      
! copy the faces from the geometry structure to the local
! frequency_output_surface structure
	
      cx=problem_surfaces(surface_number)%face_list(output_face)%cell%i
      cy=problem_surfaces(surface_number)%face_list(output_face)%cell%j
      cz=problem_surfaces(surface_number)%face_list(output_face)%cell%k
      face=problem_surfaces(surface_number)%face_list(output_face)%point

! check the side of the surface on which we want output and change if required,

      if      (face.eq.face_xmin) then      
        if (.NOT.far_field_surface(ff_surface)%output_on_outward_normal) then ! output on other side of face
          cx=cx-1
          face=face_xmax
        end if
       else if (face.eq.face_xmax) then        
        if (.NOT.far_field_surface(ff_surface)%output_on_outward_normal) then ! output on other side of face
          cx=cx+1
          face=face_xmin
        end if
      else if (face.eq.face_ymin) then         
        if (.NOT.far_field_surface(ff_surface)%output_on_outward_normal) then ! output on other side of face
          cy=cy-1
          face=face_ymax
        end if
      else if (face.eq.face_ymax) then         
        if (.NOT.far_field_surface(ff_surface)%output_on_outward_normal) then ! output on other side of face
          cy=cy+1
          face=face_ymin
        end if
      else if (face.eq.face_zmin) then         
        if (.NOT.far_field_surface(ff_surface)%output_on_outward_normal) then ! output on other side of face
          cz=cz-1
          face=face_zmax
        end if
      else if (face.eq.face_zmax) then         
        if (.NOT.far_field_surface(ff_surface)%output_on_outward_normal) then ! output on other side of face
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
        
	local_surface_output(cx,cy,cz,face)=1 

      end if ! output point belongs to this processor

! preserve the face which includes the normal information in the far_field_surface(ff_surface)%face_list	
      far_field_surface(ff_surface)%face_list(output_face)=output_face1
	       
    end do !next cell face in this surface	

  end do ! next far_field_surface
  
  CALL write_line('FINISHED: set_far_field_surfaces_in_mesh',0,output_to_screen_flag)
  
END SUBROUTINE set_far_field_surfaces_in_mesh
!
! NAME
!     initialise_far_field_surfaces
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 5/12/2012 CJS
!
!
SUBROUTINE initialise_far_field_surfaces

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

  integer	:: ff_surface
  integer	:: cx,cy,cz,face
  integer	:: number_of_faces
  integer	:: output_face

! START
  
  do ff_surface=1,n_far_field_surfaces
    
    number_of_faces=far_field_surface(ff_surface)%number_of_faces
    
    ALLOCATE( far_field_surface(ff_surface)%face_output_field_number_list(1:number_of_faces) )
    ALLOCATE( far_field_surface(ff_surface)%J(1:number_of_faces,1:3) )
    ALLOCATE( far_field_surface(ff_surface)%M(1:number_of_faces,1:3) )
    
    far_field_surface(ff_surface)%face_output_field_number_list(1:number_of_faces)=0
    far_field_surface(ff_surface)%J(1:number_of_faces,1:3)=(0d0,0d0)
    far_field_surface(ff_surface)%M(1:number_of_faces,1:3)=(0d0,0d0)
    
    do output_face=1,number_of_faces
   
      cx  =far_field_surface(ff_surface)%face_list(output_face)%cell%i
      cy  =far_field_surface(ff_surface)%face_list(output_face)%cell%j
      cz  =far_field_surface(ff_surface)%face_list(output_face)%cell%k
      face=far_field_surface(ff_surface)%face_list(output_face)%point 

      if (rank.eq.cell_face_rank(cz,face)) then
           
    	if	(face.eq.face_xmin) then      
    	  far_field_surface(ff_surface)%face_output_field_number_list(output_face)= &
        		      local_surface_output(cx  ,cy  ,cz  ,face_xmin) 
    	else if (face.eq.face_xmax) then       
    	  far_field_surface(ff_surface)%face_output_field_number_list(output_face)= &
        		      local_surface_output(cx+1,cy  ,cz  ,face_xmin)
    	else if (face.eq.face_ymin) then       
    	  far_field_surface(ff_surface)%face_output_field_number_list(output_face)= &
        		      local_surface_output(cx  ,cy  ,cz  ,face_ymin)
    	else if (face.eq.face_ymax) then       
    	  far_field_surface(ff_surface)%face_output_field_number_list(output_face)= &
        		      local_surface_output(cx  ,cy+1,cz  ,face_ymin)
    	else if (face.eq.face_zmin) then       
    	  far_field_surface(ff_surface)%face_output_field_number_list(output_face)= &
        		      local_surface_output(cx  ,cy  ,cz  ,face_zmin)
    	else if (face.eq.face_zmax) then       
    	  far_field_surface(ff_surface)%face_output_field_number_list(output_face)= &
        		      local_surface_output(cx  ,cy  ,cz+1,face_zmin)
    	end if

      end if ! cell belongs to this process
        	
    end do !next cell face in this surface	

  end do ! next far_field_surface

  RETURN

END SUBROUTINE initialise_far_field_surfaces

!
! SUBROUTINE face_output_far_field
!
! NAME
!     face_output_far_field
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
!     include dt in Fourier Integral to make the output consistent with other Fourier Transformed outputs 26/9/2013
!
!
SUBROUTINE face_output_far_field

USE TLM_general
USE mesh
USE TLM_output
USE file_information
USE output_formats
USE cell_parameters
USE constants

IMPLICIT NONE

! local variables

  integer	:: ff_surface
  integer 	:: output_surface
  integer 	:: output_face
  integer 	:: output_field_number
  integer 	:: number_of_faces
  integer 	:: face
  integer 	:: cz
  integer	:: side
  integer	:: normx,normy,normz
  real*8	:: field(6)
  
  real*8	:: frequency
  complex*16	:: ejwt

! START
  
  CALL write_line('CALLED: face_output_far_field',0,timestepping_output_to_screen_flag)

  do ff_surface=1,n_far_field_surfaces
    
    number_of_faces=far_field_surface(ff_surface)%number_of_faces
    
    frequency=far_field_surface(ff_surface)%frequency
    ejwt=exp(-j*2d0*pi*frequency*time)*dt
    
    do output_face=1,number_of_faces
         
      output_field_number=far_field_surface(ff_surface)%face_output_field_number_list(output_face)
      
      face=far_field_surface(ff_surface)%face_list(output_face)%point
      cz  =far_field_surface(ff_surface)%face_list(output_face)%cell%k
      
      normx=0
      normy=0
      normz=0
      
      if (rank.eq.cell_face_rank(cz,face)) then
! the output point is in this process so set the value  
 	     
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
            
        far_field_surface(ff_surface)%J(output_face,1)= far_field_surface(ff_surface)%J(output_face,1) &
                                    + ejwt*(normy*field(Hz)-normz*field(Hy))
        far_field_surface(ff_surface)%J(output_face,2)= far_field_surface(ff_surface)%J(output_face,2) &
                                    + ejwt*(normz*field(Hx)-normx*field(Hz))
        far_field_surface(ff_surface)%J(output_face,3)= far_field_surface(ff_surface)%J(output_face,3) &
                                    + ejwt*(normx*field(Hy)-normy*field(Hx))
      
        far_field_surface(ff_surface)%M(output_face,1)= far_field_surface(ff_surface)%M(output_face,1) &
                                    - ejwt*(normy*field(Ez)-normz*field(Ey))
        far_field_surface(ff_surface)%M(output_face,2)= far_field_surface(ff_surface)%M(output_face,2) &
                                    - ejwt*(normz*field(Ex)-normx*field(Ez))
        far_field_surface(ff_surface)%M(output_face,3)= far_field_surface(ff_surface)%M(output_face,3) &
                                    - ejwt*(normx*field(Ey)-normy*field(Ex))
	  
      end if ! output face belongs to this process
	
    end do !next cell face in this surface	  
  
  end do ! next far_field_surface
  
  CALL write_line('FINISHED: face_output_far_field',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE face_output_far_field
!
! Name write_far_field_surfaces
!     
!
! Description
!     calculate far E field and output to file
!     3D is based on Balanis, "Advanced Engineering Electromagnetics" Wiley pp285-289
!     2D is based on http://www.eecs.wsu.edu/~schneidj/ufdtd/chap14.pdf
! Comments:
!      
!   Need to check signs for 2D transform - is there a sign discrepancy between 2D and 3D?
!
! History
!
!     started 5/12/2012 CJS adapted from Fieldsolve
!     started to include 2D transform 7/05/2014 CJS 
!

SUBROUTINE write_far_field_surfaces

USE file_information
USE TLM_output
USE mesh
Use TLM_general
Use constants

IMPLICIT NONE

! variables passed to subroutine

! local_variables

  integer	:: ff_surface
  integer 	:: surface
  type(xyz)	:: face_centre
  type(xyz)	:: far_field_point
  real*8  	:: frequency,beta
  real*8  	:: rp_cos_psi
  real*8  	:: r(3),rp(3)
  real*8  	:: r0,theta,phi
  
  real*8 	:: ct,cp,st,sp
  
  real*8 	:: dA
  
  complex*16 	:: Ntheta,Nphi,Ltheta,Lphi
  complex*16 	:: Nz2D,Lz2D
  complex*16 	:: jbeta,prop,prop2
  complex*16 	:: Jsx,Jsy,Jsz,Msx,Msy,Msz
  
  complex*16,allocatable :: Etheta(:),Ephi(:)
  complex*16,allocatable :: Local_Etheta(:),Local_Ephi(:)
  
!  complex*16,allocatable :: Htheta(:),Hphi(:)
!  complex*16,allocatable :: Local_Htheta(:),Local_Hphi(:)
  
  integer 		:: n_data,count
  integer 		:: n_phi_local,phi_loop
  integer 		:: n_theta_local,theta_loop
  
  character*256 temp_filename
  character*256 numbered_filename

! function_types
  
! START

  write(*,*)'CALLED: write_far_field_surfaces'

! Far field transformation...
  
  do ff_surface=1,n_far_field_surfaces
  
    if (far_field_surface(ff_surface)%dim.EQ.2) then
  
      write(*,*)'Do 2D far field transformation'
      dA=dl                                      ! the integral is a line integral for the 2D transformation
      frequency=far_field_surface(ff_surface)%frequency
      beta=2d0*pi*frequency/c0
      jbeta=j*cmplx(beta)
      r0=1d0
      prop2=sqrt(jbeta/(8d0*pi*r0))*exp(jbeta*r0)
      
! loop over phi
     
      n_theta_local=1

      n_phi_local=int ( (far_field_surface(ff_surface)%phi_max-far_field_surface(ff_surface)%phi_min)/   &
    			 far_field_surface(ff_surface)%phi_step)+1

      n_data=n_theta_local*n_phi_local
    
      ALLOCATE( local_Etheta(1:n_data) )
      ALLOCATE( local_Ephi(1:n_data) )
      ALLOCATE( Etheta(1:n_data) )
      ALLOCATE( Ephi(1:n_data) )

      local_Etheta(1:n_data)=(0d0,0d0)
      local_Ephi(1:n_data)  =(0d0,0d0)
    
      Etheta(1:n_data)      =(0d0,0d0)
      Ephi(1:n_data)        =(0d0,0d0)
    
      count=0
 
      theta=pi/2d0

      do phi_loop=1,n_phi_local
        phi= far_field_surface(ff_surface)%phi_min+(phi_loop-1)*far_field_surface(ff_surface)%phi_step
      
        CALL rthetaphi_to_xyz_point(r0,theta,phi,far_field_point)

        r(1)=far_field_point%x
        r(2)=far_field_point%y
        r(3)=far_field_point%z*0d0   ! Note: r(3) is zero here as theta=pi/2 (90 degrees)
    
    	Nz2D=(0d0,0d0)
        Nphi  =(0d0,0d0)
        Lz2D=(0d0,0d0)
        Lphi  =(0d0,0d0)
 
        cp=cos(phi)
        sp=sin(phi)  
    
        do surface=1,far_field_surface(ff_surface)%number_of_faces
        	
          CALL get_cell_point_coordinate(far_field_surface(ff_surface)%face_list(surface),face_centre)
      
          rp(1)=face_centre%x
          rp(2)=face_centre%y
          rp_cos_psi=r(1)*rp(1)+r(2)*rp(2)
        	    
          Jsx=far_field_surface(ff_surface)%J(surface,1)
          Jsy=far_field_surface(ff_surface)%J(surface,2)
          Jsz=far_field_surface(ff_surface)%J(surface,3)
    
          Msx=far_field_surface(ff_surface)%M(surface,1)
          Msy=far_field_surface(ff_surface)%M(surface,2)
          Msz=far_field_surface(ff_surface)%M(surface,3)
    
          prop=exp(jbeta*rp_cos_psi)
    
          Nz2D=Nz2D+Jsz*prop*dA
          Nphi=Nphi+(-Jsx*sp+Jsy*cp)*prop*dA
      
          Lz2D=Lz2D+Msz*prop*dA
          Lphi=Lphi+(-Msx*sp+Msy*cp)*prop*dA
    
        end do ! next surface

        count=count+1

        local_Etheta(count)=-prop2*(Lphi-Z0*Nz2D)    ! this is minus the Ez component
        local_Ephi(count)  = prop2*(Lz2D+Z0*Nphi)	   
        		    
      end do ! next phi
    
#if defined(MPI)
 
      call MPI_REDUCE(local_Etheta, Etheta, &
                      n_data,MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD,ierror)
      call MPI_REDUCE(local_Ephi, Ephi, &
                      n_data,MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD,ierror)

#elif defined(SEQ)

      Etheta=local_Etheta
      Ephi=local_Ephi

#endif
    
      if (rank.eq.0) then
    
! loop over theta and phi
     
        count=0
 
        theta=pi/2d0
      
        n_theta_local=1

        n_phi_local=int ( (far_field_surface(ff_surface)%phi_max-far_field_surface(ff_surface)%phi_min)/   &
                           far_field_surface(ff_surface)%phi_step)+1

        do phi_loop=1,n_phi_local
          phi= far_field_surface(ff_surface)%phi_min+(phi_loop-1)*far_field_surface(ff_surface)%phi_step
	
	  count=count+1
				      
	  write(far_field_output_unit,8020)theta*180.0/pi,phi*180.0/pi,        &
					 abs(Etheta(count)),abs(Ephi(count))
8020      format(4E14.6)

        end do ! next phi
      
      end if  ! (rank.eq.0)
    
      DEALLOCATE( local_Etheta )
      DEALLOCATE( local_Ephi )
      DEALLOCATE( Etheta )
      DEALLOCATE( Ephi )
  
    else 
  
      write(*,*)'Do 3D far field transformation'
   
      dA=dl*dl
      frequency=far_field_surface(ff_surface)%frequency
      beta=2d0*pi*frequency/c0
      jbeta=j*cmplx(beta)
      r0=1d0
      prop2=jbeta*exp(jbeta*r0)/(4d0*pi*r0)
    
! loop over theta and phi

     
      n_theta_local=int ( (far_field_surface(ff_surface)%theta_max-far_field_surface(ff_surface)%theta_min)/   &
      			 far_field_surface(ff_surface)%theta_step)+1

      n_phi_local=int ( (far_field_surface(ff_surface)%phi_max-far_field_surface(ff_surface)%phi_min)/   &
    			 far_field_surface(ff_surface)%phi_step)+1

      n_data=n_theta_local*n_phi_local
    
      ALLOCATE( local_Etheta(1:n_data) )
      ALLOCATE( local_Ephi(1:n_data) )
      ALLOCATE( Etheta(1:n_data) )
      ALLOCATE( Ephi(1:n_data) )

      local_Etheta(1:n_data)=(0d0,0d0)
      local_Ephi(1:n_data)  =(0d0,0d0)
    
      Etheta(1:n_data)      =(0d0,0d0)
      Ephi(1:n_data)        =(0d0,0d0)
    
      count=0
 
      do theta_loop=1,n_theta_local
        theta= far_field_surface(ff_surface)%theta_min+(theta_loop-1)*far_field_surface(ff_surface)%theta_step

        do phi_loop=1,n_phi_local
    	  phi= far_field_surface(ff_surface)%phi_min+(phi_loop-1)*far_field_surface(ff_surface)%phi_step
        
    	  CALL rthetaphi_to_xyz_point(r0,theta,phi,far_field_point)

          r(1)=far_field_point%x
          r(2)=far_field_point%y
          r(3)=far_field_point%z
    
    	  Ntheta=(0d0,0d0)
          Nphi  =(0d0,0d0)
          Ltheta=(0d0,0d0)
          Lphi  =(0d0,0d0)

          ct=cos(theta)
          st=sin(theta)
          cp=cos(phi)
          sp=sin(phi)  
    
    	  do surface=1,far_field_surface(ff_surface)%number_of_faces
        	  
            CALL get_cell_point_coordinate(far_field_surface(ff_surface)%face_list(surface),face_centre)
       
            rp(1)=face_centre%x
            rp(2)=face_centre%y
            rp(3)=face_centre%z
            rp_cos_psi=r(1)*rp(1)+r(2)*rp(2)+r(3)*rp(3)
        	      
    	    Jsx=far_field_surface(ff_surface)%J(surface,1)
    	    Jsy=far_field_surface(ff_surface)%J(surface,2)
    	    Jsz=far_field_surface(ff_surface)%J(surface,3)
    
    	    Msx=far_field_surface(ff_surface)%M(surface,1)
    	    Msy=far_field_surface(ff_surface)%M(surface,2)
    	    Msz=far_field_surface(ff_surface)%M(surface,3)
    
    	    prop=exp(jbeta*rp_cos_psi)
    
    	    Ntheta=Ntheta+(Jsx*ct*cp+Jsy*ct*sp-Jsz*st)*prop*dA
    	    Nphi=Nphi+(-Jsx*sp+Jsy*cp)*prop*dA
        
    	    Ltheta=Ltheta+(Msx*ct*cp+Msy*ct*sp-Msz*st)*prop*dA
    	    Lphi=Lphi+(-Msx*sp+Msy*cp)*prop*dA
    
    	  end do ! next surface
	
	  count=count+1

    	  local_Etheta(count)=-prop2*(Lphi+Z0*Ntheta)
          local_Ephi(count)  = prop2*(Ltheta-Z0*Nphi)        
        		      
        end do ! next phi

      end do ! next theta
    
#if defined(MPI)
 
      call MPI_REDUCE(local_Etheta, Etheta, &
                      n_data,MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD,ierror)
      call MPI_REDUCE(local_Ephi, Ephi, &
                      n_data,MPI_DOUBLE_COMPLEX, MPI_SUM, 0,MPI_COMM_WORLD,ierror)

#elif defined(SEQ)

      Etheta=local_Etheta
      Ephi=local_Ephi

#endif
    
      if (rank.eq.0) then

! open file for far field output

        temp_filename=trim(problem_name)//far_field_output_extn
        CALL add_integer_to_filename(temp_filename,ff_surface,numbered_filename)
  
        OPEN(unit=far_field_output_unit,file=trim(numbered_filename))
        
        write(far_field_output_unit,8010)'# ',n_theta_local,n_phi_local
8010      format(A2,2I6)
    
! loop over theta and phi
     
        count=0
      
        n_theta_local=int ( (far_field_surface(ff_surface)%theta_max-far_field_surface(ff_surface)%theta_min)/   &
                           far_field_surface(ff_surface)%theta_step)+1

        n_phi_local=int ( (far_field_surface(ff_surface)%phi_max-far_field_surface(ff_surface)%phi_min)/   &
                           far_field_surface(ff_surface)%phi_step)+1
 
        do theta_loop=1,n_theta_local
          theta= far_field_surface(ff_surface)%theta_min+(theta_loop-1)*far_field_surface(ff_surface)%theta_step

          do phi_loop=1,n_phi_local
            phi= far_field_surface(ff_surface)%phi_min+(phi_loop-1)*far_field_surface(ff_surface)%phi_step
	  
	    count=count+1
	  				
	    write(far_field_output_unit,8020)theta*180.0/pi,phi*180.0/pi,        &
	                                   abs(Etheta(count)),abs(Ephi(count))

          end do ! next phi

        end do ! next theta
      
      end if  ! (rank.eq.0)
    
      DEALLOCATE( local_Etheta )
      DEALLOCATE( local_Ephi )
      DEALLOCATE( Etheta )
      DEALLOCATE( Ephi )
    
    end if ! 2D or 3D transform
    
  end do ! next far_field_surface

  write(*,*)'FINISHED: write_far_field_surfaces'
  return
  
END SUBROUTINE write_far_field_surfaces

