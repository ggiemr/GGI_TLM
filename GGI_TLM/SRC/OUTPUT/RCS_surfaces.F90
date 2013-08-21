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
! SUBROUTINE set_rcs_surfaces_in_mesh
! SUBROUTINE initialise_rcs_surfaces
! SUBROUTINE face_output_rcs
!
! NAME
!     set_rcs_surfaces_in_mesh
!
! DESCRIPTION
!
!     
! COMMENTS
!
! HISTORY
!
!     started 2/1/2013 CJS
!
!
SUBROUTINE set_rcs_surfaces_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE file_information

IMPLICIT NONE

! local variables

  integer	:: surface_number
  integer	:: number_of_faces
  integer	:: output_face
  integer	:: cx,cy,cz,face
  type(cell_point)	:: output_face1
  type(cell_point)	:: output_face2

! START
  
  CALL write_line('CALLED: set_rcs_surfaces_in_mesh',0,output_to_screen_flag)
  
  if (n_rcs_surfaces.ne.0) then
  
    surface_number=rcs_surface%surface_number
    number_of_faces=problem_surfaces(surface_number)%number_of_faces
    rcs_surface%number_of_faces=number_of_faces
      
! allocate face list for the output surface      
    ALLOCATE( rcs_surface%face_list(1:number_of_faces) )
      
    do output_face=1,number_of_faces
      
! copy the faces from the geometry structure to the local
! frequency_output_surface structure
	
      cx=problem_surfaces(surface_number)%face_list(output_face)%cell%i
      cy=problem_surfaces(surface_number)%face_list(output_face)%cell%j
      cz=problem_surfaces(surface_number)%face_list(output_face)%cell%k
      face=problem_surfaces(surface_number)%face_list(output_face)%point

! check the side of the surface on which we want output and change if required,

      if      (face.eq.face_xmin) then      
        if (.NOT.rcs_surface%output_on_outward_normal) then ! output on other side of face
          cx=cx-1
          face=face_xmax
        end if
       else if (face.eq.face_xmax) then        
        if (.NOT.rcs_surface%output_on_outward_normal) then ! output on other side of face
          cx=cx+1
          face=face_xmin
        end if
      else if (face.eq.face_ymin) then         
        if (.NOT.rcs_surface%output_on_outward_normal) then ! output on other side of face
          cy=cy-1
          face=face_ymax
        end if
      else if (face.eq.face_ymax) then         
        if (.NOT.rcs_surface%output_on_outward_normal) then ! output on other side of face
          cy=cy+1
          face=face_ymin
        end if
      else if (face.eq.face_zmin) then         
        if (.NOT.rcs_surface%output_on_outward_normal) then ! output on other side of face
          cz=cz-1
          face=face_zmax
        end if
      else if (face.eq.face_zmax) then         
        if (.NOT.rcs_surface%output_on_outward_normal) then ! output on other side of face
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

      end if ! output face belongs to this processor

! preserve the face which includes the normal information in the rcs_surface%face_list	
      rcs_surface%face_list(output_face)=output_face1
	       
    end do !next cell face in this surface	

  end if ! n_rcs_surfaces.ne.0
  
  CALL write_line('FINISHED: set_rcs_surfaces_in_mesh',0,output_to_screen_flag)
  
END SUBROUTINE set_rcs_surfaces_in_mesh
!
! NAME
!     initialise_rcs_surfaces
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 2/1/2013  CJS
!
!
SUBROUTINE initialise_rcs_surfaces

USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE mesh
USE TLM_output
USE file_information
USE constants

IMPLICIT NONE

! local variables

  integer	:: cx,cy,cz,face
  integer	:: number_of_faces
  integer	:: output_face
  
  real*8 	:: r,theta,phi,kx,ky,kz
  real*8 	:: dxdt,dxdp,dydt,dydp,dzdt,dzdp,vx,vy,vz,norm

  type(xyz)	:: xyz_point
  type(xyz)	:: K
  
  real*8	:: local_dmin,local_dmax
  real*8	:: dmin,dmax,dmrcs
  integer	:: n_data
  
  real*8	:: far_field_time_period
  real*8	:: far_field_distance

! function variables

  real*8	:: xyz_dot
  
! START
  
  if (n_rcs_surfaces.ne.0) then
  
    if (rank.eq.0) write(info_file_unit,*)'____________________________________________________'
    if (rank.eq.0) write(info_file_unit,*)''
    if (rank.eq.0) write(info_file_unit,*)'RCS surface'
    if (rank.eq.0) write(info_file_unit,*)''
  
! calculate K vector in RCS direction
    r=1d0
    theta=rcs_surface%theta
    phi=rcs_surface%phi
    
    CALL rthetaphi_to_xyz_point(r,theta,phi,K)
    
    rcs_surface%Kmrcs(1)=K%x
    rcs_surface%Kmrcs(2)=K%y
    rcs_surface%Kmrcs(3)=K%z
  
! calculate vectors perpendicular to K in the theta and phi directions
    if ((theta.eq.0d0).OR.(theta.eq.180d0)) then
      rcs_surface%Vtheta(1)=1d0
      rcs_surface%Vtheta(2)=0d0
      rcs_surface%Vtheta(3)=0d0
      rcs_surface%Vphi(1)=0d0
      rcs_surface%Vphi(2)=1d0
      rcs_surface%Vphi(3)=0d0
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
8000    format(A,3F12.6)
        STOP
      end if
      rcs_surface%Vtheta(1)=vx/norm
      rcs_surface%Vtheta(2)=vy/norm
      rcs_surface%Vtheta(3)=vz/norm

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
      rcs_surface%Vphi(1)=vx/norm
      rcs_surface%Vphi(2)=vy/norm
      rcs_surface%Vphi(3)=vz/norm

    end if
    
    number_of_faces=rcs_surface%number_of_faces
    
    if (rank.eq.0) write(info_file_unit,*)'Number of RCS surface faces=',number_of_faces
    
    ALLOCATE( rcs_surface%face_output_field_number_list(1:number_of_faces) )
    
    rcs_surface%face_output_field_number_list(1:number_of_faces)=0
        
    local_dmin=1d30
    local_dmax=-1d30

    do output_face=1,number_of_faces
   
      cx  =rcs_surface%face_list(output_face)%cell%i
      cy  =rcs_surface%face_list(output_face)%cell%j
      cz  =rcs_surface%face_list(output_face)%cell%k
      face=rcs_surface%face_list(output_face)%point 

      if (rank.eq.cell_face_rank(cz,face)) then
           
    	if	(face.eq.face_xmin) then      
    	  rcs_surface%face_output_field_number_list(output_face)= &
        		      local_surface_output(cx  ,cy  ,cz  ,face_xmin) 
    	else if (face.eq.face_xmax) then       
    	  rcs_surface%face_output_field_number_list(output_face)= &
        		      local_surface_output(cx+1,cy  ,cz  ,face_xmin)
    	else if (face.eq.face_ymin) then       
    	  rcs_surface%face_output_field_number_list(output_face)= &
        		      local_surface_output(cx  ,cy  ,cz  ,face_ymin)
    	else if (face.eq.face_ymax) then       
    	  rcs_surface%face_output_field_number_list(output_face)= &
        		      local_surface_output(cx  ,cy+1,cz  ,face_ymin)
    	else if (face.eq.face_zmin) then       
    	  rcs_surface%face_output_field_number_list(output_face)= &
        		      local_surface_output(cx  ,cy  ,cz  ,face_zmin)
    	else if (face.eq.face_zmax) then       
    	  rcs_surface%face_output_field_number_list(output_face)= &
        		      local_surface_output(cx  ,cy  ,cz+1,face_zmin)
    	end if

! work out the distance of this point along the K vector in the far field direction	    
	CALL get_cell_point_coordinate(rcs_surface%face_list(output_face),xyz_point)
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
    rcs_surface%n_far_field_points=n_rcs_points_per_cell*far_field_distance/dl+2
    rcs_surface%far_field_point_offset=n_rcs_points_per_cell*dmin/dl-1
    
    if (rank.eq.0) write(info_file_unit,*)'Far field time period=',far_field_time_period
    if (rank.eq.0) write(info_file_unit,*)'Far field distance	=',far_field_distance
    if (rank.eq.0) write(info_file_unit,*)'n_far_field_points=',rcs_surface%n_far_field_points
    if (rank.eq.0) write(info_file_unit,*)'far_field_point_offset=',rcs_surface%far_field_point_offset
    if (rank.eq.0) write(info_file_unit,*)''
    
    ALLOCATE( rcs_surface%Etheta(1:rcs_surface%n_far_field_points) )
    ALLOCATE( rcs_surface%Htheta(1:rcs_surface%n_far_field_points) )
    ALLOCATE( rcs_surface%Ephi(1:rcs_surface%n_far_field_points) )
    ALLOCATE( rcs_surface%Hphi(1:rcs_surface%n_far_field_points) )
    
    rcs_surface%Etheta(1:rcs_surface%n_far_field_points)=0d0
    rcs_surface%Htheta(1:rcs_surface%n_far_field_points)=0d0
    rcs_surface%Ephi(1:rcs_surface%n_far_field_points)  =0d0
    rcs_surface%Hphi(1:rcs_surface%n_far_field_points)  =0d0

    if (rank.eq.0) then
! rank 0 process only: open file for rcs output
  
      OPEN(unit=trcs_output_unit,file=trim(problem_name)//trcs_output_extn)
      OPEN(unit=rcs_output_unit,file=trim(problem_name)//rcs_output_extn)

    end if ! rank 0 process

  end if ! n_rcs_surfaces.ne.0

  RETURN

END SUBROUTINE initialise_rcs_surfaces

!
! SUBROUTINE face_output_rcs
!
! NAME
!     face_output_rcs
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 2/1/2013 CJS
!
!
SUBROUTINE face_output_rcs

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

! function variables

  real*8	:: xyz_dot
 
! START
  
  CALL write_line('CALLED: face_output_rcs',0,timestepping_output_to_screen_flag)
  
  if (n_rcs_surfaces.ne.0) then
  
    dA=dl*dl
   
    theta=rcs_surface%theta
    phi=rcs_surface%phi
    
    ct=cos(theta)
    st=sin(theta)
    cp=cos(phi)
    sp=sin(phi)  
    
    K%x=rcs_surface%Kmrcs(1)
    K%y=rcs_surface%Kmrcs(2)
    K%z=rcs_surface%Kmrcs(3)
    
    number_of_faces=rcs_surface%number_of_faces
    
    do output_face=1,number_of_faces
    
      output_field_number=rcs_surface%face_output_field_number_list(output_face)
               
      face=rcs_surface%face_list(output_face)%point
      cz  =rcs_surface%face_list(output_face)%cell%k
      
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
 
	CALL get_cell_point_coordinate(rcs_surface%face_list(output_face),xyz_point)
	dmrcs=xyz_dot(K,xyz_point)

! add the offset distance due to the time	    
        dmrcs=-dmrcs+time*c0
	    
! calculate the point to update in the far field array
        point=INT(dble(n_rcs_points_per_cell)*dmrcs/dl)	   
	       
! add the offset distance due to the geometry of the surface	    
        point=point-rcs_surface%far_field_point_offset
	    
        if ( (point.ge.1).AND.(point.le.rcs_surface%n_far_field_points) ) then
! calculate far field and add contributions to this point
	
! Note:include the time derivative required for the far field calculation 
! in the frequency domain later on (output.F90)     
	
          Ntheta=( Jsx*ct*cp+Jsy*ct*sp-Jsz*st)*dA/(4d0*pi)
          Nphi  =(-Jsx*sp   +Jsy*cp	  )*dA/(4d0*pi)
	 
          Ltheta=( Msx*ct*cp+Msy*ct*sp-Msz*st)*dA/(4d0*pi)
          Lphi  =(-Msx*sp   +Msy*cp	  )*dA/(4d0*pi)
	
          Etheta=-(Lphi+Z0*Ntheta)
	  Ephi  = (Ltheta-Z0*Nphi)
	
	  Htheta= (Nphi-Ltheta/Z0)
	  Hphi  =-(Ntheta+Lphi/Z0)
	    
	  
          rcs_surface%Etheta(point)=rcs_surface%Etheta(point)+Etheta
          rcs_surface%Ephi(point)  =rcs_surface%Ephi(point)  +Ephi
          rcs_surface%Htheta(point)=rcs_surface%Htheta(point)+Htheta
          rcs_surface%Hphi(point)  =rcs_surface%Hphi(point)  +Hphi
	
        else
      
          write(*,*)'Trying to put far field point outside array',point,rcs_surface%n_far_field_points

        end if
	 
      end if ! in this processor's mesh
      
    end do
    
  end if ! n_rcs_surfaces.ne.0

  
  CALL write_line('FINISHED: face_output_rcs',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE face_output_rcs
!
! Name write_rcs_surfaces
!     
!
! Description
!     calculate time domain far field and RCS then output to file
!
! Comments:
!      
!
! History
!
!     started 2/1/2013  CJS adapted from Fieldsolve
!

SUBROUTINE write_rcs_surfaces

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

  integer		:: n_frequencies,frequency_loop
  real*8		:: frequency
  integer		:: function_number
  
  complex*16		:: Einc
  real*8		:: mod_Einc
  
  complex*16		:: Fourier_const
  complex*16 		:: Etheta,Ephi,Htheta,Hphi,Etot,Htot
  
  real*8		:: RCS
  
! function_types
  
! START

  write(*,*)'CALLED: write_rcs_surfaces'

! RCS calculation...
  
  if (n_rcs_surfaces.ne.0) then
    
    n_data=rcs_surface%n_far_field_points
    
    allocate( local_Etheta(1:n_data) )
    allocate( local_Ephi(1:n_data) )
    allocate( local_Htheta(1:n_data) )
    allocate( local_Hphi(1:n_data) )
 
  
#if defined(MPI)

    call MPI_REDUCE(rcs_surface%Etheta, local_Etheta, &
                    n_data,MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(rcs_surface%Ephi, local_Ephi, &
                    n_data,MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(rcs_surface%Htheta, local_Htheta, &
                    n_data,MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(rcs_surface%Hphi, local_Hphi, &
                    n_data,MPI_DOUBLE_PRECISION, MPI_SUM, 0,MPI_COMM_WORLD,ierror)
		    
#elif defined(SEQ)

    local_Etheta=rcs_surface%Etheta
    local_Ephi  =rcs_surface%Ephi
    local_Htheta=rcs_surface%Htheta
    local_Hphi  =rcs_surface%Hphi

#endif

    if (rank.eq.0) then

! write down range response to file    
      dmin=1
      dmax=rcs_surface%n_far_field_points
      do d_local=dmin,dmax 
    
        t_local=d_local*dl/(n_rcs_points_per_cell*c0)    

        REtheta=local_Etheta(d_local)
        REphi=local_Ephi(d_local)
        RHtheta=local_Htheta(d_local)
        RHphi=local_Hphi(d_local)
      
        if (abs(REtheta).lt.1D-30) REtheta=0d0
        if (abs(REphi).lt.1D-30)   REphi=0d0
        if (abs(RHtheta).lt.1D-30) RHtheta=0d0
        if (abs(RHphi).lt.1D-30)   RHphi=0d0
      
        write(trcs_output_unit,8020)d_local,t_local,REtheta,REphi,RHtheta,RHphi
     
8020    format(I10,5E14.6)
      end do ! next timestep
   
      dA=dl*dl
      dt_local=dl/(n_rcs_points_per_cell*c0)    
    
      n_frequencies=int( (rcs_surface%fmax-rcs_surface%fmin)/rcs_surface%fstep )+1

      do frequency_loop=1,n_frequencies

        frequency=rcs_surface%fmin+(frequency_loop-1)*rcs_surface%fstep
    
! calculate the incident field at the frequency of interest for RCS calculation
        function_number=huygens_surface%excitation_function_number
        Einc=(0d0,0d0)
  
        do timestep_loop=1,n_timesteps

          t_local=(timestep_loop-1)*dt

! evaluate the function at the current time
          Fourier_const=dt*exp(-j*2d0*pi*frequency*t_local)
	  Einc=Einc+excitation_functions(function_number)%value(timestep_loop)*Fourier_const
	  
        end do ! next timestep
	
        mod_Einc=abs(Einc)     
     
! Fourier transform the far field

        Etheta=(0d0,0d0)
        Ephi=(0d0,0d0)
        Htheta=(0d0,0d0)
        Hphi=(0d0,0d0)
        dmin=1
        dmax=rcs_surface%n_far_field_points 
	
        do d_local=dmin,dmax 
	
          t_local=d_local*dl/(n_rcs_points_per_cell*c0)
! evaluate the function at the current time

          Fourier_const=dt_local*exp(-j*2d0*pi*frequency*t_local)
	  
	  Etheta=Etheta+local_Etheta(d_local)*Fourier_const
	  Ephi=Ephi+local_Ephi(d_local)*Fourier_const
	  Htheta=Htheta+local_Htheta(d_local)*Fourier_const
	  Hphi=Hphi+local_Hphi(d_local)*Fourier_const
	  
        end do ! next timestep
            
! include the time derivative required for the far field calculation in the frequency domain      
        Etheta=Etheta*j*2d0*pi*frequency/c0
        Ephi  =Ephi*  j*2d0*pi*frequency/c0
        Htheta=Htheta*j*2d0*pi*frequency/c0
        Hphi  =Hphi*  j*2d0*pi*frequency/c0
      
        Cpower=Etheta*conjg(Hphi)-Ephi*conjg(Htheta)
        Rpower=dble(Cpower)
        Ipower=Mod_Einc**2/Z0
              
        if ((Rpower.gt.0d0).AND.(Ipower.ne.0d0)) then
          RCS=10d0*log10(4d0*pi*Rpower/Ipower)
        else
          RCS=-1000d0
        end if

        write(rcs_output_unit,8040)frequency,abs(Etheta),abs(Ephi),abs(Einc),RCS
8040    format(5E14.6)
	
      end do ! next frequency
    
    end if   !(rank.eq.0)
    
    deallocate( local_Etheta )
    deallocate( local_Ephi )
    deallocate( local_Htheta )
    deallocate( local_Hphi )
    
  end if ! n_rcs_surfaces.ne.0


  write(*,*)'FINISHED: write_rcs_surfaces'
  return
  
END SUBROUTINE write_rcs_surfaces

