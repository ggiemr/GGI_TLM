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
! SUBROUTINE generate_far_field_plot
!
! NAME
!    generate_far_field_plot
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/05/2014 CJS
!
!
SUBROUTINE generate_far_field_plot

USE constants
USE post_process
USE file_information

IMPLICIT NONE

! local variables

  character(len=256)	:: filename
  
  character(len=256)	:: command
  character		:: ch
  
  integer	:: op_type
  integer,parameter :: op_type_E_theta=1
  integer,parameter :: op_type_E_phi  =2
  integer,parameter :: op_type_E_total=3
  
  logical	:: file_exists
  
  integer	:: n_theta,n_phi
  
  real*8,allocatable	:: theta(:)
  real*8,allocatable	:: phi(:)
  real*8,allocatable	:: E(:,:)
  real*8,allocatable	:: plot_value(:,:)

  real*8	:: E_theta,E_phi,min_E,max_E
  
  integer	:: theta_loop,phi_loop
 
  character*2   :: ch2
  character*3   :: scale_string
  logical       :: log_scale
  logical       :: normalise
  
  real*8	:: input_real
  real*8 	:: max_data
  real*8 	:: min_data
  real*8 	:: data_range
  
  integer	:: point,quad,n_points,n_quads

  real*8	:: cx,cy,cz
  real*8	:: rmin,rmax
  
  real*8	:: x,y,z,r
  
  real*8        :: Ptot,dtheta,dphi,Pave
  
  character*32 :: field_name

! START

! OPEN FAR FIELD FILE

5 write(*,*)
  write(*,*)'Far field field output files:'
  
  command='ls -ltr *.far_field.fout*'
  CALL system(command)

  write(*,*)'Enter the far field filename'
  read(*,*)filename
  inquire(file=trim(filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  OPEN(unit=local_file_unit,file=filename)
  
  write(*,*)'Enter the output type: 1 for E_theta, 2 for E_phi or 3 for E_total'
  read(*,*)op_type
  if ( (op_type.ne.op_type_E_theta).AND.(op_type.ne.op_type_E_phi).AND.(op_type.ne.op_type_E_total) ) then
    write(*,*)'Error: output type should be 1 for E_theta, 2 for E_phi or 3 for E_total'
    STOP
  end if

  write(record_user_inputs_unit,*)op_type,' Output type'
  
  if (op_type.EQ.op_type_E_theta) then
    field_name='Far_field_E_theta'
  else if (op_type.EQ.op_type_E_phi) then
    field_name='Far_field_E_phi'
  else
    field_name='Far_field_E_total'
  end if

! READ THE FAR FIELD FILE

  read(local_file_unit,*)ch2,n_theta,n_phi

  allocate( theta(1:n_theta) )
  allocate( phi(1:n_phi) )
  allocate( E(1:n_theta,1:n_phi) )
  allocate( plot_value(1:n_theta,1:n_phi) )
 
! loop over theta and phi
  min_E=1d30
  max_E=0d0

  do theta_loop=1,n_theta

    do phi_loop=1,n_phi

      read(local_file_unit,*)theta(theta_loop),phi(phi_loop),E_theta,E_phi
      
      if (op_type.EQ.op_type_E_theta) then
      
        E(theta_loop,phi_loop)=E_theta
	
      else if (op_type.EQ.op_type_E_phi) then
      
        E(theta_loop,phi_loop)=E_phi
      
      else
      
        E(theta_loop,phi_loop)=sqrt(E_theta*E_theta+E_phi*E_phi)
	
      end if
      
      min_E=min(min_E,E(theta_loop,phi_loop))
      max_E=max(max_E,E(theta_loop,phi_loop))
      
    end do
    
  end do
  
  close(unit=local_file_unit)
  
  write(*,*)'Maximum far field value read from file :',max_E
  write(*,*)'Minimum far field value read from file :',min_E
  
! Integrate the total radiated power so that we can normalise the field to give dBi.

  Ptot=0d0
  
  do theta_loop=1,n_theta-1

    do phi_loop=1,n_phi-1
    
      dtheta=(theta(theta_loop+1)-theta(theta_loop))*pi/180d0
      dphi=(phi(phi_loop+1)-phi(phi_loop))*pi/180d0
      
      Pave=(E(theta_loop,phi_loop)**2+E(theta_loop+1,phi_loop)**2+    &
            E(theta_loop,phi_loop+1)**2+E(theta_loop+1,phi_loop+1)**2)/4d0
      
      Ptot=Ptot+Pave*dtheta*dphi*sin( (theta(theta_loop+1)+theta(theta_loop))*pi/360d0 )
      
    end do
    
  end do
    
  Ptot=Ptot/(4d0*pi)
  
  write(*,*)'Ptot=',Ptot
  
  if (Ptot.NE.0d0) then 
    E(:,:)=E(:,:)/sqrt(Ptot)   
  end if   
  
  min_E=1d30
  max_E=0d0
  
  do theta_loop=1,n_theta
    do phi_loop=1,n_phi      
      min_E=min(min_E,E(theta_loop,phi_loop))
      max_E=max(max_E,E(theta_loop,phi_loop))      
    end do   
  end do
  
  write(*,*)'Maximum far field value scaled to isotropic radiator :',max_E
  write(*,*)'Minimum far field value scaled to isotropic radiator :',min_E
  
  write(*,*)'Enter scale for far field (log or lin)'
  read(*,'(A3)')scale_string
  write(record_user_inputs_unit,'(A3)')trim(scale_string)
  CALL convert_to_lower_case(scale_string,3)
  if (scale_string(1:3).eq.'log') then
    log_scale=.TRUE.
  else if (scale_string(1:3).eq.'lin') then
    log_scale=.FALSE.
  else
    write(*,*)'Scale should be log or lin'
    stop
  end if
  
  write(*,*)'Do you want to normalise the far field to a maximum value of 1V/m (0dB)?  (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A)')ch
  CALL convert_to_lower_case(scale_string,3)
  if (ch.eq.'y') then
    normalise=.TRUE.
    if (scale_string(1:3).eq.'log') then
      field_name=trim(field_name)//'(dB)'
    else if (scale_string(1:3).eq.'lin') then
      field_name=trim(field_name)//'(V/m)'
    end if
  else
    normalise=.FALSE.
    if (scale_string(1:3).eq.'log') then
      field_name=trim(field_name)//'(dBi)'
    else if (scale_string(1:3).eq.'lin') then
      field_name=trim(field_name)//'(V/m)'
    end if
  end if

! normalise and logscale if required
  do theta_loop=1,n_theta

    do phi_loop=1,n_phi
    
      if (normalise) then
        E(theta_loop,phi_loop)=E(theta_loop,phi_loop)/max_E
      end if
      
      if (log_scale) then
        E(theta_loop,phi_loop)=20d0*log10(E(theta_loop,phi_loop))
      end if
        
    end do
  
  end do
  
  if (normalise) then
    min_E=min_E/max_E
    max_E=1d0
  end if
  
  if (log_scale) then
    max_E=20d0*log10(max_E)
    min_E=20d0*log10(min_E)
  end if
 
  write(*,*)'Maximum E field value:',max_E
  write(*,*)'Minimum E field value:',min_E
  
  write(*,*)'Enter Maximum value to use in plots or 0 to use maximum in data'
  read(*,*)input_real
  write(record_user_inputs_unit,*)input_real
  if (input_real.ne.0D0) then
    max_E=input_real
  end if
  
  write(*,*)'Enter Minimum value to use in plots or 0 to use minimum in data'
  read(*,*)input_real
  write(record_user_inputs_unit,*)input_real
  if (input_real.ne.0D0) then
    min_E=input_real
  end if
  
  if ((max_E-min_E).le.0D0) then
    write(*,*)'Maximum data value should be greater than the minimum data value'
    STOP
  end if
  
  write(*,*)'Number of theta values read:',n_theta
  write(*,*)'Number of phi values read  :',n_phi
  
  write(*,*)'Maximum field value for plotting:',max_E
  write(*,*)'Minimum field value for plotting:',min_E
   
  data_range=(max_E-min_E)
  if (data_range.eq.0d0) data_range=1d0
  
  write(*,*)'Enter the cente coordinates for the far field surface'
  read(*,*)cx,cy,cz
  write(record_user_inputs_unit,*)cx,cy,cz,' cente coordinates for the far field surface'
  
  write(*,*)'Enter the max radius for the far field surface'
  read(*,*)rmax
  write(record_user_inputs_unit,*)rmax,' max radius for the far field surface'

  write(*,*)'Enter the far field plot name (without .vtk extension)'
  read(*,*)filename
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  OPEN(unit=local_file_unit,file=trim(filename)//'.vtk')
  
  n_quads=(n_theta-1)*(n_phi-1)
  n_points=n_quads*4
    
! WRITE HEADER INFORMATION    
! write vtk header      
    write(local_file_unit,'(A)')'# vtk DataFile Version 2.0'
    write(local_file_unit,'(A)')trim(filename)
    write(local_file_unit,'(A)')'ASCII'
    write(local_file_unit,'(A)')'DATASET POLYDATA'
    write(local_file_unit,'(A,I10,A)')'POINTS',n_points,' float'
    
! WRITE POINT DATA 
    point=0

    do theta_loop=1,n_theta-1

      do phi_loop=1,n_phi-1
    
! write 4 points for each quad    
        point=point+1
        CALL get_xyz_far_field(E(theta_loop  ,phi_loop  ),theta(theta_loop  ),phi(phi_loop  ),min_E,max_E,	&
	                       x,y,z,cx,cy,cz,rmax,plot_value(theta_loop  ,phi_loop  ),log_scale)
        write(local_file_unit,8005)x,y,z
      
        point=point+1
        CALL get_xyz_far_field(E(theta_loop+1,phi_loop  ),theta(theta_loop+1),phi(phi_loop  ),min_E,max_E,	&
	                       x,y,z,cx,cy,cz,rmax,plot_value(theta_loop+1,phi_loop  ),log_scale)
        write(local_file_unit,8005)x,y,z
      
        point=point+1
        CALL get_xyz_far_field(E(theta_loop+1,phi_loop+1),theta(theta_loop+1),phi(phi_loop+1),min_E,max_E,	&
	                       x,y,z,cx,cy,cz,rmax,plot_value(theta_loop+1,phi_loop+1),log_scale)
        write(local_file_unit,8005)x,y,z
      
        point=point+1
        CALL get_xyz_far_field(E(theta_loop  ,phi_loop+1),theta(theta_loop  ),phi(phi_loop+1),min_E,max_E,	&
	                       x,y,z,cx,cy,cz,rmax,plot_value(theta_loop  ,phi_loop+1),log_scale)
        write(local_file_unit,8005)x,y,z
     
8005  format(3E14.5)
    	  
      end do ! next phi_loop
    	  
    end do ! next theta_loop
    
    if (point.ne.n_points) then
      write(*,*)'Error in generate_far_field_plot'
      write(*,*)'Problem with the number of output points, n_points=',n_points
      write(*,*)'Number of coordinates written=',point
      STOP
    end if
  
! write quad data
    write(local_file_unit,'(A,2I10)')'POLYGONS',n_quads,n_quads*5

    point=0
    do quad=1,n_quads
    
      write(local_file_unit,8010)4,point,point+1,point+2,point+3
      point=point+4

8010  format(I3,4I8)
      
    end do ! next quad
  
! WRITE DATA ASSOCIATED WITH POINTS

! write point based data
    write(local_file_unit,'(A,I10)')'POINT_DATA ',n_points
    write(local_file_unit,'(A,A,A)')'SCALARS ',trim(field_name),' float 1'
    write(local_file_unit,'(A)')'LOOKUP_TABLE field_on_cells_table'

    point=0

    do theta_loop=1,n_theta-1

      do phi_loop=1,n_phi-1
    
! write 4 points for each quad    
        write(local_file_unit,8020)plot_value(theta_loop  ,phi_loop  )
        write(local_file_unit,8020)plot_value(theta_loop+1,phi_loop  )
        write(local_file_unit,8020)plot_value(theta_loop+1,phi_loop+1)
        write(local_file_unit,8020)plot_value(theta_loop  ,phi_loop+1)
      
8020  format(E14.5)
    	  
    end do ! next phi
    	  
  end do ! next theta
  
  CLOSE(unit=local_file_unit)
  
  deallocate( theta )
  deallocate( phi )
  deallocate( E ) 
  deallocate( plot_value ) 
  
  RETURN
  
9000 continue
  write(*,*)'Error opening output file'
  write(*,*)'Filename='
  write(*,*)trim(filename)
  stop

  RETURN
  
END SUBROUTINE generate_far_field_plot

!
! ___________________________________________________________________
!
!
  SUBROUTINE get_xyz_far_field(E,theta,phi,min_E,max_E,x,y,z,cx,cy,cz,rmax,plot_value,log_scale)
  
USE constants  
  
  IMPLICIT NONE
  
  real*8	:: E,theta,phi,min_E,max_E,x,y,z,cx,cy,cz,rmax,plot_value
  logical	:: log_scale
  
! local variables
  real*8	:: r,t,p

! START

  plot_value=E

! keep data in range...
  
  if (plot_value.gt.max_E) plot_value=max_E
  if (plot_value.lt.min_E) plot_value=min_E
  
  if (abs(plot_value).lt.1D-30) plot_value=0d0
  
  if (log_scale) then
! radius is zero at min_E and rmax at max_E
    r=rmax*(plot_value-min_E)/(max_E-min_E)
  else
! on the linear scale the radius is zero for zero field and rmax at maximum field
    r=rmax*plot_value/max_E
  end if
  
  t=theta*pi/180d0
  p=phi*pi/180d0
  
  x=r*sin(t)*cos(p)+cx
  y=r*sin(t)*sin(p)+cy
  z=r*cos(t)       +cz

  END SUBROUTINE get_xyz_far_field
