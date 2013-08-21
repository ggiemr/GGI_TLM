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
! SUBROUTINE create_frequency_domain_animation
!
! NAME
!    create_frequency_domain_animation
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 28/01/2013 CJS
!
!
SUBROUTINE create_frequency_domain_animation

USE post_process
USE file_information
USE constants

IMPLICIT NONE

! local variables

  character(len=256)	:: filename
  
  integer logscale
  parameter (logscale=1)
  
  integer linscale
  parameter (linscale=2)
  
  character*256 base_filename
  character*256 temp_filename
  character*256 temp_filename2
  character*256 temp_filename3
  character*256 temp_filename4
  character*256 frame_filename
  
  integer n_surfaces,surface_loop,surface
  integer output_surface
  
  TYPE(surface_animation_data),allocatable	:: surface_animation(:)
  
  integer			:: n_points,n_quads,n_frames
  
  real*8 			:: max_data
  real*8 			:: min_data
  real*8 			:: data_range
  
  integer max_point_number

  integer max_quad_number
  
  character*256 :: ipline

  character(len=256)	:: command
  
  logical	:: file_exists
  
  integer 	:: ip_point,ip_quad
  integer 	:: point,quad
  integer 	:: point1,point2,point3,point4
  real*8 		:: x,y,z
  
  integer 	:: frame
  integer 	:: i_max_data,i_min_data,i_data
  
  integer 	:: i_red,i_green,i_blue
  
  integer i
  character*3 	:: scale_string
  integer     	:: scale
  real*8     	:: input_real
  real*8     	:: point_data
  
  character 	:: ch
  logical 	:: scale_data
  
  integer 	:: cx,cy,cz
  real*8	:: re,im,mag,phase

! START

! STAGE 1. Open data file

5 write(*,*)
  write(*,*)'Surface field output files:'
  
  command='ls -ltr *.frequency_output_surface.fout'
  CALL system(command)

  write(*,*)'Enter the surface field filename'
  read(*,*)filename
  inquire(file=trim(filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  OPEN(unit=local_file_unit,file=filename)
   
  base_filename=trim(filename)
    
! STAGE 2. Read frequency_output_surface file

!# NUMBER OF FREQUENCY OUTPUT SURFACES:
!2
!# START SURFACE FIELD FILE TEMPLATE, OUTPUT SURFACE NUMBER:    1
!# Number of points:'
!20
!# Number of faces:'
!80
!# Number of frames:'
!24

  write(*,*)'Reading frequency_output_surface file:'
  
  read(local_file_unit,'(A80)'),ipline
  read(local_file_unit,*)n_surfaces
  
  ALLOCATE ( surface_animation(1:n_surfaces) )
#
! read header information  
  do surface=1,n_surfaces
  
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,*)n_points
    surface_animation(surface)%n_points=n_points
     
    write(*,*)'Number of points=',surface_animation(surface)%n_points
  
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,*)n_quads
    surface_animation(surface)%n_quads=n_quads
     
    write(*,*)'Number of quads=',surface_animation(surface)%n_quads
  
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,*)n_frames
    surface_animation(surface)%n_frames=n_frames
     
    write(*,*)'Number of frames=',surface_animation(surface)%n_frames

    write(*,*)'Allocating data'
       
    ALLOCATE ( surface_animation(surface)%points(1:n_points,1:3) ) 
    ALLOCATE ( surface_animation(surface)%quads(1:n_quads,1:4) ) 
    ALLOCATE ( surface_animation(surface)%complex_data(1:n_quads) ) 
  
    surface_animation(surface)%points(1:n_points,1:3)=0D0
    surface_animation(surface)%quads(1:n_quads,1:4)=0
    surface_animation(surface)%complex_data(1:n_quads)=0D0
  
    write(*,*)'Reading points'

    do point=1,n_points
      read(local_file_unit,*)x,y,z
      surface_animation(surface)%points(point,1)=x
      surface_animation(surface)%points(point,2)=y
      surface_animation(surface)%points(point,3)=z
    end do
    
! read quads    
    write(*,*)'Reading quads'

    do quad=1,n_quads
      read(local_file_unit,*)ip_quad,point1,point2,point3,point4
      surface_animation(surface)%quads(quad,1)=point1
      surface_animation(surface)%quads(quad,2)=point2
      surface_animation(surface)%quads(quad,3)=point3
      surface_animation(surface)%quads(quad,4)=point4
    end do
    
    surface_animation(surface)%max_data=-1e30
    surface_animation(surface)%min_data=1e30
    
! read data    
    write(*,*)'Reading data'
    do quad=1,n_quads
    
      read(local_file_unit,*)cx,cy,cz,re,im,mag
      
      surface_animation(surface)%complex_data(quad)=cmplx(re,im)
      
      surface_animation(surface)%max_data=max(surface_animation(surface)%max_data,mag)
      surface_animation(surface)%min_data=min(surface_animation(surface)%min_data,mag)
      
    end do
    
    write(*,*)'Finished reading fieldsolve geometry data'
      
  end do ! read next surface
 

100 CONTINUE  ! start of surface animation output loop

  write(*,*)' '
  write(*,*)'Number of surfaces read:',n_surfaces
  write(*,*)' '
  
! read the number of the surface to write to file  

  write(*,*)'Enter the number of the surface to output field animation data or 0 to finish'
  read(*,*)output_surface
  
  if (output_surface.eq.0) then
    
    write(record_user_inputs_unit,*)output_surface,' Output surface'
    do surface=1,n_surfaces
      if ( allocated( surface_animation(surface)%points ) ) DEALLOCATE ( surface_animation(surface)%points ) 
      if ( allocated( surface_animation(surface)%quads ) ) DEALLOCATE ( surface_animation(surface)%quads ) 
      if ( allocated( surface_animation(surface)%complex_data ) ) DEALLOCATE ( surface_animation(surface)%complex_data ) 
    end do
      
    if ( allocated( surface_animation ) ) DEALLOCATE ( surface_animation )
    
    RETURN
  
  end if ! return to main post_process 
  
  if ( (output_surface.lt.0).OR.(output_surface.gt.n_surfaces) )then
    write(*,*)'Output surface must be between 0 and ',n_surfaces
    GOTO 100
  end if

  write(*,*)'Output field animation on surface number',output_surface
  write(record_user_inputs_unit,*)output_surface,' Output surface'

! STAGE 4. Work out scaling values for data
  max_data=surface_animation(output_surface)%max_data
  min_data=surface_animation(output_surface)%min_data
  write(*,*)'Maximum data value:',max_data
  write(*,*)'Minimum data value:',min_data
      
  if (scale.eq.logscale) then
    min_data=log(max_data/1000.0)
    max_data=log(max_data)
    write(*,*)'Maximum data value for plotting log scale:',max_data
    write(*,*)'Minimum data value for plotting log scale:',min_data
  end if
  
  data_range=(max_data-min_data)
  if (data_range.eq.0d0) data_range=1d0
  
! add surface number to base filename

  temp_filename=trim(base_filename)//'_surface'

  CALL add_integer_to_filename(temp_filename,output_surface,temp_filename2)
  temp_filename3=trim(temp_filename2)//'_frame'
    
  n_points=surface_animation(output_surface)%n_points
  n_quads=surface_animation(output_surface)%n_quads
  n_frames=surface_animation(output_surface)%n_frames
  
  do frame=1,n_frames
  
    write(*,*)'Writing frame data',frame
    
! add frame number to base filename

    CALL add_integer_to_filename(temp_filename3,frame,temp_filename4)

    frame_filename=trim(temp_filename4)//'.vtk'
    
    write(*,*)'Opening file:',trim(frame_filename)
    
    open(UNIT=animation_output_unit,FILE=trim(frame_filename),ERR=9010)
    
! write header information    
! write vtk header      
    write(animation_output_unit,'(A)')'# vtk DataFile Version 2.0'
    write(animation_output_unit,'(A)')'trim(frame_filename)'
    write(animation_output_unit,'(A)')'ASCII'
    write(animation_output_unit,'(A)')'DATASET POLYDATA'
    write(animation_output_unit,'(A,I10,A)')'POINTS',n_points,' float'
    
! write point data 
    point=0
    do quad=1,n_quads
    
! write 4 points for each quad    
      point=point+1
      write(animation_output_unit,8000)surface_animation(output_surface)%points(point,1),	&
                                       surface_animation(output_surface)%points(point,2),	&
				       surface_animation(output_surface)%points(point,3)
      point=point+1
      write(animation_output_unit,8000)surface_animation(output_surface)%points(point,1),	&
                                       surface_animation(output_surface)%points(point,2),	&
				       surface_animation(output_surface)%points(point,3)
      point=point+1
      write(animation_output_unit,8000)surface_animation(output_surface)%points(point,1),	&
                                       surface_animation(output_surface)%points(point,2),	&
				       surface_animation(output_surface)%points(point,3)
      point=point+1
      write(animation_output_unit,8000)surface_animation(output_surface)%points(point,1),	&
                                       surface_animation(output_surface)%points(point,2),	&
				       surface_animation(output_surface)%points(point,3)
      
8000  format(3E14.5)
    	  
    end do ! next quad
  
! write quad data
    write(animation_output_unit,'(A,2I10)')'POLYGONS',n_quads,n_quads*5

    point=0
    do quad=1,n_quads
    
      write(animation_output_unit,8010)4,point,point+1,point+2,point+3
      point=point+4

8010  format(I3,4I8)
      
    end do ! next quad

  
! STAGE 6. write data associated with points

! write point based data
    write(animation_output_unit,'(A,I10)')'POINT_DATA ',n_points
    write(animation_output_unit,'(A)')'SCALARS Field_on_cells float 1'
    write(animation_output_unit,'(A)')'LOOKUP_TABLE field_on_cells_table'

    point=0
    do quad=1,n_quads
    
      re=dble(surface_animation(output_surface)%complex_data(quad))
      im=dimag(surface_animation(output_surface)%complex_data(quad))
      
      phase=2d0*pi*dble(frame)/dble(n_frames)
      
      point_data=re*cos(phase)-im*sin(phase)
      
! write 4 data points for each quad    
      write(animation_output_unit,8020)point_data
      write(animation_output_unit,8020)point_data
      write(animation_output_unit,8020)point_data
      write(animation_output_unit,8020)point_data
      
8020  format(E14.5)
    	  
    end do ! next quad

    close(UNIT=animation_output_unit)

  end do ! next frame
  
! STAGE 7. Close fieldsolve file and deallocate memory

  
  CLOSE(unit=local_file_unit)
  
  GOTO 100 ! next animation output
  
  
  RETURN
  
9000 continue
  write(*,*)'Error opening output file'
  write(*,*)'Filename='
  write(*,*)base_filename
  stop
9010 continue
  write(*,*)'Error opening frame file'
  write(*,*)'Filename='
  write(*,*)frame_filename
  stop

  RETURN
  

  
END SUBROUTINE create_frequency_domain_animation
