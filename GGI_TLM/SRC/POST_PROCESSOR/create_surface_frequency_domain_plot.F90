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
! SUBROUTINE create_surface_frequency_domain_plot
!
! NAME
!    create_surface_frequency_domain_plot
!
! DESCRIPTION
!     
!     
! COMMENTS
!  This uses the surface animation structures with a single frame only to show magnitude data     
!
! HISTORY
!
!     started 21/08/2017 CJS 
!
!
SUBROUTINE create_surface_frequency_domain_plot

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
  character*256 output_filename
  character*256 plot_filename
  
  integer n_surfaces,surface_loop,surface
  integer output_surface
  integer template_surface
  integer output_surface_xyz
  
  TYPE(surface_animation_data),allocatable	:: surface_animation(:)
  
  integer			:: n_points,n_quads,n_frames
  
  real*8 			:: max_data
  real*8 			:: min_data
  real*8 			:: data_range
  real*8 			:: scale_factor
  
  real*8			:: field(3)
  real*8			:: centre(3)
  real*8			:: cell_size
  real*8			:: vector(3)
  real*8			:: length,diameter
  real*8			:: plot_length
  real*8			:: arrow_points(2,3)
  
  integer max_point_number

  integer max_quad_number
  
  character*256 :: ipline

  character(len=256)	:: command
  
  logical	:: file_exists
  
  integer 	:: ip_point,ip_quad
  integer 	:: point,quad
  integer 	:: point1,point2,point3,point4
  real*8 	:: x,y,z
  
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
  
  logical	:: found_surface
  integer 	:: number_of_cone_surfaces
  logical	:: scale_vector
! START

  number_of_cone_surfaces=8

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
    ALLOCATE ( surface_animation(surface)%magnitude_data(1:n_quads) ) 
  
    surface_animation(surface)%points(1:n_points,1:3)=0D0
    surface_animation(surface)%quads(1:n_quads,1:4)=0
    surface_animation(surface)%complex_data(1:n_quads)=0D0
    surface_animation(surface)%magnitude_data(1:n_quads)=0D0
  
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

    write(*,*)'Maximum data value:',surface_animation(surface)%max_data
    write(*,*)'Minimum data value:',surface_animation(surface)%min_data
    
    write(*,*)'Finished reading fieldsolve geometry data'
      
  end do ! read next surface
 
  write(*,*)' '
  write(*,*)'Number of surfaces read:',n_surfaces
  write(*,*)' '
  
! read the number of the surface to write to file  

100 CONTINUE  ! read surface numbers for the animation

  write(*,*)'Enter the number of the surface for the field plot data'
  read(*,*)output_surface_xyz
 
  write(record_user_inputs_unit,*)output_surface_xyz,' Output surface  number'
    
! Work out range of data

  max_data=-1e30
  min_data=1d30

  surface=output_surface_xyz
  n_quads=surface_animation(surface)%n_quads
  n_points=surface_animation(surface)%n_points
  n_frames=surface_animation(surface)%n_frames
  max_data=max(surface_animation(surface)%max_data,max_data)
  min_data=min(surface_animation(surface)%min_data,min_data)
 
  write(*,*)'Maximum data value:',max_data
  write(*,*)'Minimum data value:',min_data
  
  data_range=(max_data-min_data)
  if (data_range.eq.0d0) data_range=1d0

  write(*,*)'Enter the scale factor for the field plot data'
  read(*,*)scale_factor
 
  write(record_user_inputs_unit,*)scale_factor,' Scale factor'

  write(*,*)'Enter the filename for the plot'
  read(*,'(A256)')output_filename
  write(record_user_inputs_unit,'(A)')trim(output_filename) 
  
! add surface number to base filename

  plot_filename=trim(output_filename)//'.vtk'

  write(*,*)'Writing plot data'
        
  open(UNIT=animation_output_unit,FILE=trim(plot_filename),ERR=9010)
    
! write header information    
! write vtk header      
  write(animation_output_unit,'(A)')'# vtk DataFile Version 2.0'
  write(animation_output_unit,'(A)')trim(plot_filename)
  write(animation_output_unit,'(A)')'ASCII'
  write(animation_output_unit,'(A)')'DATASET POLYDATA'
  write(animation_output_unit,'(A,I10,A)')'POINTS',n_points,' float'
    
! write point data 
  point=0
  do quad=1,n_quads
    
! write 4 points for each quad    
    point=point+1
    write(animation_output_unit,8000)surface_animation(surface)%points(point,1),       &
                                     surface_animation(surface)%points(point,2),       &
                                     surface_animation(surface)%points(point,3)
    point=point+1
    write(animation_output_unit,8000)surface_animation(surface)%points(point,1),       &
                                     surface_animation(surface)%points(point,2),       &
                                     surface_animation(surface)%points(point,3)
    point=point+1
    write(animation_output_unit,8000)surface_animation(surface)%points(point,1),       &
                                     surface_animation(surface)%points(point,2),       &
                                     surface_animation(surface)%points(point,3)
    point=point+1
    write(animation_output_unit,8000)surface_animation(surface)%points(point,1),       &
                                     surface_animation(surface)%points(point,2),       &
                                     surface_animation(surface)%points(point,3)
    
8000  format(3E14.5)
    	  
  end do ! next quad
  
! write quad data
   write(animation_output_unit,'(A,2I10)')'POLYGONS',n_quads,n_quads*5

  point=0
  do quad=1,n_quads
    
    write(animation_output_unit,8010)4,point,point+1,point+2,point+3
    point=point+4

8010  format(I3,4I12)
      
  end do ! next quad
  
! STAGE 6. write data associated with points

! write point based data
  write(animation_output_unit,'(A,I10)')'POINT_DATA ',n_points
  write(animation_output_unit,'(A)')'SCALARS Field_on_cells float 1'
  write(animation_output_unit,'(A)')'LOOKUP_TABLE field_on_cells_table'

  point=0
  do quad=1,n_quads

    point_data=abs(surface_animation(surface)%complex_data(quad))*scale_factor
    
! write 4 data points for each quad    
    write(animation_output_unit,8020)point_data
    write(animation_output_unit,8020)point_data
    write(animation_output_unit,8020)point_data
    write(animation_output_unit,8020)point_data
      
8020  format(E14.5)
    	  
  end do ! next quad

  close(UNIT=animation_output_unit)
  
! STAGE 7. Close fieldsolve file and deallocate memory

  CLOSE(unit=local_file_unit)
    
! Deallocate memory for animation  
    
  do surface=1,n_surfaces
    if ( allocated( surface_animation(surface)%points ) ) DEALLOCATE ( surface_animation(surface)%points ) 
    if ( allocated( surface_animation(surface)%quads ) ) DEALLOCATE ( surface_animation(surface)%quads ) 
    if ( allocated( surface_animation(surface)%complex_data ) ) DEALLOCATE ( surface_animation(surface)%complex_data ) 
    if ( allocated( surface_animation(surface)%magnitude_data ) ) DEALLOCATE ( surface_animation(surface)%magnitude_data ) 
  end do
    
  if ( allocated( surface_animation ) ) DEALLOCATE ( surface_animation )  
  
  RETURN
  
9000 continue
  write(*,*)'Error opening output file'
  write(*,*)'Filename='
  write(*,*)base_filename
  stop
9010 continue
  write(*,*)'Error opening plot file'
  write(*,*)'Filename='
  write(*,*)plot_filename
  stop

  RETURN
  

  
END SUBROUTINE create_surface_frequency_domain_plot
