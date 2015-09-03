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
! SUBROUTINE create_poynting_vector_plot
!
! NAME
!    create_poynting_vector_plot
!
! DESCRIPTION
!    Read frequency domain vector E and H fields and calculate the Poynting vector
!    Then plot the resulting vector field
!    Now extended to plot real and imaginary E and H fioeld vectors
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 12/6/2015 CJS
!     plot real and imaginary E and H fioeld vectors 3/9/2015 CJS
!
SUBROUTINE create_poynting_vector_plot

USE post_process
USE file_information
USE constants

IMPLICIT NONE

! local variables

  character(len=256)	:: filename
    
  character*256 base_filename
  character*256 output_filename
  character*256 temp_filename
  character*256 temp_filename2
  character*256 temp_filename3
  character*256 temp_filename4
  character*256 frame_filename
  
  integer n_volumes,surface_loop,surface
  integer output_volume
  integer template_volume
  integer output_volume_xyz(6)
  
  TYPE(surface_animation_data),allocatable	:: volume_animation(:)
  
  integer			:: n_points,n_quads,n_cells,n_frames
  integer			:: n_points_arrow,n_quads_arrow
  integer			:: n_points_cell,n_quads_cell
  
  real*8 			:: max_E
  real*8 			:: min_E
  real*8 			:: max_H
  real*8 			:: min_H
  real*8 			:: max_data
  real*8 			:: data_range
  
  complex*16			:: Efield(3)
  complex*16			:: Hfield(3)
  complex*16			:: PV(3)
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
  integer 	:: point,quad,cell
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
  
  logical	:: found_volume
  integer 	:: number_of_cone_surfaces
  logical	:: scale_vector
  
  integer	:: optype
  integer,parameter :: avg_power=1
  integer,parameter :: reactive_power=2
  integer,parameter :: E_re=3
  integer,parameter :: E_im=4
  integer,parameter :: H_re=5
  integer,parameter :: H_im=6
  
  integer	:: plot_type
  integer,parameter :: vector_plot=1
  integer,parameter :: magnitude_plot=2
  
  character*3 :: ch3
  logical	:: logscale_flag
  
! START

  number_of_cone_surfaces=8

! STAGE 1. Open data file

5 write(*,*)
  write(*,*)'Volume field output files:'
  
  command='ls -ltr *.frequency_output_volume.fout'
  CALL system(command)

  write(*,*)'Enter the volume field filename'
  read(*,*)filename
  inquire(file=trim(filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')trim(filename)
  
  OPEN(unit=local_file_unit,file=filename)
   
  base_filename=trim(filename)
    
! STAGE 2. Read frequency_output_volume file

!# NUMBER OF FREQUENCY OUTPUT VOLUMES:
!2
!# START SURFACE FIELD FILE TEMPLATE, OUTPUT VOLUME NUMBER:    1
!# Number of points:'
!20
!# Number of faces:'
!80
!# Number of frames:'
!24

  write(*,*)'Reading frequency_output_volume file:'
  
  read(local_file_unit,'(A80)'),ipline
  read(local_file_unit,*)n_volumes
  
  ALLOCATE ( volume_animation(1:n_volumes) )
#
! read header information  
  do surface=1,n_volumes
  
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,*)n_points
    volume_animation(surface)%n_points=n_points
     
    write(*,*)'Number of points=',volume_animation(surface)%n_points
  
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,*)n_quads
    volume_animation(surface)%n_quads=n_quads
     
    write(*,*)'Number of quads=',volume_animation(surface)%n_quads
  
    read(local_file_unit,'(A80)'),ipline
    read(local_file_unit,*)n_frames
    volume_animation(surface)%n_frames=n_frames
     
    write(*,*)'Number of frames=',volume_animation(surface)%n_frames
    
    n_cells=n_points/8
    if (n_cells.ne.n_quads/6) then
      write(*,*)'Error in create_poynting_vector_plot'
      write(*,*)'n_cells.ne.n_quads/6'
      write(*,*)'n_points=',n_points,' n_points/8=',n_points/8
      write(*,*)'n_quads =',n_quads, ' n_quads/6 =',n_quads/6
      STOP
    end if

    write(*,*)'Allocating data'
       
    ALLOCATE ( volume_animation(surface)%points(1:n_points,1:3) ) 
    ALLOCATE ( volume_animation(surface)%quads(1:n_quads,1:4) ) 
    ALLOCATE ( volume_animation(surface)%complex_data(1:n_cells) ) 
    ALLOCATE ( volume_animation(surface)%magnitude_data(1:n_cells) ) 
  
    volume_animation(surface)%points(1:n_points,1:3)=0D0
    volume_animation(surface)%quads(1:n_quads,1:4)=0
    volume_animation(surface)%complex_data(1:n_cells)=0D0
    volume_animation(surface)%magnitude_data(1:n_cells)=0D0
  
    write(*,*)'Reading points'

    do point=1,n_points
      read(local_file_unit,*)x,y,z
      volume_animation(surface)%points(point,1)=x
      volume_animation(surface)%points(point,2)=y
      volume_animation(surface)%points(point,3)=z
    end do
    
! read quads    
    write(*,*)'Reading quads'

    do quad=1,n_quads
      read(local_file_unit,*)ip_quad,point1,point2,point3,point4
      volume_animation(surface)%quads(quad,1)=point1
      volume_animation(surface)%quads(quad,2)=point2
      volume_animation(surface)%quads(quad,3)=point3
      volume_animation(surface)%quads(quad,4)=point4
    end do
    
    volume_animation(surface)%max_data=-1e30
    volume_animation(surface)%min_data=1e30
    
! read data    
    write(*,*)'Reading data'
    do quad=1,n_cells
    
      read(local_file_unit,*)cx,cy,cz,re,im,mag
      
      volume_animation(surface)%complex_data(quad)=cmplx(re,im)
      
      volume_animation(surface)%max_data=max(volume_animation(surface)%max_data,mag)
      volume_animation(surface)%min_data=min(volume_animation(surface)%min_data,mag)
      
    end do

    write(*,*)'Maximum data value:',volume_animation(surface)%max_data
    write(*,*)'Minimum data value:',volume_animation(surface)%min_data
    
    write(*,*)'Finished reading geometry data'
      
  end do ! read next volume
 
  write(*,*)' '
  write(*,*)'Number of volumes read:',n_volumes
  write(*,*)' '
  
! read the number of the volume to write to file  

100 CONTINUE  ! read volume numbers for the animation
  
  do i=1,6
  
    if (i.eq.1) then
      write(*,*)'Enter the number of the volume for the Efield in the x direction or 0 if this component does not exixst'
    else if(i.eq.2) then
      write(*,*)'Enter the number of the volume for the Efield in the y direction or 0 if this component does not exixst'
    else if(i.eq.3) then
      write(*,*)'Enter the number of the volume for the Efield in the z direction or 0 if this component does not exixst'
    else if(i.eq.4) then
      write(*,*)'Enter the number of the volume for the Hfield in the x direction or 0 if this component does not exixst'
    else if(i.eq.5) then
      write(*,*)'Enter the number of the volume for the Hfield in the y direction or 0 if this component does not exixst'
    else if(i.eq.6) then
      write(*,*)'Enter the number of the volume for the Hfield in the z direction or 0 if this component does not exixst'
    end if 
    
    read(*,*)output_volume_xyz(i)
    if ( (output_volume_xyz(i).lt.0).OR.(output_volume_xyz(i).gt.n_volumes) )then
      write(*,*)'Output volume must be between 0 and ',n_volumes
      GOTO 100
    end if
    write(record_user_inputs_unit,*)output_volume_xyz(i),' Output volume component number',i
    
  end do

! Check that we have the same amount of data for each field component
! Also work out scaling values for data

  max_E=-1e30
  min_E=1d30
  max_H=-1e30
  min_H=1d30

  found_volume=.FALSE.
  do i=1,6
    surface=output_volume_xyz(i)
    
    write(*,*)'checking i=',i,' surface=',surface
    
    if ( (surface.ne.0).AND.(found_volume) ) then
! An output volume has already been found so just check that the data is consistent 
    
      if (  (n_quads.NE.volume_animation(surface)%n_quads).OR.	&
            (n_points.NE.volume_animation(surface)%n_points).OR.	&
            (n_frames.NE.volume_animation(surface)%n_frames ) ) then
	write(*,*)'Error in create_poynting_vector_plot'
	write(*,*)'Component volumes have different amounts of data'
	STOP
      end if

    else if (surface.ne.0) then
! this is the first component of the volume data found    

      if (.NOT.found_volume) then 
! this is the first component of the volume data found so record some data from this surface for checking
        template_volume=surface
        n_quads=volume_animation(surface)%n_quads
        n_points=volume_animation(surface)%n_points
        n_frames=volume_animation(surface)%n_frames
        found_volume=.TRUE.
      end if
    end if

! get min and max data from ALL volumes   
    if (surface.ne.0) then    
      
      if (i.LE.3) then
        max_E=max(volume_animation(surface)%max_data,max_E)
        min_E=min(volume_animation(surface)%min_data,min_E)
      else
        max_H=max(volume_animation(surface)%max_data,max_H)
        min_H=min(volume_animation(surface)%min_data,min_H)  
      end if
            
    end if
    
  end do

  write(*,*)'Maximum E field value:',max_E
  write(*,*)'Minimum E field value:',min_E
  write(*,*)'Maximum H field value:',max_H
  write(*,*)'Minimum H field value:',min_H
  
    
  write(*,*)
  write(*,*)'Do you want to plot Average power flow (1), Reactive power (2), Real E (3), Imag E (4), Real H (5) or Imag H (6) ?'

200 CONTINUE  
  read(*,*)optype
  if ( (optype.LT.avg_power).AND.(optype.GT.H_im) )then
    write(*,*)'Expecting either 1-6 as a response...'
    GOTO 200
  end if
  write(record_user_inputs_unit,*)optype,&
   ' Average power flow (1), Reactive power (2), Real E (3), Imag E (4), Real H (5) or Imag H (6)'
    
  write(*,*)
  write(*,*)'Do you want a vector plot (1) or magnitude plot (2) ?'

300 CONTINUE  
  read(*,*)plot_type
  if ( (plot_type.NE.vector_plot).AND.(plot_type.NE.magnitude_plot) )then
    write(*,*)'Expecting either 1 or 2 as a response...'
    GOTO 300
  end if
  write(record_user_inputs_unit,*)plot_type,' vector plot (1) or magnitude plot (2) '
    
  write(*,*)
  write(*,*)'Linear or log scale output? (lin/log)'

  read(*,'(A3)')ch3
  if (ch3(1:3).EQ.'log') then
    logscale_flag=.TRUE.
  else if (ch3(1:3).EQ.'lin') then
    logscale_flag=.FALSE.  
  else
    write(*,*)'Shoule be lin or log'
    STOP
  end if
  write(record_user_inputs_unit,'(A3)')ch3

! work out the maximum length of the vector to plot  

  data_range=0d0
  
  do cell=1,n_cells

! get the complex E and H fields for this cell  
      
    Efield(:)=(0d0,0d0)
    Hfield(:)=(0d0,0d0)
    PV(:)=(0d0,0d0)
      
    do i=1,6
      
      surface=output_volume_xyz(i)
      if (surface.ne.0) then
    
	if (i.LE.3) then
	  Efield(i)=volume_animation(surface)%complex_data(cell)
	else
          Hfield(i-3)=volume_animation(surface)%complex_data(cell)
	end if
	      
      end if

    end do  ! next field component

! get the plot vector
    if ( (optype.EQ.avg_power).OR.(optype.EQ.reactive_power) )then
! Calculate the complex Poynting vector P=ExH*      
      PV(1)=Efield(2)*conjg(Hfield(3))-Efield(3)*conjg(Hfield(2))
      PV(2)=Efield(3)*conjg(Hfield(1))-Efield(1)*conjg(Hfield(3))
      PV(3)=Efield(1)*conjg(Hfield(2))-Efield(2)*conjg(Hfield(1))
    else if ( optype.EQ.E_re )then
      PV(1)=cmplx(Dble(Efield(1)))
      PV(2)=cmplx(Dble(Efield(2)))
      PV(3)=cmplx(Dble(Efield(3)))
    else if ( optype.EQ.H_re )then
      PV(1)=cmplx(Dble(Hfield(1)))
      PV(2)=cmplx(Dble(Hfield(2)))
      PV(3)=cmplx(Dble(Hfield(3))) 
    else if ( optype.EQ.E_im )then
      PV(1)=cmplx(Dble(-j*Efield(1)))
      PV(2)=cmplx(Dble(-j*Efield(2)))
      PV(3)=cmplx(Dble(-j*Efield(3)))
    else if ( optype.EQ.H_im )then
      PV(1)=cmplx(Dble(-j*Hfield(1)))
      PV(2)=cmplx(Dble(-j*Hfield(2)))
      PV(3)=cmplx(Dble(-j*Hfield(3))) 
    end if
     
    if (optype.eq.avg_power) then
      data_range=max(data_range,real(PV(1))/2,real(PV(2))/2,real(PV(3))/2)
    else
      data_range=max(data_range,imag(PV(1))/2,imag(PV(2))/2,imag(PV(3))/2)
    end if
    
  end do
  
! assume that max H and ,max E are at the same point to estimate the data range
  if (data_range.eq.0d0) data_range=1d0
  
  write(*,*)'Mamimum vector component length to plot=',data_range
  
  write(*,*)'Do you want to scale the size of the vector? (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A)')ch  
  CALL convert_to_lower_case(ch,1)
  scale_vector=.FALSE.
  if (ch.eq.'y') then
    scale_vector=.TRUE.
  end if

  write(*,*)'Enter the filename for the Poynting vector output plot'
  read(*,'(A256)')output_filename
  write(record_user_inputs_unit,'(A)')trim(output_filename)
  
! add surface number to base filename

  temp_filename=trim(output_filename)

  n_points=volume_animation(template_volume)%n_points
  n_quads=volume_animation(template_volume)%n_quads
  n_frames=1   ! no e(jwt) dependence in the Poynting vector - only write a single frame
  
  n_cells=n_points/8   ! bug fix 22/11/2013 CJS 
  
  n_quads_arrow=n_cells*number_of_cone_surfaces
  n_points_arrow=n_quads_arrow*4
  
  n_quads_cell=n_cells*6
  n_points_cell=n_quads_cell*4
    
  n_frames=1  
    
  do frame=1,n_frames

    frame_filename=trim(temp_filename)//'.vtk'
    
    write(*,*)'Opening file:',trim(frame_filename)
    
    open(UNIT=animation_output_unit,FILE=trim(frame_filename),ERR=9010)
    
! write header information    
! write vtk header      
    write(animation_output_unit,'(A)')'# vtk DataFile Version 2.0'
    write(animation_output_unit,'(A)')'trim(frame_filename)'
    write(animation_output_unit,'(A)')'ASCII'
    write(animation_output_unit,'(A)')'DATASET POLYDATA'
    
    if (plot_type.eq.vector_plot) then
      write(animation_output_unit,'(A,I10,A)')'POINTS',n_points_arrow,' float'
    else if (plot_type.eq.magnitude_plot) then
      write(animation_output_unit,'(A,I10,A)')'POINTS',n_points_cell,' float'
    end if
    
! write cell data 
    point=0
    do cell=1,n_cells

! get the complex E and H fields for this cell  
      
      Efield(:)=(0d0,0d0)
      Hfield(:)=(0d0,0d0)
      PV(:)=(0d0,0d0)
      
      do i=1,6
      
        surface=output_volume_xyz(i)
        if (surface.ne.0) then
    
	  if (i.LE.3) then
	    Efield(i)=volume_animation(surface)%complex_data(cell)
	  else
	    Hfield(i-3)=volume_animation(surface)%complex_data(cell)
	  end if
	      
        end if
	
      end do  ! next field component
      
      do i=1,3
        centre(i)=0.125d0*( volume_animation(template_volume)%points(point+1,i)+	&
                            volume_animation(template_volume)%points(point+2,i)+	&
                            volume_animation(template_volume)%points(point+3,i)+	&
                            volume_animation(template_volume)%points(point+4,i)+	&
                            volume_animation(template_volume)%points(point+5,i)+	&
                            volume_animation(template_volume)%points(point+6,i)+	&
                            volume_animation(template_volume)%points(point+7,i)+	&
                            volume_animation(template_volume)%points(point+8,i) ) 
      end do
                    
      cell_size =sqrt( ( volume_animation(template_volume)%points(point+1,1)		&
                        -volume_animation(template_volume)%points(point+2,1) )**2+	&
                       ( volume_animation(template_volume)%points(point+1,2)		&
                        -volume_animation(template_volume)%points(point+2,2) )**2+	&
                       ( volume_animation(template_volume)%points(point+1,3)		&
                        -volume_animation(template_volume)%points(point+2,3) )**2 )

! get the plot vector
    if ( (optype.EQ.avg_power).OR.(optype.EQ.reactive_power) )then
! Calculate the complex Poynting vector P=ExH*      
      PV(1)=Efield(2)*conjg(Hfield(3))-Efield(3)*conjg(Hfield(2))
      PV(2)=Efield(3)*conjg(Hfield(1))-Efield(1)*conjg(Hfield(3))
      PV(3)=Efield(1)*conjg(Hfield(2))-Efield(2)*conjg(Hfield(1))
    else if ( optype.EQ.E_re )then
      PV(1)=cmplx(Dble(Efield(1)))
      PV(2)=cmplx(Dble(Efield(2)))
      PV(3)=cmplx(Dble(Efield(3)))
    else if ( optype.EQ.H_re )then
      PV(1)=cmplx(Dble(Hfield(1)))
      PV(2)=cmplx(Dble(Hfield(2)))
      PV(3)=cmplx(Dble(Hfield(3))) 
    else if ( optype.EQ.E_im )then
      PV(1)=cmplx(Dble(-j*Efield(1)))
      PV(2)=cmplx(Dble(-j*Efield(2)))
      PV(3)=cmplx(Dble(-j*Efield(3)))
    else if ( optype.EQ.H_im )then
      PV(1)=cmplx(Dble(-j*Hfield(1)))
      PV(2)=cmplx(Dble(-j*Hfield(2)))
      PV(3)=cmplx(Dble(-j*Hfield(3))) 
    end if

! calculate the vector to plot
      do i=1,3
	
! must decide what to plot here...

        re=real(PV(i))
	im=imag(PV(i))
		
        if (optype.eq.avg_power) then
          vector(i)=0.5d0*re  ! Real power
	else if (optype.eq.reactive_power) then
          vector(i)=0.5d0*im  ! Reactive power
	else
          vector(i)=re
	end if
	
      end do

      length=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
            
      if (length.eq.0d0) then
! we must give the plotting routine a non-zero length vector so just send something very small...
 
        vector(1)=data_range/1000d0
	length=data_range/1000d0
 
      end if
      
      if (logscale_flag) then
        volume_animation(template_volume)%magnitude_data(cell)=10d0*log10(length)
      else
        volume_animation(template_volume)%magnitude_data(cell)=length
      end if
      
      if (plot_type.eq.vector_plot) then
      
        if (scale_vector) then
          vector(:)=vector(:)*cell_size/data_range
          diameter=(cell_size*length/data_range)/4d0
        else
          vector(:)=vector(:)*cell_size/length
          diameter=cell_size/4d0
        end if
       
        arrow_points(1,1:3)=centre(1:3)-vector(1:3)/2d0
        arrow_points(2,1:3)=centre(1:3)+vector(1:3)/2d0
      
        CALL write_cone_points_vtk(arrow_points(1,1),arrow_points(1,2),arrow_points(1,3),	&
                                   arrow_points(2,1),arrow_points(2,2),arrow_points(2,3),	&
	  		         diameter,number_of_cone_surfaces,animation_output_unit)
	
      else if (plot_type.eq.magnitude_plot) then
      
       
        CALL write_cube_points_vtk(volume_animation(template_volume)%points(point+1,1),&
	                           volume_animation(template_volume)%points(point+1,2),&
	                           volume_animation(template_volume)%points(point+1,3),&
	                           volume_animation(template_volume)%points(point+2,1),&
	                           volume_animation(template_volume)%points(point+2,2),&
	                           volume_animation(template_volume)%points(point+2,3),&
	                           volume_animation(template_volume)%points(point+3,1),&
	                           volume_animation(template_volume)%points(point+3,2),&
	                           volume_animation(template_volume)%points(point+3,3),&
	                           volume_animation(template_volume)%points(point+4,1),&
	                           volume_animation(template_volume)%points(point+4,2),&
	                           volume_animation(template_volume)%points(point+4,3),&
	                           volume_animation(template_volume)%points(point+5,1),&
	                           volume_animation(template_volume)%points(point+5,2),&
	                           volume_animation(template_volume)%points(point+5,3),&
	                           volume_animation(template_volume)%points(point+6,1),&
	                           volume_animation(template_volume)%points(point+6,2),&
	                           volume_animation(template_volume)%points(point+6,3),&
	                           volume_animation(template_volume)%points(point+7,1),&
	                           volume_animation(template_volume)%points(point+7,2),&
	                           volume_animation(template_volume)%points(point+7,3),&
	                           volume_animation(template_volume)%points(point+8,1),&
	                           volume_animation(template_volume)%points(point+8,2),&
	                           volume_animation(template_volume)%points(point+8,3),&
				   animation_output_unit)
      
      end if
      
      point=point+8
      
8000  format(3E14.5)
    	  
    end do ! next cell
  
! write quad data

    if (plot_type.eq.vector_plot) then

      write(animation_output_unit,'(A,2I12)')'POLYGONS',n_quads_arrow,n_quads_arrow*5

      point=0
      do quad=1,n_quads_arrow
    
        write(animation_output_unit,8010)4,point,point+1,point+2,point+3
        point=point+4

8010    format(I3,4I12)
      
      end do ! next cell

! write point based data
      write(animation_output_unit,'(A,I12)')'POINT_DATA ',n_points_arrow
      write(animation_output_unit,'(A)')'SCALARS Field_on_cells float 1'
      write(animation_output_unit,'(A)')'LOOKUP_TABLE field_on_cells_table'

      point=0
      do cell=1,n_cells
      
        point_data=volume_animation(template_volume)%magnitude_data(cell)
      
        do i=1,number_of_cone_surfaces*4
          write(animation_output_unit,8020)point_data
        end do
      
8020    format(E14.5)
    	  
      end do ! next quad
      
    else if (plot_type.eq.magnitude_plot) then

      write(animation_output_unit,'(A,2I12)')'POLYGONS',n_quads_cell,n_quads_cell*5

      point=0
      do quad=1,n_quads_cell
    
        write(animation_output_unit,8010)4,point,point+1,point+2,point+3
        point=point+4
      
      end do ! next cell

! write point based data
      write(animation_output_unit,'(A,I12)')'POINT_DATA ',n_points_cell
      write(animation_output_unit,'(A)')'SCALARS Field_on_cells float 1'
      write(animation_output_unit,'(A)')'LOOKUP_TABLE field_on_cells_table'

      point=0
      do cell=1,n_cells
      
        point_data=volume_animation(template_volume)%magnitude_data(cell)
      
        do i=1,24    ! 4 points on each of 6 faces
          write(animation_output_unit,8020)point_data
        end do
    	  
      end do ! next quad

    end if

    close(UNIT=animation_output_unit)

  end do ! next frame
  
! STAGE 7. Close fieldsolve file and deallocate memory

  CLOSE(unit=local_file_unit)
    
! Deallocate memory for animation  
    
  do surface=1,n_volumes
    if ( allocated( volume_animation(surface)%points ) ) DEALLOCATE ( volume_animation(surface)%points ) 
    if ( allocated( volume_animation(surface)%quads ) ) DEALLOCATE ( volume_animation(surface)%quads ) 
    if ( allocated( volume_animation(surface)%complex_data ) ) DEALLOCATE ( volume_animation(surface)%complex_data ) 
    if ( allocated( volume_animation(surface)%magnitude_data ) ) DEALLOCATE ( volume_animation(surface)%magnitude_data ) 
  end do
    
  if ( allocated( volume_animation ) ) DEALLOCATE ( volume_animation )  
  
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
  

  
END SUBROUTINE create_poynting_vector_plot
