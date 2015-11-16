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
! SUBROUTINE Convert_Frequency_Domain_Volume_Field_Format
!
! NAME
!    Convert_Frequency_Domain_Volume_Field_Format
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 23/10/2014 CJS
!     Fix problem when number of cells is not the same for all volumes... CJS 16/1/2015
!     5/11/2015  CJS allow the process to work on compressed files
!
SUBROUTINE Convert_Frequency_Domain_Volume_Field_Format

USE post_process
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer :: run_number  
  
integer :: n_volumes
integer :: n_points
integer :: n_cells
integer :: n_quads
integer :: n_frames

integer :: volume,point,quad
integer :: ip_quad,point1,point2,point3,point4
integer :: cx,cy,cz
real*8	:: x,y,z
real*8	:: re,im

integer :: i

character(len=256) :: filename
  
character(len=256)	:: command
character		:: ch
  
logical	:: file_exists

TYPE(surface_animation_data),allocatable	:: volume_animation(:)

! compression stuff
  logical	:: compression_flag
  integer	:: len_filename
  character(len=256)	:: filename2
  character*3   :: extn

! START
! Open data file

5 write(*,*)
  write(*,*)'Volume field output files:'
  
  command='ls -ltr *.frequency_output_volume.fout*'
  CALL system(command)

  write(*,*)'Enter the volume field filename'
  read(*,'(A)')filename
  inquire(file=filename,exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')trim(filename)

! check for .gz extension which indicates a compressed file
  compression_flag=.FALSE.
  len_filename=LEN(trim(filename))
  extn=filename(len_filename-2:len_filename)
  if (extn.EQ.'.gz') then
    compression_flag=.TRUE.
  end if
  
  if (compression_flag) then
! We must give the filename without the .gz extension here
    filename2=filename(1:len_filename-3)
    CALL open_output_file_read(local_file_unit,filename2,compression_flag)
  else
    OPEN(unit=local_file_unit,file=filename)
  end if
      
! STAGE 2. Read frequency_output_volume file

!# NUMBER OF FREQUENCY OUTPUT VOLUMES:
!2
!# START VOLUME FIELD FILE TEMPLATE, OUTPUT VOLUME NUMBER:    1
!# Number of points:'
!20
!# Number of faces:'
!80
!# Number of frames:'
!24

  write(*,*)'Reading frequency_output_volume file:'
  
  read(local_file_unit,'(A80)')
  read(local_file_unit,*)n_volumes
  
  ALLOCATE ( volume_animation(1:n_volumes) )
#
! read header information  
  do volume=1,n_volumes
  
    read(local_file_unit,'(A80)')
    read(local_file_unit,'(A80)')
    read(local_file_unit,*)n_points
    volume_animation(volume)%n_points=n_points
     
    write(*,*)'Number of points=',volume_animation(volume)%n_points
  
    read(local_file_unit,'(A80)')
    read(local_file_unit,*)n_quads
    volume_animation(volume)%n_quads=n_quads
     
    write(*,*)'Number of quads=',volume_animation(volume)%n_quads
  
    read(local_file_unit,'(A80)')
    read(local_file_unit,*)n_frames
    volume_animation(volume)%n_frames=n_frames
     
    write(*,*)'Number of frames=',volume_animation(volume)%n_frames
    
    n_cells=n_points/8
    if (n_cells.ne.n_quads/6) then
      write(*,*)'Error in create_vector_volume_frequency_domain_animation'
      write(*,*)'n_cells.ne.n_quads/6'
      write(*,*)'n_points=',n_points,' n_points/8=',n_points/8
      write(*,*)'n_quads =',n_quads, ' n_quads/6 =',n_quads/6
      STOP
    end if

    write(*,*)'Allocating data'
       
    ALLOCATE ( volume_animation(volume)%points(1:n_points,1:3) ) 
    ALLOCATE ( volume_animation(volume)%quads(1:n_quads,1:4) ) 
    ALLOCATE ( volume_animation(volume)%complex_data(1:n_cells) ) 
    ALLOCATE ( volume_animation(volume)%magnitude_data(1:n_cells) ) 
  
    volume_animation(volume)%points(1:n_points,1:3)=0D0
    volume_animation(volume)%quads(1:n_quads,1:4)=0
    volume_animation(volume)%complex_data(1:n_cells)=0D0
    volume_animation(volume)%magnitude_data(1:n_cells)=0D0
  
    write(*,*)'Reading points'

    do point=1,n_points
      read(local_file_unit,*)x,y,z
      volume_animation(volume)%points(point,1)=x
      volume_animation(volume)%points(point,2)=y
      volume_animation(volume)%points(point,3)=z
    end do
    
! read quads    
    write(*,*)'Reading quads'

    do quad=1,n_quads
      read(local_file_unit,*)ip_quad,point1,point2,point3,point4
      volume_animation(volume)%quads(quad,1)=point1
      volume_animation(volume)%quads(quad,2)=point2
      volume_animation(volume)%quads(quad,3)=point3
      volume_animation(volume)%quads(quad,4)=point4
    end do
    
    volume_animation(volume)%max_data=-1e30
    volume_animation(volume)%min_data=1e30
    
! read data    
    write(*,*)'Reading data'
    do quad=1,n_cells
    
      read(local_file_unit,*)cx,cy,cz,re,im
      
      volume_animation(volume)%complex_data(quad)=cmplx(re,im)
      
    end do
    
    write(*,*)'Finished reading fieldsolve geometry data'
      
  end do ! read next volume
 
  close (local_file_unit)
 
  write(*,*)' '
  write(*,*)'Number of volumes read:',n_volumes
  write(*,*)' '
  
2000 CONTINUE
! get the volume number to write to file

  write(*,*)'Enter the number of the frequency output volume to write to file or 0 to end'
  
  read(*,*)volume
  write(record_user_inputs_unit,*)volume,'   ! frequency output volume to write to file'
  
  if ( (volume.GT.0).AND.(volume.LE.n_volumes) ) then
! write the data to file

! open a file for the re-formatted data
  
    write(*,*)'Enter the filename for the formatted data'
    read(*,'(A)')filename
    write(record_user_inputs_unit,'(A)')trim(filename)
    
    write(*,*)'Do you want to append new data to an existing data file? (y/n)'
    read(*,'(A)')ch
    write(record_user_inputs_unit,'(A,A)')ch,' ! Append to existing file?'
  
    if ( (ch.EQ.'y').OR.(ch.EQ.'Y') ) then

! open file and read to end of file if it exists
      inquire(file=trim(filename),exist=file_exists)
      if (file_exists) then
    
        OPEN(unit=local_file_unit,file=filename,status='OLD',access='APPEND',position='APPEND')
      
      else
    
! file does not exist so open new file    
        OPEN(unit=local_file_unit,file=filename)
 
      end if

    else  

! open new file 
      OPEN(unit=local_file_unit,file=filename)
    
    end if

    n_points=volume_animation(volume)%n_points
    n_quads=volume_animation(volume)%n_quads
    n_frames=volume_animation(volume)%n_frames
    n_cells=n_points/8

! loop over the output cells  
    do quad=1,n_cells
    
      re=real( volume_animation(volume)%complex_data(quad) )
      im=imag( volume_animation(volume)%complex_data(quad) )

! Get the cell centre coordinate of the cell as the average of the 8 corner point coordinates
      point=(quad-1)*8+1
      x=0d0
      y=0d0
      z=0d0
      do i=point,point+7
        x=x+volume_animation(volume)%points(i,1)
        y=y+volume_animation(volume)%points(i,2)
        z=z+volume_animation(volume)%points(i,3)
      end do
      
      x=x/8d0
      y=y/8d0
      z=z/8d0

       write(local_file_unit,8000)x,y,z,re,im
       
8000  format(5E16.7)
  
    end do ! next cell
 
    if (compression_flag) then
     CALL close_output_file(local_file_unit,filename2,compression_flag)
    else
      CLOSE(unit=local_file_unit)
    end if
  
    GOTO 2000 ! see whether we need to output another of the volumes
 
  end if
 
! Deallocate memory
    
! Deallocate memory for animation  
    
  do volume=1,n_volumes
    if ( allocated( volume_animation(volume)%points ) ) DEALLOCATE ( volume_animation(volume)%points ) 
    if ( allocated( volume_animation(volume)%quads ) ) DEALLOCATE ( volume_animation(volume)%quads ) 
    if ( allocated( volume_animation(volume)%complex_data ) ) DEALLOCATE ( volume_animation(volume)%complex_data ) 
    if ( allocated( volume_animation(volume)%magnitude_data ) ) DEALLOCATE ( volume_animation(volume)%magnitude_data ) 
  end do
    
  if ( allocated( volume_animation ) ) DEALLOCATE ( volume_animation )  

  RETURN
  
END SUBROUTINE Convert_Frequency_Domain_Volume_Field_Format
