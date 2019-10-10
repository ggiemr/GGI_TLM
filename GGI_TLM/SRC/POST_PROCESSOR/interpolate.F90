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
! SUBROUTINE interpolate
!
! NAME
!    interpolate
!
! DESCRIPTION
!    read a data file with a set of x values x1 or alternatively set a 
!    uniform set of n_samples sampling points x1 from x1min to x1max
!    read a second file with a set of real or complex function values f(x) at a set of x values x2
!    interpolate the values f(x) to the set of values x1
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 1/03/2018 CJS
!
!
SUBROUTINE interpolate


USE post_process
USE file_information

IMPLICIT NONE

! local variables

  real*8	:: interpolate_value
  
  real*8	:: x1min,x1max
  real*8	:: x2min,x2max
  real*8	:: xmin,xmax
  real*8	:: opxmin,opxmax
  real*8	:: x
  
  integer	:: n_samples1,n_samples2
  integer       :: n_ignore
  integer       :: x1_col,x2_col,re_col,im_col,max_col
  
  real*8,allocatable :: data_line(:)
  real*8        :: xscale
  
  integer	:: xloop
  integer	:: loop
  
 
  character(len=256)	:: x1filename
  character(len=256)	:: x2filename
  character(len=256)	:: opfilename
  
  integer	:: sample_min,sample_max,last_sample,sample
  
  real*8,allocatable	:: x1(:)
  real*8,allocatable	:: x2(:)
  real*8,allocatable	:: re(:)
  real*8,allocatable	:: im(:)
  
  complex*16	:: value_1,value_2
         
  complex*16	:: interpolated_value
    
  integer	:: i
  
! START
  
  write(*,*)' '

  write(*,*)'Please enter the filename with required x values (x1) '
  write(*,*)'or press enter to define the x1 samples from n_x1_samples, x1min and x1max'
  
  read(*,'(A)')x1filename
  
  write(record_user_inputs_unit,'(A)')trim(x1filename)
  
  if (len(trim(x1filename)).EQ.0) then

! set the x1 values from start, end and step values
    write(*,*)'Please enter the number of x1 data samples'  
    read(*,*)n_samples1
    
    if (n_samples1.LT.1) then
      write(*,*)'ERROR in interpolate: number of samples should be at least 1'
      STOP 1    
    end if
    
    write(record_user_inputs_unit,'(I12,A33)')n_samples1,'   ! number of x1 samples '

    write(*,*)'Please enter the minimum x1 value'  
    read(*,*)x1min
    write(record_user_inputs_unit,'(ES16.6,A33)')x1min,'   ! minimum value for x1 data  '

    write(*,*)'Please enter the maximum x1 value'  
    read(*,*)x1max
    write(record_user_inputs_unit,'(ES16.6,A33)')x1max,'   ! maximum value for x1 data  '
    
    write(*,*)'Number of x1 samples:',n_samples1
    write(*,*)'Minimum x1 value is: ',x1min
    write(*,*)'Maximum x1 value is: ',x1max

    ALLOCATE( x1(1:n_samples1) )
    
    if (n_samples1.GT.1) then
      do i=1,n_samples1
        x1(i)=x1min+dble(i-1)*(x1max-x1min)/dble(n_samples1-1)
      end do
    else if (n_samples1.EQ.1) then
      x1(1)=x1min
    end if
    
  else
  
! read x1 values from a file
  
    write(post_process_info_unit,*)'Filename for x1 data:',trim(x1filename)
  
    OPEN(unit=local_file_unit,file=x1filename,err=9000)

! read the x1 values only

    write(*,*)'Please enter the number of lines to ignore at the top of the x1 file'  
    read(*,*)n_ignore
    write(record_user_inputs_unit,'(I4,A32)')n_ignore,'   ! number of lines to ignore  '

    write(*,*)'Please enter the column for x1 data'  
    read(*,*)x1_col
    write(record_user_inputs_unit,'(I4,A25)')x1_col,'   ! column for x1 data  '
  
! Allocate the array for reading the data file
    max_col=x1_col
    ALLOCATE( data_line(1:max_col) )
  
    xscale=1d0
    do loop=1,2      ! read loop, the first to count the samples and the second to read them
  
      do i=1,n_ignore
        read(local_file_unit,*)
      end do
    
      n_samples1=0
      x1min=1D30
      x1max=-1D30

100   CONTINUE

      read(local_file_unit,*,end=110,err=120)(data_line(i),i=1,max_col)
      n_samples1=n_samples1+1
    
      x1min=min(x1min,data_line(x1_col)*xscale)
      x1max=max(x1max,data_line(x1_col)*xscale)
      
      if (loop.EQ.2) then
    
        x1(n_samples1)=data_line(x1_col)*xscale
      
      end if
    
120   CONTINUE
    
      GOTO 100

110   CONTINUE ! jump here when the file has been read
    
      write(*,*)'Loop:',loop
      write(*,*)'Number of samples read:',n_samples1
      write(*,*)'Minimum x1 value is: ',x1min
      write(*,*)'Maximum x1 value is: ',x1max

      if (loop.EQ.1) then

        write(*,*)'Please enter the scaling factor for x1 data'  
        read(*,*)xscale
        write(record_user_inputs_unit,'(ES16.6,A33)')xscale,'   ! scaling factor for x1 data  '
      
        ALLOCATE( x1(1:n_samples1) )
        rewind(local_file_unit)
      
      end if
 
    end do ! next read loop
  
    DEALLOCATE( data_line )

    CLOSE(unit=local_file_unit)
    
  end if
  
! we now have n_samples1 time samples in the array x1
  
! Read the existing x values (x2) and function values, f(x2)

  write(*,*)'Please enter the filename for x=x2 values and function data at x2'
  read(*,*)x2filename
  write(record_user_inputs_unit,'(A)')trim(x2filename)
  write(post_process_info_unit,*)'Filename for function data at x2:',trim(x2filename)
  
! read the x2 values and function values
  
  OPEN(unit=local_file_unit,file=x2filename,err=9010)

  write(*,*)'Please enter the number of lines to ignore at the top of the x2 file'  
  read(*,*)n_ignore
  write(record_user_inputs_unit,'(I4,A32)')n_ignore,'   ! number of lines to ignore  '

  write(*,*)'Please enter the column for x2 data'  
  read(*,*)x2_col
  write(record_user_inputs_unit,'(I4,A25)')x2_col,'   ! column for x2 data  '

  write(*,*)'Please enter the column for real part data'  
  read(*,*)re_col
  write(record_user_inputs_unit,'(I4,A25)')re_col,'   ! column for Re data  '

  write(*,*)'Please enter the column for imaginary part data or 0 if there is none'  
  read(*,*)im_col
  write(record_user_inputs_unit,'(I4,A25)')im_col,'   ! column for Im data  '
  
! Allocate the array for reading the data file
  max_col=max(x2_col,re_col,im_col)
  ALLOCATE( data_line(1:max_col) )
  
  xscale=1d0
  do loop=1,2      ! read loop, the first to count the samples and the second to read them
  
    do i=1,n_ignore
      read(local_file_unit,*)
    end do
    
    n_samples2=0
    x2min=1D30
    x2max=-1D30
    
200 CONTINUE

    read(local_file_unit,*,end=210,err=220)(data_line(i),i=1,max_col)
    n_samples2=n_samples2+1
    
    x2min=min(x2min,data_line(x2_col)*xscale)
    x2max=max(x2max,data_line(x2_col)*xscale)
      
    if (loop.EQ.2) then
    
      x2(n_samples2)=data_line(x2_col)*xscale
      
      re(n_samples2)=data_line(re_col)
      
      if (im_col.NE.0) then
        im(n_samples2)=data_line(im_col)
      else
        im(n_samples2)=0d0
      end if
      
    end if

220 CONTINUE 
    
    GOTO 200

210 CONTINUE ! jump here when the file has been read

    write(*,*)'Loop:',loop
    write(*,*)'Number of samples read:',n_samples2
    write(*,*)'Minimum x2 value is: ',x2min
    write(*,*)'Maximum x2 value is: ',x2max
    
    if (loop.EQ.1) then

      write(*,*)'Please enter the scaling factor for x2 data'  
      read(*,*)xscale
      write(record_user_inputs_unit,'(ES18.8,A33)')xscale,'   ! scaling factor for x2 data  '
      
      ALLOCATE( x2(1:n_samples2) )
      ALLOCATE( re(1:n_samples2) )
      ALLOCATE( im(1:n_samples2) )
      rewind(local_file_unit)
      
    end if
 
  end do ! next read loop
  
  DEALLOCATE( data_line )

  CLOSE(unit=local_file_unit)
  
  xmin=max(x1min,x2min)
  xmax=min(x1max,x2max)
  
  write(*,*)'Overlapping x range is ',xmin,' to ',xmax
  
  write(*,*)'Please enter xmin for the output data'  
  read(*,*)opxmin
  write(record_user_inputs_unit,'(ES18.8,A24)')opxmin,'   ! output xmin value  '
  
  write(*,*)'Please enter xmax for the output data'  
  read(*,*)opxmax
  write(record_user_inputs_unit,'(ES18.8,A24)')opxmax,'   ! output xmax value  '

  write(*,*)'Please enter the filename for the output data'
  
  read(*,*)opfilename
  write(record_user_inputs_unit,'(A)')trim(opfilename)
  write(post_process_info_unit,*)'Filename for output data:',trim(opfilename)
  
  OPEN(unit=local_file_unit,file=opfilename)

! loop over samples of file 1 to work out sample_min and sample_max 

  sample_min=0
  sample_max=0
  
  do sample=1,n_samples1
  
    if ( ( x1(sample).ge.opxmin ).AND.( sample_min.eq.0 ) ) sample_min=sample
    if ( ( x1(sample).ge.opxmax ).AND.( sample_max.eq.0 ) ) sample_max=sample
    
  end do
  
  if (sample_min.eq.0) sample_min=1
  if (sample_max.eq.0) sample_max=n_samples1
  
  write(*,*)'n_samples1 =',n_samples1
  write(*,*)'Sample_min=',sample_min,' xmin=',x1(sample_min)
  write(*,*)'Sample_max=',sample_max,' xmax=',x1(sample_max)
       
! loop over samples in the x range and interpolate values

  last_sample=1

  do sample=sample_min,sample_max
  
    x=x1(sample)
    
! find the x values which lie either side of x and interpolate to give the function value from file 2
    do i=last_sample,n_samples2-1
    	 
      if ( (x2(i).le.x).AND.	&
           (x2(i+1).gt.x)  ) then
           
        value_1=cmplx(re(i),im(i))
        value_2=cmplx(re(i+1),im(i+1))
	   
    	interpolated_value=value_1+( (value_2-value_1)/(x2(i+1)-x2(i)) )*(x-x2(i))
	    
    	last_sample=i
    	
    	GOTO 1000
	
      end if
      
    end do  ! next sample of file 2

    write(*,*)'Sample not found in file 2, x=',x
    write(*,*)'First x=',x2(last_sample)
    write(*,*)'Last x=',x2(n_samples2-1)

    STOP

1000  CONTINUE
  
  write(local_file_unit,'(ES20.10,2ES16.6)')x,real(interpolated_value),imag(interpolated_value)

  end do ! next x1 sample

  CLOSE(local_file_unit)
   
  DEALLOCATE( x1 )
  DEALLOCATE( x2 )
  DEALLOCATE( re )
  DEALLOCATE( im )
  
  RETURN
  
9000 write(*,*)'ERROR in interpolate: cannot open the file:',trim(x1filename)
     STOP
  
9010 write(*,*)'ERROR in interpolate: cannot open the file:',trim(x2filename)
     STOP
      
END SUBROUTINE interpolate
