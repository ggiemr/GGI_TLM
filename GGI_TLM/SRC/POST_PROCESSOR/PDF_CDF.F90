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
! SUBROUTINE PDF_CDF
!
! NAME
!    PDF_CDF
!
! DESCRIPTION
!     
!     
! COMMENTS
!     Insertion sorting algorithm now replaced by heapsort which is much faster for large data sets
!
! HISTORY
!
!     started 8/08/2012 CJS
!     CDF calculation based on sorting the sample data 15/7/2014
!
SUBROUTINE PDF_CDF()

USE post_process
USE file_information

IMPLICIT NONE

! local variables

  character(len=256)	:: filename

  integer		:: sample,n_samples
  
  logical		:: file_exists
	   
  integer   		:: n_ignore
  integer   		:: n_function

  real*8,allocatable	:: read_data(:)
  real*8,allocatable	:: f(:)
  
  real*8		::fmax,fmin
 
  integer 		:: nbins
  integer		:: bx1
       
  integer,allocatable 	:: binx1(:)
  real*8,allocatable 	:: x(:)
  real*8,allocatable 	:: pdfx(:)
  real*8 dx
  
  real*8 		:: sum
  real*8 		:: mean,mu,variance,sigma
  
  integer 		:: loop,i
  
  character		:: ch
  integer		:: conversion_type,max_conversion_type
! START

! get the filename for the data
5 CONTINUE
  write(*,*)'Enter the filename for the input data'
  read(*,'(A256)')filename
  inquire(file=trim(filename),exist=file_exists)
  if (.NOT.file_exists) then
    write(*,*)'Error file does not exist'
    GOTO 5
  end if
  write(record_user_inputs_unit,'(A)')trim(filename)
      
  write(post_process_info_unit,*)'	Data filename:',trim(filename)

! open and read the file
  
  OPEN(unit=local_file_unit,file=filename)
  
  write(*,*)'Enter the number of lines to ignore at the top of the data file'
  read(*,*)n_ignore
  write(record_user_inputs_unit,*)n_ignore,' number of lines to ignore at the top of the data file'
    
  write(*,*)'Enter the column number for function data'
  read(*,*)n_function
  write(record_user_inputs_unit,*)n_function,' column for function data'
  write(post_process_info_unit,*)'	Process comun number:',n_function
  
  conversion_type=0
  
  write(*,*)'Is data conversion required? (y/n)'
  read(*,'(A)')ch
  write(record_user_inputs_unit,'(A,A)')ch,' ! Is data conversion required?'
  
  if ( (ch.EQ.'y').OR.(ch.EQ.'Y') ) then
    
    write(*,*)'Enter the conversion type:'
    write(*,*)'0: No conversion'
    write(*,*)'1: Convert from dB to field amplitude'
    write(*,*)'2: Convert from field amplitude to dB'
    
    max_conversion_type=2

    read(*,*)conversion_type
    if ( (conversion_type.LT.0).OR.(conversion_type.GT.max_conversion_type) ) then 
      write(*,*)'Conversion type should be in the range 1:',max_conversion_type
      STOP
    end if

    write(record_user_inputs_unit,*)conversion_type,' data conversion type'
  
    if (conversion_type.EQ.0) then
      write(post_process_info_unit,*)'	No data conversion is applied'
    else if (conversion_type.EQ.1) then
      write(post_process_info_unit,*)'	Convert from dB to field amplitude' 
    else if (conversion_type.EQ.2) then
      write(post_process_info_unit,*)'	Convert from field amplitude to dB' 
    end if

  end if
    
  ALLOCATE( read_data(1:n_function) )
       
  do loop=1,2

! read lines to ignore
    do i=1,n_ignore
      read(local_file_unit,*,end=1000,ERR=9000)
    end do

    sample=0
    
10  continue

      read(local_file_unit,*,end=1000,ERR=9000)(read_data(i),i=1,n_function)
      sample=sample+1
      if (loop.eq.2) then  
      
        if (conversion_type.eq.0) then
! No data conversion required	
    	  f(sample)=read_data(n_function)
	 
	else if (conversion_type.eq.1) then
! Convert from 	dB to field amplitude
    	  f(sample)=10d0**(read_data(n_function)/20d0)
	  
	else if (conversion_type.eq.2) then
! Convert from field amplitude to dB	 
          f(sample)=20d0*log10(abs(read_data(n_function)))
	  
        end if
	
      end if
     
    goto 10  ! read next sample
    
1000 continue 

    n_samples=sample
     
    if (loop.eq.1) then
      allocate ( f(1:n_samples) )
    end if
     
    rewind(unit=local_file_unit)
     
  end do ! next loop
      
  close(UNIT=local_file_unit)
      
  DEALLOCATE( read_data )
      
  write(*,*)'Number of data samples read:',n_samples

  fmax=f(1)
  fmin=f(1)
  
  do sample=2,n_samples
    fmax=max(fmax,f(sample))
    fmin=min(fmin,f(sample))
  end do
  
  write(*,*)'Minimum data value=',fmin
  write(*,*)'Maximum data value=',fmax
  
! calculate the mean and variance of the data set

! calculate mean
  mean=0.0       
  do sample=1,n_samples
    mean=mean+f(sample)
  end do
  mean=mean/n_samples
  write(*,2010)'mean f=     ',mean

! calculate variance       
  variance=0.0       
  do sample=1,n_samples
    variance=variance+((f(sample)-mean)**2)
  end do
  variance=variance/n_samples
  write(*,2010)'variance f= ',variance
  write(*,2010)'standard deviation f = ',sqrt(variance)

! get the number of bins for the distribution data

  write(*,*)'Enter the number of bins for the distribution data'
  read(*,*)nbins
  write(record_user_inputs_unit,*)nbins,' number of bins for the distribution data'

! open a file for the distribution data
  
  write(*,*)'Enter the filename for the distribution data (PDF and CDF)'
  read(*,'(A256)')filename
  write(record_user_inputs_unit,'(A)')trim(filename)
  
! open the file to write
  
  OPEN(unit=local_file_unit,file=filename)
  
! calculate the PDF and the CDF  og binned data
  
  ALLOCATE ( x(1:nbins) )
  ALLOCATE ( binx1(1:nbins) )
  ALLOCATE ( pdfx(1:nbins) )
  
  do i=1,nbins
    binx1(i)=0
  end do

  do sample=1,n_samples       
    
    bx1=int(dble(nbins)*(f(sample)-fmin)/(fmax-fmin))+1
    if ((bx1.ge.1).and.(bx1.le.nbins))then
      binx1(bx1)=binx1(bx1)+1
    else
      write(*,*)'Bin out of range',bx1,nbins
    end if
   	    
  end do
  
  dx=(fmax-fmin)/nbins
  write(*,*)'dx=',dx
  
  do i=1,nbins
    x(i)=fmin+dx*(dble(i)-0.5d0)
    pdfx(i)=real(binx1(i))/(dble(n_samples)*dx)
  end do
  
! integrate the PDF
  sum=0
  write(local_file_unit,2000)0,fmin,0d0,0d0
  do i=1,nbins  	    
    sum=sum+pdfx(i)*dx/2d0
    write(local_file_unit,2000)i-1,x(i),pdfx(i),sum
    sum=sum+pdfx(i)*dx/2d0
  end do
  write(local_file_unit,2000)nbins,fmax,0d0,1d0
2000   format(I6,3E16.8)	 

  close(unit=local_file_unit)

  write(*,*)' '
  write(*,*)'integral pdf(x)dx=',sum   
           
! calculate mean
  mean=0.0       
  do i=1,nbins
    mean=mean+pdfx(i)*x(i)*dx
  end do
  write(*,2010)'mean f=     ',mean
2010 format(A,E14.6)

! calculate variance       
  variance=0.0       
  do i=1,nbins
    variance=variance+((x(i)-mean)**2)*pdfx(i)*dx
  end do
  write(*,2010)'variance f= ',variance
  write(*,2010)'sigma f   = ',sqrt(variance)
       
  DEALLOCATE ( x )
  DEALLOCATE ( binx1 )
  DEALLOCATE ( pdfx )
  
! NEW PROCESS TO GIVE HIGHER RESOLUTION RESULTS FOR PDF

  nbins=n_samples

! open a file for the distribution data
  
  write(*,*)'Enter the filename for the distribution data : CDF only (new process)'
  read(*,'(A256)')filename
  write(record_user_inputs_unit,'(A)')trim(filename)
  
! open the file to write
  
  OPEN(unit=local_file_unit,file=filename)
  
! calculate the PDF by sorting of the sample data

!  CALL sample_sort_insertion(f,n_samples)
  CALL sample_sort_heapsort(f,n_samples)
  
  do sample=1,n_samples       
    
    write(local_file_unit,*)f(sample),real(sample)/real(n_samples),1d0-real(sample)/real(n_samples)

  end do
  
  close(unit=local_file_unit)

  RETURN
     
9000 CALL write_line('Error reading data file',0,.TRUE.)
     CALL write_line('filename:',0,.TRUE.)
     write(*,*)trim(filename)
     STOP 
  
END SUBROUTINE PDF_CDF
!
! ___________________________________________________________________________
!
!
 SUBROUTINE sample_sort_insertion(f,n)
 
IMPLICIT NONE 
 
  integer n
  real*8 f(1:n)
  real*8 fj
 
  integer i,j
 
! START
 
  do j=2,n
  
    fj=f(j)
    i=j-1
    
    do while ( (i.GT.0).AND.(f(i).GT.fj) )
    
      f(i+1)=f(i)
      i=i-1 
      
    end do
  
    f(i+1)=fj
    
  end do
 
  RETURN
 
 END SUBROUTINE sample_sort_insertion
 
!
! ___________________________________________________________________________
!
!
  SUBROUTINE sample_sort_heapsort(f,n)
 
  IMPLICIT NONE 
 
  integer n
  real*8 f(1:n)
 
  integer i,end
  
  integer j
 
! START

  CALL build_heap(f,n)
    
  end=n
  
!   write(*,*)'START'
!  CALL check_heap(f,n,end)
  
  do i=1,n-1
   
! the maximum value is at the head of the heap. 
! Swap this with the end value in the heap and 
! reduce the size of the heap by 1

    CALL swap(f,n,1,end)
    
    end=end-1

! the new value at the top of the heap now needs to be filtered down to it's correct level
! this process ensures that the heap property is maintained
    CALL sift_down(f,n,end,1)
  
!    write(*,*)'STAGE',i
!    CALL check_heap(f,n,end)
     
  end do
 
  RETURN
 
  END SUBROUTINE sample_sort_heapsort
!
! ___________________________________________________________________________
!
!
  SUBROUTINE check_heap(f,n,end)
 
  IMPLICIT NONE 
 
  integer n,end
  real*8 f(1:n)
  
  integer i,j,parent
 
! START

  do i=2,end
  
    parent=i/2
  
    if (f(parent).LT.f(i)) then
      write(*,*)'Heap failed, i=',i,' i/2=',i/2
      write(*,*)'f(i/2)=',f(i/2)
      write(*,*)'f(i)=',f(i)
      
      write(*,*)
 
      do j=1,n
  
         write(*,*)j,f(j),j/2

      end do
      
      STOP
    end if
  
  end do
 
  RETURN
 
  END SUBROUTINE check_heap
!
! ___________________________________________________________________________
!
!
  SUBROUTINE build_heap(f,n)
 
  IMPLICIT NONE 
 
  integer n
  real*8 f(1:n)
  
  integer i
 
! START

  do i=2,n
  
    CALL sift_up(f,n,i)
  
  end do
 
  RETURN
 
  END SUBROUTINE build_heap
!
! ___________________________________________________________________________
!
!
  RECURSIVE SUBROUTINE sift_up(f,n,i)
 
IMPLICIT NONE 
 
  integer n,i
  real*8 f(1:n)
  
  integer parent
 
! START

  if (i.eq.1) RETURN
  
  parent=i/2
  
  if (f(i).GT.f(parent)) then
    CALL swap(f,n,i,parent)
    CALL sift_up(f,n,parent)
  else
    RETURN
  end if
 
  END SUBROUTINE sift_up
!
! ___________________________________________________________________________
!
!
  RECURSIVE SUBROUTINE sift_down(f,n,end,parent)
 
  IMPLICIT NONE 
 
  integer n,end,parent
  real*8 f(1:n)
  
  integer child1,child2,cmax
  real*8 fcmax
 
! START

!  if (parent.GE.end) RETURN   ! Check not now required due to checks on the existance of children
  
  child1=parent*2
  
  if (child1.GT.end) RETURN  ! no children in the heap
  
  child2=child1+1
  
! work out the maximum child value which exists in the heap
  cmax=child1
  fcmax=f(child1)
  if ( (child2.LE.end).AND.(f(child2).GT.fcmax) ) then
    cmax=child2
    fcmax=f(child2)
  end if
  
! check whether the parent value is larger than the maximum child value and of so continue to sift_down
  if (fcmax.GT.f(parent)) then
  
    CALL swap(f,n,cmax,parent)
    CALL sift_down(f,n,end,cmax)
    
  end if
  
  END SUBROUTINE sift_down
!
! ___________________________________________________________________________
!
!
  SUBROUTINE swap(f,n,i,j)
 
  IMPLICIT NONE 
 
  integer n,i,j
  real*8 f(1:n)
    
  real*8 fswap
 
! START
    
  fswap=f(i)
  f(i)=f(j)
  f(j)=fswap
 
  RETURN
 
  END SUBROUTINE swap

