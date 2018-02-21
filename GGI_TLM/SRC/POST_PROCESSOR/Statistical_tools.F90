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
! SUBROUTINE Statistical_tools
!
! NAME
!    Statistical_tools
!
! DESCRIPTION
!     
!     
! COMMENTS
!     Some statistical tools for investigating 
!
! HISTORY
!
!     Some statistical tools for investigating 28/11/2014
!
SUBROUTINE Statistical_tools()

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
  
  real*8		:: fmax,fmin
 
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
    write(*,*)'3: Convert from field to field magnitude'
    
    max_conversion_type=3

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
    else if (conversion_type.EQ.3) then
      write(post_process_info_unit,*)'	Convert from field to field magnitude' 
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
	  
	else if (conversion_type.eq.3) then
! Convert from field to field magnitude	 
          f(sample)=(abs(read_data(n_function)))
	  
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
  
  fmin=fmin-abs(fmin)*0.000001d0
  fmax=fmax+abs(fmax)*0.000001d0


  write(*,*)'Enter the number of bins for the distribution data'
  read(*,*)nbins
  write(record_user_inputs_unit,*)nbins,' number of bins for the distribution data'

! open a file for the distribution data
  
  write(*,*)'Enter the filename for the distribution data (PDF and CDF)'
  read(*,'(A256)')filename
  write(record_user_inputs_unit,'(A)')trim(filename)
  
! open the file to write
  
  OPEN(unit=local_file_unit,file=filename)
  
! calculate the PDF and the CDF 1:of binned data
  
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
  write(*,2010),'mean f=     ',mean
2010 format(A,E14.6)

! calculate variance       
  variance=0.0       
  do i=1,nbins
    variance=variance+((x(i)-mean)**2)*pdfx(i)*dx
  end do
  write(*,2010),'variance f= ',variance
  write(*,2010),'sigma f   = ',sqrt(variance)
       
  DEALLOCATE ( x )
  DEALLOCATE ( binx1 )
  DEALLOCATE ( pdfx )

  RETURN
     
9000 CALL write_line('Error reading data file',0,.TRUE.)
     CALL write_line('filename:',0,.TRUE.)
     write(*,*)trim(filename)
     STOP 
      
END SUBROUTINE Statistical_tools
