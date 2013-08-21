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
       program near_field_convolution

USE constants

IMPLICIT NONE
       
       
       real*8,allocatable :: ft1(:),ft_theta(:),ft_phi(:),t1(:),f1(:)
       complex*16,allocatable :: ff1(:),ff_theta(:),ff_phi(:)
       complex*16,allocatable :: ff_theta2(:),ff_phi2(:)
       real*8 dt1,time
       real*8 df
       real*8 t_in,ft_in
       real*8 rs
       integer n_in
       integer loop
      
       real*8 fmin,fmax,fstep
       real*8 fscale,tscale
       integer nt
       
       character*256 ipfile1
       character*256 ipfile2
       character*256 opfile1
       character*256 opfile2
       integer op_number1
       integer tloop
       integer sample,nsamples,sample1
       integer fsample,nfsamples
       real*8 f,w
       real*8 f_in,theta,phi,rcs,et1,et2,ep1,ep2
       
       real*8 time_shift
       
       complex*16 Etheta,Ephi
       complex*16 integral1,integral2,ftc,cout,ftc1,ftc2
       
! START       
       
! read input file     
       
       open(unit=99,file='progress')
       write(99,'(A)')'STARTED: GGI_near_field_convolution'
       close(unit=99)

       print*,'Enter the time domain source input file'
       read(*,'(A256)'),ipfile1       
       open(unit=10,file=ipfile1)
       print*,'Enter the output number to transform'
       read(*,*),op_number1     
       
       print*,'Enter the sphere near field filename (from RCS_sphere code)'
       read(*,'(A256)'),ipfile2       
       open(unit=11,file=ipfile2)
       
! read ipfile1             

       do loop=1,2

         sample=0
       
10       continue

           read(10,*,end=1000)t_in,n_in,ft_in
           if (n_in.eq.op_number1) then ! use this output
	     sample=sample+1
	     if (loop.eq.2) then  
	       t1(sample)=t_in
	       ft1(sample)=ft_in
	     end if
	   end if
	 
         goto 10  ! read next sample
        
1000     continue 
   
         nsamples=sample
	 
         if (loop.eq.1) then
	   allocate ( t1(1:nsamples) )
	   allocate ( ft1(1:nsamples) )
	   allocate ( ft_theta(1:nsamples) )
	   allocate ( ft_phi(1:nsamples) )
	 end if
	 
	 rewind(unit=10)
	 
       end do ! next loop
             
       dt1=t1(2)-t1(1)
       print*,nsamples,' samples read from input file1'
       print*,'dt=',dt1
                     
       close(unit=10)
       
       do loop=1,2

         sample=0
       
20       continue

           read(11,*,end=2000)f_in,theta,phi,rcs,et1,et2,ep1,ep2
	   sample=sample+1
	   if (loop.eq.2) then  
	     f1(sample)=f_in
	     ff_theta(sample)=cmplx(et1,et2)
	     ff_phi(sample)=cmplx(ep1,ep2)
	   end if
	 
         goto 20  ! read next sample
        
2000     continue 
   
         nfsamples=sample
	 
         if (loop.eq.1) then
	   allocate ( f1(1:nfsamples) )
	   allocate ( ff1(1:nfsamples) )
	   allocate ( ff_theta(1:nfsamples) )
	   allocate ( ff_phi(1:nfsamples) )
	   allocate ( ff_theta2(1:nfsamples) )
	   allocate ( ff_phi2(1:nfsamples) )
	 end if
	 
	 rewind(unit=11)
	 
       end do ! next loop
       
       close(unit=11)
              
       print*,' '
       print*,'Enter the file name for the fourier transformed pulse'
       read(*,'(A256)'),opfile1       
       open(unit=12,file=opfile1)
       print*,' '
       print*,'Enter the output file name for the time domain near field'
       read(*,'(A256)'),opfile2       
       open(unit=13,file=opfile2)
       print*,'Enter the time shift to be applied to the time domain near field'
       read(*,*),time_shift 

! transform pulse to frequency domain and calculate 
! response in the frequency domain
! loop over frequency
       do fsample=1,nfsamples
       
         w=2.0*pi*f1(fsample)
	 
         integral1=(0.0,0.0)
	 
	 do sample=1,nsamples
	   ftc=exp(-j*w*t1(sample))*dt1
	   integral1=integral1+ft1(sample)*ftc
	 end do
	 	 
	 ff1(fsample)=integral1
	 
	 ff_theta2(fsample)=ff_theta(fsample)*integral1
	 ff_phi2(fsample)=ff_phi(fsample)*integral1
	 		 
	 cout=ff1(fsample)
	 cout=ff_theta2(fsample)
	 
	 write(12,8000)f1(fsample),real(cout),aimag(cout),abs(cout)
8000     format(4E16.8)

       end do
      
       close(unit=12)
       
! Transform back to the time domain and write the output pulse      
! loop over time

       do sample=1,nsamples
!       do sample=-nsamples,nsamples
 
         if (sample.lt.0) then
	   time=-t1(-sample)
	 else if (sample.eq.0) then
	   time=0.0
	 else 
	   time=t1(sample)
	 end if
 
         integral1=(0.0,0.0)
         integral2=(0.0,0.0)
	 
         do fsample=1,nfsamples
	         
           if (fsample.eq.1) then
! special case for close to d.c.
             df=(f1(1)+f1(2))/2.0
           else
             df=(f1(2)-f1(1))
	   end if
	   
           w=2.0*pi*f1(fsample)
	   ftc1=exp(j*w*time)*df
!	   ftc2=exp(-j*w*time)*df
	   
! note: add negative frequencies...
	   integral1=integral1+ff_theta2(fsample)*ftc1
	   integral1=integral1+conjg(ff_theta2(fsample)*ftc1)
	   
	   integral2=integral2+ff_phi2(fsample)*ftc1
	   integral2=integral2+conjg(ff_phi2(fsample))*ftc2
	   
	 end do
		 
	 cout=integral1
	 
	 write(13,8010)time+time_shift ,real(integral1),imag(integral1)
8010     format(3E16.8)

       end do
       
       close(unit=13)
       
       deallocate ( t1 )
       deallocate ( f1 )
       deallocate ( ft1 )
       deallocate ( ff1 )
       deallocate ( ft_theta )
       deallocate ( ft_phi )
       deallocate ( ff_theta )
       deallocate ( ff_phi )
       deallocate ( ff_theta2 )
       deallocate ( ff_phi2 )
       
       open(unit=99,file='progress')
       write(99,'(A)')'FINISHED: GGI_near_field_convolution'
       close(unit=99)
       
       end

