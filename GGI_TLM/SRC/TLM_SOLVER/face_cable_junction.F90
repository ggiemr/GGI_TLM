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
! SUBROUTINE face_cable_junction
! SUBROUTINE face_junction
!
! NAME
!     face_cable_junction
!
! DESCRIPTION
!
!     Cable update on a face - at the moment this is an interface to the Fieldsolve subroutine face_junction()
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 24/09/2012 CJS
!             17/12/2012 CJS	Started to implement junction loads
!             23/04/2013 CJS	Fix bug relating to P matrix dimensions
!
!
SUBROUTINE face_cable_junction(cable_face_number)

USE TLM_General
USE Cables

IMPLICIT NONE

  integer	:: cable_face_number

! local variables

  integer maxwires,nw		      
  integer N_int
  integer N_ext(3)  
       
  real*8,allocatable	::  P(:,:,:)
  real*8,allocatable	::  Tv(:,:,:)
  real*8,allocatable	::  Ti(:,:,:)
  
  real*8,allocatable	:: NOT_SC(:,:)
  integer,allocatable	:: BC(:)
  
  real*8,allocatable	:: Yf(:,:,:)
  real*8,allocatable	:: Zlink(:,:,:)
  real*8,allocatable	:: Vf(:,:)
  
  real*8,allocatable	:: Vw(:,:)
  real*8,allocatable	:: Iw(:,:)
  
  real*8,allocatable	:: Vlink(:,:)
  
  real*8,allocatable	:: Temp_vector(:)
  
  real*8 Itemp

  integer	:: face
  integer	:: segment
  integer	:: row,col
  integer	:: filter

! START

  CALL write_line('CALLED: face_cable_junction',0,timestepping_output_to_screen_flag)
  
  do face=1,3
  
    N_ext(face)=face_junction_list(cable_face_number)%n_external_conductors(face)
    
  end do ! next face
  
!!!!!!  n_ext(3)=0  ! this if there are no internal imedances 
 
  N_int=face_junction_list(cable_face_number)%n_internal_connection_nodes
  
  maxwires=max(N_ext(1),N_ext(2),N_ext(3),N_int)

! allocate arrays
  ALLOCATE( P(maxwires,maxwires,3) )
  ALLOCATE( Tv(maxwires,maxwires,3) )
  ALLOCATE( Ti(maxwires,maxwires,3) )
  ALLOCATE( NOT_SC(maxwires,3) )
  ALLOCATE( BC(maxwires) )
  ALLOCATE( Yf(maxwires,maxwires,3) )
  ALLOCATE( Zlink(maxwires,maxwires,3) )
  ALLOCATE( Vf(maxwires,3) )
  ALLOCATE( Vw(maxwires,3) )
  ALLOCATE( Iw(maxwires,3) )
  ALLOCATE( Temp_vector(maxwires) )

! reset arrays

  P(:,:,:)=0d0
  Tv(:,:,:)=0d0
  Ti(:,:,:)=0d0
  NOT_SC(:,:)=0d0
  BC(:)=0
  Yf(:,:,:)=0d0
  Zlink(:,:,:)=0d0
  Vf(:,:)=0d0
  Vw(:,:)=0d0
  Iw(:,:)=0d0
  
  do row=1,N_int
    BC(row)=face_junction_list(cable_face_number)%BC(row)
  end do

  do face=1,2
  
    segment=face_junction_list(cable_face_number)%segment_list(face)
    
    if (segment.ne.0) then
! there is a cable on this face so fill the appropriate data from the bundle segment data structure

      if (N_ext(face).NE.bundle_segment_list(segment)%n_conductors) then
        write(*,*)'Error in face_cable_junction'
	write(*,*)'N_external discrepancy, face=',face
	STOP
      end if

! Copy matrices and vectors with any required type conversion from integer to real*8      

      do row=1,N_ext(face)
      
        if (bundle_segment_list(segment)%SC(row).EQ.0) then
          NOT_SC(row,face)=1d0
        end if
	
        do col=1,N_ext(face)
	
          Tv(row,col,face)=dble(bundle_segment_list(segment)%Tv(row,col))
          Ti(row,col,face)=dble(bundle_segment_list(segment)%Ti(row,col))
          Zlink(row,col,face)=bundle_segment_list(segment)%Zlink(row,col)
          Yf(row,col,face)=bundle_segment_list(segment)%Ylink(row,col)
	
	end do ! next col
	
      end do ! next row
     
! Copy P matrix with required type conversion from integer to real*8      
      do col=1,N_ext(face)
        	
        do row=1,N_int
	
          P(row,col,face) =dble(face_junction_list(cable_face_number)%P_matrix_list(face)%P(row,col))
 		
	end do ! next col
	
      end do ! next row
      
! Work out the cable voltage vector, Vf including cable voltage source where required  
      
      do row=1,N_ext(face)
    
        Vf(row,face)=2d0*bundle_segment_list(segment)%Vlink(row)
		
      end do ! next row
      
    else ! segment=0 so no cable on this face

      N_ext(face)=0

    end if  ! segment=0?
    
  end do
  
! START of process for junction impedances.
  face=3
  
! Copy matrices and vectors with any required type conversion from integer to real*8      
  do row=1,N_ext(face)
    
    NOT_SC(row,face)=0d0

    do col=1,N_ext(face)

      if (row.ne.col) then
        Tv(row,col,face)=0d0
        Ti(row,col,face)=0d0
      else
        Tv(row,col,face)=1d0
        Ti(row,col,face)=1d0
      end if
      Yf(row,col,face)=face_junction_list(cable_face_number)%Yf(row,col)

    end do ! next col

  end do ! next row
      
! Copy P matrix with required type conversion from integer to real*8      
  do col=1,N_ext(face)
  	    
    do row=1,N_int

      P(row,col,face) =dble(face_junction_list(cable_face_number)%P_matrix_list(face)%P(row,col))
  	    
    end do ! next col

  end do ! next row
   	 
! Voltage sources related to the impedance filter functions  

! Evaluate Impedance filters
  do row=1,face_junction_list(cable_face_number)%n_internal_impedance_filters

    filter=row
    CALL timeshift_Zfilter(face_junction_list(cable_face_number)%Zfilter_data(filter))
    CALL evaluate_Zfilter (face_junction_list(cable_face_number)%Zfilter(filter),		&
    			   face_junction_list(cable_face_number)%Zfilter_data(filter), &
    			   0d0)
    Vf(row,face)=Vf(row,face)+face_junction_list(cable_face_number)%Zfilter_data(filter)%f

  end do
  
! END of internal junction impedance process

  CALL face_junction(cable_face_number,maxwires,N_int,N_ext,P,Tv,Ti,BC,NOT_SC,Yf,Vf,Vw,Iw  )

! calculate the voltage across the lines on both sides of the face
	   
  do face=1,2
	   
    segment=face_junction_list(cable_face_number)%segment_list(face)
    
    if (segment.ne.0) then
    
! there is a cable on this face so calculate the scattered voltages (Vr)=(Vi)-[Z]*(Iw)

      nw=N_ext(face)

      call dmatvmul(Zlink(1,1,face),nw,nw,Iw(1,face),nw,Temp_vector,maxwires)
      do row=1,nw
        bundle_segment_list(segment)%Vlink(row)=bundle_segment_list(segment)%Vlink(row)-Temp_vector(row)
      end do
      
! save the cable current vector
      do row=1,nw
        bundle_segment_list(segment)%Iw_face(row)=Iw(row,face)
      end do
	       	       
    end if ! there is a bundle segment on this face
	          
  end do ! next face
  
! Internal junction impedance process  
  face=3
 
! Update Internal Impedance filters
  do row=1,face_junction_list(cable_face_number)%n_internal_impedance_filters

    filter=row
    Itemp=-Iw(row,face)
    CALL evaluate_Zfilter (face_junction_list(cable_face_number)%Zfilter(filter),	   &
    			   face_junction_list(cable_face_number)%Zfilter_data(filter), &
    			   Itemp)

  end do
  

! deallocate arrays
  DEALLOCATE( P )
  DEALLOCATE( Tv )
  DEALLOCATE( Zlink )
  DEALLOCATE( NOT_SC )
  DEALLOCATE( BC )
  DEALLOCATE( Yf )
  DEALLOCATE( Vf )
  DEALLOCATE( Vw )
  DEALLOCATE( Iw )
  DEALLOCATE( Temp_vector )

 
  CALL write_line('FINISHED: face_cable_junction',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE face_cable_junction
!
! Name face_junction
!     
!
! Description
!     TLM MTL general TLM face junction formulation
!
! Comments:
!     Equation numbers refer to the GGI_TLM_cable_model_theory document
!
!     This subroutine has been taken directly from Fieldsolve and needs to be re-written
!
! History
!
!     started 4/10/09 CJS
!     termination internal impedances included 24/02/2011 CJS
!     made TLM specific for testing, original saved as: face_junction.F90_saved_25_2_2011
!
       SUBROUTINE face_junction(cell,maxwires,Nint_in,Next,Pin,Tvin,Tiin,BC,SC,Yf,Vf,Vw,Iw  )
USE cell_parameters		        
USE file_information

IMPLICIT NONE

! variables passed to subroutine
       integer cell                    
       integer maxwires                    
       integer Nint_in
       integer Next(3)
       real*8  Pin(maxwires,maxwires,3)
       real*8  Tvin(maxwires,maxwires,3)
       real*8  Tiin(maxwires,maxwires,3)
       integer  BC(maxwires)
       
       real*8  SC(maxwires,3)
       
       real*8  Yf(maxwires,maxwires,3)
       real*8  Vf(maxwires,3)
       
       real*8  Vw(maxwires,3)
       real*8  Iw(maxwires,3)
       
! local_variables
              
       integer nwires,nwires2
       integer Nmax,Nint
       
       real*8 PBC(maxwires,maxwires)
       real*8 S(maxwires,maxwires,3)
       real*8 P(maxwires,maxwires,3)
       real*8 Q(maxwires,maxwires,3)
       real*8 PI(maxwires,maxwires,3)
       real*8 Pnew(maxwires,maxwires,3)

       real*8 W(maxwires,3)
       real*8 U(maxwires,3)
       
       real*8 Tii(maxwires,maxwires)
       real*8 PvT(maxwires,maxwires)
       
       real*8 Q2(maxwires*3,maxwires*3)
       real*8 PI2(maxwires*3,maxwires*3)
       
       real*8 W2(maxwires*2)
       
       real*8 J2(3*maxwires,3*maxwires)
       real*8 Ji2(3*maxwires,3*maxwires)
       real*8 Y2(3*maxwires,3*maxwires)
       
       real*8 A2(3*maxwires,3*maxwires)
       real*8 Ai2(3*maxwires,3*maxwires)
       real*8 X2(3*maxwires)
       
       real*8 V2(3*maxwires)
       
       real*8 Vint(maxwires)
       real*8 V(maxwires,3)
       real*8 IW2(maxwires,3)
       
       real*8 SCZ(maxwires,maxwires,3)
       real*8 SCT(maxwires,maxwires,3)
       
       real*8 TM1(maxwires,maxwires)
       real*8 TM2(maxwires,maxwires)
       
       real*8 TM12(3*maxwires,3*maxwires)
       real*8 TM22(3*maxwires,3*maxwires)
       
       real*8 JiYQ2(3*maxwires,3*maxwires)
       
       real*8 JiYW2(3*maxwires)
       real*8 I2(3*maxwires)
       
       real*8 TV12(3*maxwires)
       real*8 TV22(3*maxwires)

       integer i,j,k,face,nw
       integer row,col
       
       integer face_min,face_max,face_internal
       
       integer Nbc,n_count,Ne            
       
       integer row_offset  

! function_types
      
! START      
       
! reset matrices and vectors

       face_min=1
       face_max=2
       face_internal=3
                   
! Overall dimension of matrix system

       Nmax=Next(1)+Next(2)+Next(3)

! Calucalte PBC

       PBC(:,:)=0d0
       
       Nint=0
       Nbc=0
       do i=1,Nint_in
         if (BC(i).ne.-1) then
           Nint=Nint+1
	 else
	   Nbc=1
	 end if
       end do
     
       Nint=Nint+Nbc
       
       N_count=0 ! counts internal connection nodes which are not subject to BCs
       do i=1,Nint_in
         if (BC(i).ne.-1) then
           N_count=N_count+1
	   Pbc(N_count,i)=1d0
	 else
	   Pbc(Nint,i)=1d0
	 end if      
       end do
       
! Calculate the S matrix

       S(:,:,:)=0d0
       if (Nbc.eq.1) then
         S(Nint,1:Next(1),1)=1d0
         S(Nint,1:Next(2),2)=1d0
         S(Nint,1:Next(3),3)=0d0 ! **** no reference return conductor for load impedances
       end if
       
! calculate the revised PI matrix       
       P(:,:,:)=0d0
       call dmatmul(Pbc,Nint,Nint_in,Pin(1,1,1),Nint_in,Next(1),P(1,1,1),maxwires)       
       call dmatmul(Pbc,Nint,Nint_in,Pin(1,1,2),Nint_in,Next(2),P(1,1,2),maxwires)       
       call dmatmul(Pbc,Nint,Nint_in,Pin(1,1,3),Nint_in,Next(3),P(1,1,3),maxwires)       

! calculate Pnew = P-S (equation 5.1.1)

       Pnew(:,:,:)=P(:,:,:)-S(:,:,:)
       
! calculate the elements of matrix Q for each face, Q=[Tv][Pv]T
       do face=1,3
	 
         CALL dtranspose(P(1,1,face),Nint_in,Next(face),   &
	                 PvT,maxwires)			 
	 CALL dmatmul(Tvin(1,1,face),Next(face),Next(face),PvT,Next(face),Nint_in,Q(1,1,face),maxwires)       
	 
       end do	
       
! calculate the elements of matrix PI for each face, PI=[Pnew][Ti]^-1
! the use of the Pnew matrix is the only real difference between this analysis and that for the cell centre
! junction
       do face=1,3
	 
         call dsvd_invert(Tiin(1,1,face),Next(face),Next(face),Tii,maxwires) 
	 
	 CALL dmatmul(Pnew(1,1,face),Nint_in,Next(face),Tii,Next(face),Next(face),PI(1,1,face),maxwires)       

       end do	

! assemble the full Q matrix and PI matrix from sub matrices, equation 3.3.24 and 3.3.22      
       
       PI2(:,:)=0d0	 
       do i=1,Nint
         do j=1,Next(1)
           row=i
	   col=j
	   PI2(row,col)=PI(i,j,1)
         end do
       end do
       do i=1,Nint
         do j=1,Next(2)
           row=i
	   col=j+Next(1)
	   PI2(row,col)=PI(i,j,2)
         end do
       end do
       
       do i=1,Nint
         do j=1,Next(3)
           row=i
	   col=j+Next(1)+Next(2)
	   PI2(row,col)=PI(i,j,3)
         end do
       end do
       
       Q2(:,:)=0d0	 
       do i=1,Next(1)
         do j=1,Nint
           row=i
	   col=j
	   Q2(row,col)=Q(i,j,1)
         end do
       end do
       do i=1,Next(2)
         do j=1,Nint
           row=i+Next(1)
	   col=j
	   Q2(row,col)=Q(i,j,2)
         end do
       end do
       
       do i=1,Next(3)
         do j=1,Nint
           row=i+Next(1)+Next(2)
	   col=j
	   Q2(row,col)=Q(i,j,3)
         end do
       end do
		
! W vector , equation 3.3.23
       W(:,:)=0d0	
			 
       W(1:Next(face_min),face_min)= Vf(1:Next(face_min),face_min) 		 
       W(1:Next(face_max),face_max)= Vf(1:Next(face_max),face_max) 
       W(1:Next(face_internal),face_internal)= Vf(1:Next(face_internal),face_internal)
				    		 
       W2(:)=0d0        
       do i=1,Next(1)
         row=i
	 W2(row)=W(i,1)
       end do
       
       do i=1,Next(2)
         row=i+Next(1)
	 W2(row)=W(i,2)
       end do
       
       do i=1,Next(3)
         row=i+Next(1)+Next(2)
	 W2(row)=W(i,3)
       end do
		
! J and Y matrices	(equations 3.3.24, 3.3.23) 				
       J2(:,:)=0d0
       Y2(:,:)=0d0
       
       SCZ(:,:,:)=0d0
       SCT(:,:,:)=0d0

!! no coupling to field in TLM - may comment out the next 4 lines for testing TLM...      
!       SCZ(1:Next(1),1,1)=SC(1:Next(1),1)*Zxin(1)
!       SCZ(1:Next(2),1,2)=SC(1:Next(2),2)*Zxin(2)       
!       SCT(1,1:Next(1),1)=SC(1:Next(1),1)
!       SCT(1,1:Next(2),2)=SC(1:Next(2),2)
!       
!       call dmatmul(Yf(1,1,1),Next(1),Next(1),SCZ(1,1,1),Next(1),1,TM1,maxwires)       
!       call dmatmul(TM1,Next(1),1,SCT(1,1,1),1,Next(1),TM2,maxwires)       

       TM2(:,:)=0d0
       
       do i=1,next(1)
	 row=i
         do j=1,next(1)	
	   col=j
	   J2(row,col)=TM2(i,j)
	   Y2(row,col)=Yf(i,j,1)
	 end do
	 J2(row,row)=J2(row,row)+1d0
       end do 
		   
       call dmatmul(Yf(1,1,2),Next(2),Next(2),SCZ(1,1,2),Next(2),1,TM1,maxwires)       
       call dmatmul(TM1,Next(2),1,SCT(1,1,2),1,Next(2),TM2,maxwires)       
       
       do i=1,next(2)
	 row=i+next(1)
         do j=1,next(2)	
	   col=j+next(1)
	   J2(row,col)=TM2(i,j)
	   Y2(row,col)=Yf(i,j,2)
	 end do
	 J2(row,row)=J2(row,row)+1d0
       end do 

! **** internal impedance stuff ****
       
       do i=1,next(3)
	 row=i+next(1)+next(2)
         do j=1,next(3)	
	   col=j+next(1)+next(2)
	   J2(row,col)=0d0
	   Y2(row,col)=Yf(i,j,3)
	 end do
	 J2(row,row)=J2(row,row)+1d0
       end do 

! ****  end of internal impedance stuff *****
       
! build the full matrix system (equation 3.3.27)      
       
       Ne=Next(1)+Next(2)+Next(3)
       
       call dsvd_invert(J2,Ne,Ne,Ji2,maxwires*3) 

! calculate LHS vector (X)      
       call dmatvmul(Y2  ,Ne,Ne ,W2  ,Ne,TV12,maxwires*3) 
       call dmatvmul(Ji2,Ne,Ne  ,TV12,Ne,JiYW2,maxwires*3) 
       call dmatvmul(PI2,Nint,Ne,JiYW2,Ne,X2  ,maxwires*3) 

! calculate Matrix [A]

       call dmatmul(Ji2,Ne,Ne,Y2 ,Ne,Ne  ,TM22,maxwires*3)       
       call dmatmul(TM22,Ne,Ne,Q2,Ne,Ne, JiYQ2 ,maxwires*3)       
       call dmatmul(PI2 ,Nint,Ne,JiYQ2,Ne,Ne  ,A2,maxwires*3)       
       
! Note: reduce the dimension of the system if we have a grounded node       
       call dsvd_invert(A2,Nint-Nbc,Nint-Nbc,Ai2,maxwires*3) 
       
       V2(:)=0d0
       call dmatvmul(Ai2,Nint-Nbc,Nint-Nbc,X2,Nint-Nbc,V2,maxwires*3) 
      
       if (Nbc.eq.1) V2(Nint)=0d0
       
       Vint(:)=0d0
       Vint(1:Nint)=V2(1:Nint)		

! solve for individual conductor voltages (equation 3.3.18)	
       call dmatvmul(Q(1,1,1),Next(1),Nint,Vint,Nint,Vw(1,1),maxwires) 
       call dmatvmul(Q(1,1,2),Next(2),Nint,Vint,Nint,Vw(1,2),maxwires)
       call dmatvmul(Q(1,1,3),Next(3),Nint,Vint,Nint,Vw(1,3),maxwires)
       
! solve for individual conductor currents (equation 3.3.25)
       call dmatvmul(JiYQ2,Ne,Ne  ,V2,Ne,TV22,maxwires*3)

       do row=1,Ne
         I2(row)=JiYW2(row)-TV22(row)
       end do
        
! allocate currents to their correct incident directions        
       row_offset=0
       
       do face=1,3
         do i=1,Next(face)
           row=i+row_offset
	   Iw(i,face)=I2(row)
         end do
	 row_offset=row_offset+Next(face)
      end do
       
      RETURN
		  

     END SUBROUTINE face_junction
       
