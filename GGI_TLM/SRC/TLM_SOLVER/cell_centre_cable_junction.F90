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
! SUBROUTINE cell_centre_cable_junction
! SUBROUTINE junction
!
! NAME
!     cell_centre_cable_junction
!
! DESCRIPTION
!
!     Cable update at a cell centre - at the moment this is an interface to the Fieldsolve subroutine junction()
!     Do the connection process on the inductive stub here for convenience
!
! COMMENTS
!     
!
! HISTORY
!
!     started 24/09/2012 CJS
!             17/12/2012 CJS	implement junction loads
!
!
SUBROUTINE cell_centre_cable_junction(cell_centre_junction, Uxin,Uyin,Uzin ,Yxin,Yyin,Yzin , Ix,Iy,Iz)

USE Cables
USE TLM_General
USE TLM_Excitation
USE cell_parameters		        
USE filter_types		        
USE filter_functions		        

IMPLICIT NONE

  integer	:: cell_centre_junction
  
  real*8	:: Uxin,Uyin,Uzin
  real*8	:: Yxin,Yyin,Yzin
  
  real*8	:: Ix,Iy,Iz

! local variables

  integer maxwires,nw		      
  integer N_int
  integer N_ext(7)  
  
  real*8	:: Zxin,Zyin,Zzin
       
  real*8,allocatable	::  P(:,:,:)
  real*8,allocatable	::  Tv(:,:,:)
  real*8,allocatable	::  Ti(:,:,:)
  
  real*8,allocatable	:: NOT_SC(:,:)
  
  real*8,allocatable	:: Yf(:,:,:)
  real*8,allocatable	:: Zlink(:,:,:)
  real*8,allocatable	:: ZLstub(:,:,:)
  real*8,allocatable	:: Vf(:,:)
  
  real*8,allocatable	:: Vw(:,:)
  real*8,allocatable	:: Iw(:,:)
  
  real*8,allocatable	:: Temp_vector(:)
    
  real*8 Vtx
  real*8 Vty
  real*8 Vtz

  integer	:: face
  integer	:: segment
  integer	:: row,col
  integer	:: excitation_function
  
  integer	:: filter_data
  integer	:: filter
  
  real*8	:: Itemp
  
! START
  
  CALL write_line('CALLED: cell_centre_cable_junction',0,timestepping_output_to_screen_flag)
  
  Zxin=1d0/Yxin
  Zyin=1d0/Yyin
  Zzin=1d0/Yzin
  
  do face=1,7
    N_ext(face)=cell_centre_junction_list(cell_centre_junction)%n_external_conductors(face)
  end do ! next face
  
  N_int=cell_centre_junction_list(cell_centre_junction)%n_internal_connection_nodes
  
  maxwires=max(N_ext(1),N_ext(2),N_ext(3),N_ext(4),N_ext(5),N_ext(6),N_ext(7),N_int)
  
! allocate arrays
  ALLOCATE( P(maxwires,maxwires,7) )
  ALLOCATE( Tv(maxwires,maxwires,7) )
  ALLOCATE( Ti(maxwires,maxwires,7) )
  ALLOCATE( NOT_SC(maxwires,7) )
  ALLOCATE( Yf(maxwires,maxwires,7) )
  ALLOCATE( Zlink(maxwires,maxwires,7) )
  ALLOCATE( ZLstub(maxwires,maxwires,7) )
  ALLOCATE( Vf(maxwires,7) )
  ALLOCATE( Vw(maxwires,7) )
  ALLOCATE( Iw(maxwires,7) )
  ALLOCATE( Temp_vector(maxwires) )
  
! reset arrays

  P(:,:,:)=0d0
  Tv(:,:,:)=0d0
  NOT_SC(:,:)=0d0
  Yf(:,:,:)=0d0
  Zlink(:,:,:)=0d0
  ZLstub(:,:,:)=0d0
  Vf(:,:)=0d0
  Vw(:,:)=0d0
  Iw(:,:)=0d0
  Temp_vector(:)=0d0
  
! START process for connecting bundle segments  

  do face=1,6
  
    segment=cell_centre_junction_list(cell_centre_junction)%segment_list(face)
    
    if (segment.ne.0) then
! there is a cable on this face so fill the appropriate data from the bundle segment data structure

      if (N_ext(face).NE.bundle_segment_list(segment)%n_conductors) then
        write(*,*)'Error in cell_centre_cable_junction'
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
          Yf(row,col,face)=bundle_segment_list(segment)%Yf(row,col)
          Zlink(row,col,face)=bundle_segment_list(segment)%Zlink(row,col)
          ZLstub(row,col,face)=bundle_segment_list(segment)%ZLstub(row,col)
	
	end do ! next col
	
      end do ! next row
      
! Copy P matrix with required type conversion from integer to real*8      
      do col=1,N_ext(face)
        	
        do row=1,N_int
	
	  P(row,col,face) =dble(cell_centre_junction_list(cell_centre_junction)%P_matrix_list(face)%P(row,col))
		
	end do ! next col
	
      end do ! next row
      
! Work out the cable voltage vector, Vf including cable voltage source where required
  
      do row=1,N_ext(face)
    
        Vf(row,face)=2d0*bundle_segment_list(segment)%Vlink(row)+2d0*bundle_segment_list(segment)%VLstub(row)
	
	excitation_function=bundle_segment_list(segment)%excitation_function(row)
	
	if (excitation_function.ne.0) then
! add source voltage

          Vf(row,face)=Vf(row,face)+excitation_functions(excitation_function)%value(timestep)	&
	                           *bundle_segment_list(segment)%direction_sign_list(row)
	  
        end if ! source voltage present
	
      end do ! next row
      
! Add voltage sources related to the impedance filter functions  

! Evaluate Impedance filters
      if (bundle_segment_list(segment)%n_filter_data.NE.0) then

! loop over filters allocating memory for the particular filter order
        filter_data=0
        do row=1,N_ext(face)
          do col=1,N_ext(face)
      
            filter=bundle_segment_list(segment)%filter_number(row,col)
	    if (filter.ne.0) then ! there is a filter on this matrix element
	    
	      filter_data=filter_data+1
              CALL timeshift_Zfilter(bundle_segment_list(segment)%Zfilter_data(filter_data))
              CALL evaluate_Zfilter (bundle_segment_list(segment)%Zfilter(filter),	 &
	        		     bundle_segment_list(segment)%Zfilter_data(filter_data),	 &
			  	     0d0)
				    		   
	      Vf(row,face)=Vf(row,face)+bundle_segment_list(segment)%Zfilter_data(filter_data)%f
	      
	    end if
	    
          end do
        end do
 
      end if  ! n_filter_data.NE.0
      
    else ! no cable segment on this face

      N_ext(face)=0

    end if  ! segment=0 so no cable on this face
    
  end do ! next face
  
! END process for connecting bundle segments  
  
! START of process for junction impedances.
  face=7
  
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
      Yf(row,col,face)=cell_centre_junction_list(cell_centre_junction)%Yf(row,col)

    end do ! next col

  end do ! next row
      
! Copy P matrix with required type conversion from integer to real*8      
  do col=1,N_ext(face)
  	    
    do row=1,N_int

      P(row,col,face) =dble(cell_centre_junction_list(cell_centre_junction)%P_matrix_list(face)%P(row,col))
  	    
    end do ! next col

  end do ! next row
   	 
! Voltage sources related to the impedance filter functions  

! Evaluate Impedance filters
  do row=1,cell_centre_junction_list(cell_centre_junction)%n_internal_impedance_filters

    filter=row
    CALL timeshift_Zfilter(cell_centre_junction_list(cell_centre_junction)%Zfilter_data(filter))
    CALL evaluate_Zfilter (cell_centre_junction_list(cell_centre_junction)%Zfilter(filter),		&
    			   cell_centre_junction_list(cell_centre_junction)%Zfilter_data(filter), &
    			   0d0)
    Vf(row,face)=Vf(row,face)+cell_centre_junction_list(cell_centre_junction)%Zfilter_data(filter)%f

  end do
  
! END of internal junction impedance process
  
  Vtx=0d0
  Vty=0d0
  Vtz=0d0
  Ix=0d0
  Iy=0d0
  Iz=0d0
  
! calculate the cable voltages and currents at the junction
  CALL junction(cell_centre_junction,maxwires,N_int,N_ext,P,Tv,Ti,NOT_SC,Yf,Vf,Vw,Iw,  &	     
                Zxin,Uxin,Zyin,Uyin,Zzin,Uzin,     &
		Vtx,Ix,Vty,Iy,Vtz,Iz		    )

! calculate the voltage across the lines and stubs in all directions
	   
  do face=1,6
	   
    segment=cell_centre_junction_list(cell_centre_junction)%segment_list(face)
    
    if (segment.ne.0) then
    
! there is a cable on this face so calculate the scattered voltages (Vr)=(Vi)-[Z]*(Iw)

      nw=N_ext(face)

      call dmatvmul(Zlink(1,1,face),nw,nw,Iw(1,face),nw,Temp_vector,maxwires)
!      call dmatvmul(bundle_segment_list(segment)%Zlink,nw,nw,Iw(1,face),nw,Temp_vector,maxwires)
      do row=1,nw
        bundle_segment_list(segment)%Vlink(row)=bundle_segment_list(segment)%Vlink(row)-Temp_vector(row)
      end do
      
      call dmatvmul(ZLstub(1,1,face),nw,nw,Iw(1,face),nw,Temp_vector,maxwires)
!      call dmatvmul(bundle_segment_list(segment)%ZLstub,nw,nw,Iw(1,face),nw,Temp_vector,maxwires)
      do row=1,nw
        bundle_segment_list(segment)%VLstub(row)=bundle_segment_list(segment)%VLstub(row)-Temp_vector(row)
      end do
    
! Do the connection process on the inductive stub here for convenience
      do row=1,nw
        bundle_segment_list(segment)%VLstub(row)=-bundle_segment_list(segment)%VLstub(row)
      end do
      
! save the cable current vector
      do row=1,nw
        bundle_segment_list(segment)%Iw_centre(row)=Iw(row,face)	&
	                                           *bundle_segment_list(segment)%direction_sign_list(row)
      end do

      if (bundle_segment_list(segment)%n_filter_data.NE.0) then

! loop over filters allocating memory for the particular filter order
        filter_data=0
        do row=1,N_ext(face)
   	  do col=1,N_ext(face)
    
   	    filter=bundle_segment_list(segment)%filter_number(row,col)
            if (filter.ne.0) then ! there is a filter on this matrix element
       
	      filter_data=filter_data+1
   	      Itemp=-Iw(col,face)
   	      CALL evaluate_Zfilter(bundle_segment_list(segment)%Zfilter(filter),        &
        	  		  bundle_segment_list(segment)%Zfilter_data(filter_data),      &
        	  		  Itemp)
            
            end if
          
   	  end do
        end do
 
      end if  ! n_filter_data.NE.0
	       	       
    end if ! there is a bundle segment on this face
	          
  end do ! next face
  
! Internal junction impedance process  
  face=7
 
! Update Internal Impedance filters
  do row=1,cell_centre_junction_list(cell_centre_junction)%n_internal_impedance_filters

    filter=row
    Itemp=-Iw(row,face)
    CALL evaluate_Zfilter (cell_centre_junction_list(cell_centre_junction)%Zfilter(filter),		&
    			   cell_centre_junction_list(cell_centre_junction)%Zfilter_data(filter), &
    			   Itemp)

  end do
    
! deallocate arrays
  DEALLOCATE( P )
  DEALLOCATE( Tv )
  DEALLOCATE( Ti )
  DEALLOCATE( NOT_SC )
  DEALLOCATE( Yf )
  DEALLOCATE( Zlink )
  DEALLOCATE( ZLstub )
  DEALLOCATE( Vf )
  DEALLOCATE( Vw )
  DEALLOCATE( Iw )
  DEALLOCATE( Temp_vector )
 
  CALL write_line('FINISHED: cell_centre_cable_junction',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE cell_centre_cable_junction
!
! Name junction
!     
!
! Description
!     TLM MTL general TLM node junction formulation
!     This subroutine has been taken directly from Fieldsolve and needs to be re-written
!
! Comments:
!     Equation numbers refer to the GGI_TLM_cable_model_theory document
!
! History
!
!     started 17/03/09 CJS
!     Include shielded cables into the formulation 9/11/09
!     junction internal impedances included 8/03/2011 CJS
!     made TLM specific for testing, original saved as: junction.F90_saved_8_3_2011
!
       SUBROUTINE junction(cell,maxwires,Nint_in,Next,P,Tv,Ti,SC,Yf,Vf,Vw,Iw,  &		
                           Zxin,Uxin,Zyin,Uyin,Zzin,Uzin,     &
		           Vtx,Ix,Vty,Iy,Vtz,Iz                )
USE cell_parameters		        

IMPLICIT NONE

! variables passed to subroutine
       integer cell                    
       integer maxwires                    
       integer Nint_in
       integer Next(7)
       real*8 P(maxwires,maxwires,7)
       real*8 Tv(maxwires,maxwires,7)
       real*8 Ti(maxwires,maxwires,7)
       
       real*8 SC(maxwires,7)
       
       real*8 Yf(maxwires,maxwires,7)
       real*8 Vf(maxwires,7)
       
       real*8 Vw(maxwires,7)
       real*8 Iw(maxwires,7)

       real*8 Uxin
       real*8 Uyin
       real*8 Uzin
       real*8 Zxin
       real*8 Zyin
       real*8 Zzin
       
       real*8 Ix,Iy,Iz
         
       real*8 Vtx
       real*8 Vty
       real*8 Vtz
       
! local_variables
              
       integer nwires,nwires2
       integer Nmax,Nint
       
       real*8 Tii(maxwires,maxwires)
       real*8 PvT(maxwires,maxwires)
       
       real*8 Q(maxwires,maxwires,7)
       real*8 PI(maxwires,maxwires,7)

       real*8 W(maxwires,7)
       real*8 U(maxwires,7)
       
       real*8 Q2(maxwires*7,maxwires*7)
       real*8 PI2(maxwires*7,maxwires*7)
       
       real*8 W2(maxwires*7)
       
       real*8 J2(7*maxwires,7*maxwires)
       real*8 Ji2(7*maxwires,7*maxwires)
       real*8 Y2(7*maxwires,7*maxwires)
       
       real*8 A2(7*maxwires,7*maxwires)
       real*8 Ai2(7*maxwires,7*maxwires)
       real*8 X2(7*maxwires)
       
       real*8 V2(7*maxwires)
       
       real*8 Vint(maxwires)
       real*8 V(maxwires,7)
       real*8 IW2(maxwires,7)
       
       real*8 SCZ(maxwires,maxwires,7)
       real*8 SCT(maxwires,maxwires,7)
       
       real*8 TM1(maxwires,maxwires)
       real*8 TM2(maxwires,maxwires)
       
       real*8 TM12(7*maxwires,7*maxwires)
       real*8 TM22(7*maxwires,7*maxwires)
       real*8 JiYQ2(7*maxwires,7*maxwires)
       
       real*8 JiYW2(7*maxwires)
       real*8 I2(7*maxwires)
       
       real*8 TV1(maxwires)
       real*8 TV2(maxwires)
       real*8 TV3(maxwires)
       real*8 TV4(maxwires)
       
       real*8 TV12(maxwires*7)
       real*8 TV22(maxwires*7)
       real*8 TV32(maxwires*7)
       real*8 TV42(maxwires*7)

       integer i,j,k,face,nw
       integer row,col
       integer row_offset,col_offset
       integer row_offset2,col_offset2
              
       integer Nbc,n_count,Ne              

       integer direction,face_min,face_max
       integer faceloop,face_list(7)
       integer face_internal
       parameter (face_internal=7)

! function_types
      
! START      
       
! reset matrices and vectors
       
       face_list(1)=face_xmin
       face_list(2)=face_xmax
       face_list(3)=face_ymin
       face_list(4)=face_ymax
       face_list(5)=face_zmin
       face_list(6)=face_zmax
       face_list(7)=face_internal
             
! Overall diension of matrix system

       Nmax=Next(1)+Next(2)+Next(3)+Next(4)+Next(5)+Next(6)+Next(7)
     
       Nint=Nint_in
       
       PI(:,:,:)=P(:,:,:)
       
! calculate the elements of matrix Q for each face, Q=[Tv][Pv]T (equation 3.3.18)
       do faceloop=1,7
         face=face_list(faceloop)
	 
         CALL dtranspose(P(1,1,face),Nint_in,Next(face),   &
	                 PvT,maxwires)			 
	 CALL dmatmul(Tv(1,1,face),Next(face),Next(face),PvT,Next(face),Nint_in,Q(1,1,face),maxwires)       
	 
       end do	
       
! calculate the elements of matrix PI for each face, PI=[P][Ti]^-1 (equation 3.3.24)
       do faceloop=1,7
         face=face_list(faceloop)
	 
         call dsvd_invert(Ti(1,1,face),Next(face),Next(face),Tii,maxwires) 
	 
	 CALL dmatmul(P(1,1,face),Nint_in,Next(face),Tii,Next(face),Next(face),PI(1,1,face),maxwires)       

       end do	

! assemble the full Q matrix and PI matrix from sub matrices, equation 3.3.24 and 3.3.22      
       PI2(:,:)=0d0	 
       
       col_offset=0
       
       do faceloop=1,7
         face=face_list(faceloop)
         do i=1,Nint
           do j=1,Next(face)
             row=i
	     col=j+col_offset
	     PI2(row,col)=PI(i,j,face)
           end do
         end do
	 col_offset=col_offset+Next(face)
       end do
       
       Q2(:,:)=0d0	 
       
       row_offset=0
       
       do faceloop=1,7
         face=face_list(faceloop)
         do i=1,Next(face)
           do j=1,Nint
             row=i+row_offset
	     col=j
	     Q2(row,col)=Q(i,j,face)
           end do
         end do
	 row_offset=row_offset+Next(face)
      end do
       		
! W vector, equation 3.3.23
       W(:,:)=0d0	
		
       W(1:Next(face_xmin),face_xmin)= Vf(1:Next(face_xmin),face_xmin)   &
        			  -SC(1:Next(face_xmin),face_xmin)*Uxin/2d0   
       W(1:Next(face_xmax),face_xmax)= Vf(1:Next(face_xmax),face_xmax)   &
        			  +SC(1:Next(face_xmax),face_xmax)*Uxin/2d0
       W(1:Next(face_ymin),face_ymin)= Vf(1:Next(face_ymin),face_ymin)   &
        			  -SC(1:Next(face_ymin),face_ymin)*Uyin/2d0   
       W(1:Next(face_ymax),face_ymax)= Vf(1:Next(face_ymax),face_ymax)   &
        			  +SC(1:Next(face_ymax),face_ymax)*Uyin/2d0
       W(1:Next(face_zmin),face_zmin)= Vf(1:Next(face_zmin),face_zmin)   &
        			  -SC(1:Next(face_zmin),face_zmin)*Uzin/2d0   
       W(1:Next(face_zmax),face_zmax)= Vf(1:Next(face_zmax),face_zmax)   &
        			  +SC(1:Next(face_zmax),face_zmax)*Uzin/2d0

! **** note SC(:,face_internal) should be zero. no coupling to external field of internal impedances   				  
       W(1:Next(face_internal),face_internal)= Vf(1:Next(face_internal),face_internal)
       
       row_offset=0
       col_offset=0
       
       W2(:)=0d0        
       do faceloop=1,7
         face=face_list(faceloop)
         do i=1,Next(face)
           row=i+row_offset
	   col=j
	   W2(row)=W(i,face)
         end do
	 row_offset=row_offset+Next(face)
       end do
       		
       		 		
! J and Y matrices (equations 3.3.24, 3.3.23) 			
       J2(:,:)=0d0
       Y2(:,:)=0d0
       
       SCZ(:,:,:)=0d0
       SCT(:,:,:)=0d0
       
       SCZ(1:Next(face_xmin),1,face_xmin)=SC(1:Next(face_xmin),face_xmin)*Zxin/4d0
       SCZ(1:Next(face_xmax),1,face_xmax)=SC(1:Next(face_xmax),face_xmax)*Zxin/4d0
       SCZ(1:Next(face_ymin),1,face_ymin)=SC(1:Next(face_ymin),face_ymin)*Zyin/4d0
       SCZ(1:Next(face_ymax),1,face_ymax)=SC(1:Next(face_ymax),face_ymax)*Zyin/4d0
       SCZ(1:Next(face_zmin),1,face_zmin)=SC(1:Next(face_zmin),face_zmin)*Zzin/4d0
       SCZ(1:Next(face_zmax),1,face_zmax)=SC(1:Next(face_zmax),face_zmax)*Zzin/4d0
       
       SCZ(1:Next(face_internal),1,face_internal)=0d0
       
       do faceloop=1,7
         face=face_list(faceloop)
         SCT(1,1:Next(face),face)=SC(1:Next(face),face)
       end do
       
       do direction=1,3
         if (direction.eq.1) then
	   face_min=face_xmin
	   face_max=face_xmax
	   row_offset=0
	   col_offset=0
	   row_offset2=Next(face_xmin)
	   col_offset2=Next(face_xmin)
	 else if (direction.eq.2) then
	   face_min=face_ymin
	   face_max=face_ymax
	   row_offset=Next(face_xmin)+Next(face_xmax)
	   col_offset=Next(face_xmin)+Next(face_xmax)
	   row_offset2=row_offset+Next(face_ymin)
	   col_offset2=col_offset+Next(face_ymin)
	 else
	   face_min=face_zmin
	   face_max=face_zmax
	   row_offset=Next(face_xmin)+Next(face_xmax)+Next(face_ymin)+Next(face_ymax)
	   col_offset=Next(face_xmin)+Next(face_xmax)+Next(face_ymin)+Next(face_ymax)
	   row_offset2=row_offset+Next(face_zmin)
	   col_offset2=col_offset+Next(face_zmin)
	 end if
       
         call dmatmul(Yf(1,1,face_min),Next(face_min),Next(face_min),SCZ(1,1,face_min),Next(face_min),1,TM1,maxwires)       
         call dmatmul(TM1,Next(face_min),1,SCT(1,1,face_min),1,Next(face_min),TM2,maxwires)       
       
         do i=1,Next(face_min)
	   row=i+row_offset
           do j=1,Next(face_min)	
	     col=j+col_offset
	     J2(row,col)=TM2(i,j)
	     Y2(row,col)=Yf(i,j,face_min)
	   end do
	   J2(row,row)=J2(row,row)+1d0
         end do 
       
         call dmatmul(Yf(1,1,face_min),Next(face_min),Next(face_min),SCZ(1,1,face_min),Next(face_min),1,TM1,maxwires)       
         call dmatmul(TM1,Next(face_min),1,SCT(1,1,face_max),1,Next(face_max),TM2,maxwires)       
       
         do i=1,Next(face_min)
	   row=i+row_offset
           do j=1,Next(face_max)	
	     col=j+col_offset2
	     J2(row,col)=-TM2(i,j)
	   end do
         end do 
       

         call dmatmul(Yf(1,1,face_max),Next(face_max),Next(face_max),SCZ(1,1,face_max),Next(face_max),1,TM1,maxwires)	    
         call dmatmul(TM1,Next(face_max),1,SCT(1,1,face_min),1,Next(face_min),TM2,maxwires)	  
       
         do i=1,Next(face_max)
	   row=i+row_offset2
           do j=1,Next(face_min)      
	     col=j+col_offset
	     J2(row,col)=-TM2(i,j)	 
	   end do     
         end do 
	         
         call dmatmul(Yf(1,1,face_max),Next(face_max),Next(face_max),SCZ(1,1,face_max),Next(face_max),1,TM1,maxwires)	    
         call dmatmul(TM1,Next(face_max),1,SCT(1,1,face_max),1,Next(face_max),TM2,maxwires)	  
       
         do i=1,Next(face_max)
	   row=i+row_offset2
           do j=1,Next(face_max)      
	     col=j+col_offset2
  	     J2(row,col)=TM2(i,j)
	     Y2(row,col)=Yf(i,j,face_max)
	   end do
	   J2(row,row)=J2(row,row)+1d0
         end do 
              
       end do ! next direction

! **** internal impedance stuff ****

       if (Next(face_internal).ne.0) then

	 row_offset=Next(face_xmin)+Next(face_xmax)+Next(face_ymin)+Next(face_ymax)+Next(face_zmin)+Next(face_zmax)
	 col_offset=Next(face_xmin)+Next(face_xmax)+Next(face_ymin)+Next(face_ymax)+Next(face_zmin)+Next(face_zmax)

         call dmatmul(Yf(1,1,face_internal),Next(face_internal),	&
	              Next(face_internal),SCZ(1,1,face_internal),	&
		      Next(face_internal),1,TM1,maxwires)       
         call dmatmul(TM1,Next(face_internal),1,SCT(1,1,face_internal),1,Next(face_internal),TM2,maxwires)       
       
         do i=1,Next(face_internal)
	   row=i+row_offset
           do j=1,Next(face_internal)	
	     col=j+col_offset
	     J2(row,col)=TM2(i,j)
	     Y2(row,col)=Yf(i,j,face_internal)
	   end do
	   J2(row,row)=J2(row,row)+1d0
         end do 
       	  
       end if

! ****  end of internal impedance stuff *****

! build the full matrix system (equation 3.3.27)      
       Ne=Nmax
       
       call dsvd_invert(J2,Ne,Ne,Ji2,maxwires*7) 

! calculate LHS vector (X)      
       call dmatvmul(Y2  ,Ne,Ne ,W2  ,Ne,TV12,maxwires*7) 
       call dmatvmul(Ji2,Ne,Ne  ,TV12,Ne,JiYW2,maxwires*7) 
       call dmatvmul(PI2,Nint,Ne,JiYW2,Ne,X2  ,maxwires*7) 

! calculate Matrix [A]

       call dmatmul(Ji2,Ne,Ne,Y2 ,Ne,Ne  ,TM22,maxwires*7)       
       call dmatmul(TM22,Ne,Ne,Q2,Ne,Ne, JiYQ2 ,maxwires*7)       
       call dmatmul(PI2 ,Nint,Ne,JiYQ2,Ne,Ne  ,A2,maxwires*7)       

! solve for the internal connection node voltages       
       call dsvd_invert(A2,Nint,Nint,Ai2,maxwires*7) 
       
       V2(:)=0d0
       call dmatvmul(Ai2,Nint,Nint,X2,Nint,V2,maxwires*7) 
       
       Vint(:)=0d0
       Vint(1:Nint)=V2(1:Nint)		

! solve for individual conductor voltages (equation 3.3.18)		
       do faceloop=1,7
         face=face_list(faceloop)
         call dmatvmul(Q(1,1,face),Next(face),Nint,Vint,Nint,Vw(1,face),maxwires) 
       end do
       
! solve for individual conductor currents (equation 3.3.25)
       call dmatvmul(JiYQ2,Ne,Ne  ,V2,Ne,TV22,maxwires*7)
        
       do row=1,Ne
         I2(row)=JiYW2(row)-TV22(row)
       end do

! allocate currents to their correct incident directions        
       row_offset=0
       
       do faceloop=1,7
         face=face_list(faceloop)
         do i=1,Next(face)
           row=i+row_offset
	   Iw(i,face)=I2(row)
         end do
	 row_offset=row_offset+Next(face)
      end do

! calculate the transformer current ,Ix (equation 3.3.8)      		
       Ix=0d0
       do i=1,Next(face_xmin)
         Ix=Ix-Iw(i,face_xmin)*SC(i,face_xmin)/2.0d0
       end do
       do i=1,Next(face_xmax)
         Ix=Ix+Iw(i,face_xmax)*SC(i,face_xmax)/2.0d0
       end do

! transformer voltages on the field side in the x, y and z directions (equation 3.3.7)               
       Vtx=Uxin-Zxin*Ix
       		
       Iy=0d0
       do i=1,Next(face_ymin)
         Iy=Iy-Iw(i,face_ymin)*SC(i,face_ymin)/2.0d0
       end do
       do i=1,Next(face_ymax)
         Iy=Iy+Iw(i,face_ymax)*SC(i,face_ymax)/2.0d0
       end do
              
       Vty=Uyin-Zyin*Iy
		       		
       Iz=0d0
       do i=1,Next(face_zmin)
         Iz=Iz-Iw(i,face_zmin)*SC(i,face_zmin)/2.0d0
       end do
       do i=1,Next(face_zmax)
         Iz=Iz+Iw(i,face_zmax)*SC(i,face_zmax)/2.0d0
       end do
              
       Vtz=Uzin-Zzin*Iz
         
       return      

       END SUBROUTINE junction
      
      
