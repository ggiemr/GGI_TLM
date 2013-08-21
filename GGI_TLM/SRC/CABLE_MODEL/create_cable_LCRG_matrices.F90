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
!SUBROUTINE create_cable_LCRG_matrices
!
! NAME
!     SUBROUTINE create_cable_LCRG_matrices
!
! DESCRIPTION
!     calculate the internal LC matrices of each cable i.e. the L and C matrices
!     of the portion of the cable inside the shield
!     Also calculate the Z domain filters for frequency dependent cable parameters
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/09/2012 CJS
!     include coaxial cable 19/11/2012 CJS
!     include coaxial cable transfer impedance 29/11/2012 CJS
!     include frequency domain filters 14/12/2012 CJS
!     frequency warping included 12/02/2013 CJS
!     include ribbon cable 24/04/2013 CJS
!
SUBROUTINE create_cable_LCRG_matrices()

USE TLM_general
USE Cables
USE pul_data
USE constants
USE filter_types
USE filter_functions

IMPLICIT NONE

! local variables

  integer	:: cable_geometry_number
  integer	:: n_conductors,n_external,n_shielded
  
  real*8	:: r_inner,r_shield,epsr
  real*8	:: r_dielectric,separation,ribbon_width,ribbon_centre_width
  real*8	:: Lw,Cw
  real*8	:: Lt,Rt
  
  integer	:: n_params
  integer	:: n_filters
  integer	:: filter
  
  TYPE(Zfilter) :: Zfilter1
  TYPE(Zfilter) :: Zfilter2
  real*8	:: Z_f

  integer	:: row,col,i

! START

  CALL write_line('CALLED: create_cable_LCRG_matrices',0,output_to_screen_flag)

! calculate Z domain filters  

  do cable_geometry_number=1,n_cable_geometries
  
    n_filters=cable_geometry_list(cable_geometry_number)%n_filters
    
    do filter=1,n_filters

! scale filters from per metre value to per cell_segment value    
       
        cable_geometry_list(cable_geometry_number)%Sfilter(filter)%a%coeff(:)=	&
	      cable_geometry_list(cable_geometry_number)%Sfilter(filter)%a%coeff(:)*dl/2d0
	
! Calculate Z domain susceptibility admittance filter coefficients by bilinear transformation
!        Zfilter1=s_to_z(cable_geometry_list(cable_geometry_number)%Sfilter(filter),dt) 
        Zfilter1=s_to_z_warp(cable_geometry_list(cable_geometry_number)%Sfilter(filter),dt,bicubic_warp_flag,frequency_scale) 
          
        write(*,*)'Zfilter1'
        CALL write_Zfilter(Zfilter1,0)
    
! Calculate fast impedance and slow impedance filter

        CALL Z_fast_slow_docomposition( Zfilter1 ,Z_f  , Zfilter2 )
          
        write(*,*)'Zs_eps_f',Z_f
        write(*,*)'Zfilter2'
        CALL write_Zfilter(Zfilter2,0)

        cable_geometry_list(cable_geometry_number)%Zfilter(filter)=Zfilter2
        cable_geometry_list(cable_geometry_number)%Z_f(filter)=Z_f
    
    end do
    
  end do ! next cable geometry

! Calculate matrices associated with the cable geometry
  do cable_geometry_number=1,n_cable_geometries
  
    if ((cable_geometry_list(cable_geometry_number)%cable_geometry_type.EQ.cable_geometry_type_cylindrical).OR.	      &
        (cable_geometry_list(cable_geometry_number)%cable_geometry_type.EQ.cable_geometry_type_FD_cylindrical) ) then
    
      n_params=3
      if (cable_geometry_list(cable_geometry_number)%n_parameters.ne.n_params) GOTO 9010     
    
      n_conductors=1
      n_external=1
      n_shielded=0
      
      cable_geometry_list(cable_geometry_number)%n_conductors         =n_conductors      
      cable_geometry_list(cable_geometry_number)%n_shielded_conductors=n_shielded
      cable_geometry_list(cable_geometry_number)%n_external_conductors=n_external
      
      ALLOCATE( cable_geometry_list(cable_geometry_number)%external_conductor_xc(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%external_conductor_yc(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%external_conductor_radius(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%external_dielectric_radius(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%external_dielectric_permittivity(1:n_external) )
           
      cable_geometry_list(cable_geometry_number)%external_conductor_xc(1)=0d0
      cable_geometry_list(cable_geometry_number)%external_conductor_yc(1)=0d0
      
      cable_geometry_list(cable_geometry_number)%external_conductor_radius(1)=	&
      				cable_geometry_list(cable_geometry_number)%parameters(1)
				
      cable_geometry_list(cable_geometry_number)%external_dielectric_radius(1)=	&
      				cable_geometry_list(cable_geometry_number)%parameters(2)
				
      cable_geometry_list(cable_geometry_number)%external_dielectric_permittivity(1)=	&
      				cable_geometry_list(cable_geometry_number)%parameters(3)

      cable_geometry_list(cable_geometry_number)%cable_offset_radius=1.01d0*		&
         cable_geometry_list(cable_geometry_number)%external_dielectric_radius(1) 
 				
      ALLOCATE( cable_geometry_list(cable_geometry_number)%L_internal(1:n_conductors,1:n_conductors) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%C_internal(1:n_conductors,1:n_conductors) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%R_internal(1:n_conductors,1:n_conductors) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%Tv(1:n_conductors,1:n_conductors) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%Ti(1:n_conductors,1:n_conductors) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%SC(1:n_conductors) )

      cable_geometry_list(cable_geometry_number)%L_internal(1:n_conductors,1:n_conductors)=0d0
      cable_geometry_list(cable_geometry_number)%C_internal(1:n_conductors,1:n_conductors)=0d0
      cable_geometry_list(cable_geometry_number)%R_internal(1:n_conductors,1:n_conductors)=0d0
      cable_geometry_list(cable_geometry_number)%Tv(n_conductors,n_conductors)=1
      cable_geometry_list(cable_geometry_number)%Ti(n_conductors,n_conductors)=1
      cable_geometry_list(cable_geometry_number)%SC(n_conductors)=0

      ALLOCATE( cable_geometry_list(cable_geometry_number)%filter_number(1:n_conductors,1:n_conductors) )

      if (cable_geometry_list(cable_geometry_number)%cable_geometry_type.EQ.cable_geometry_type_cylindrical) then

        n_filters=0
        if (cable_geometry_list(cable_geometry_number)%n_filters.ne.n_filters) GOTO 9020
        cable_geometry_list(cable_geometry_number)%filter_number(1:n_conductors,1:n_conductors)=0
      
      else if (cable_geometry_list(cable_geometry_number)%cable_geometry_type.EQ.cable_geometry_type_FD_cylindrical) then
    
! set cable impedance filters
        n_filters=1
        if (cable_geometry_list(cable_geometry_number)%n_filters.ne.n_filters) GOTO 9020
        cable_geometry_list(cable_geometry_number)%filter_number(1,1)=1

      end if
      
    else if ((cable_geometry_list(cable_geometry_number)%cable_geometry_type.EQ.cable_geometry_type_coaxial).OR.   &
             (cable_geometry_list(cable_geometry_number)%cable_geometry_type.EQ.cable_geometry_type_FD_coaxial) ) then
    
      if (cable_geometry_list(cable_geometry_number)%cable_geometry_type.EQ.cable_geometry_type_coaxial) then

        n_params=7
        if (cable_geometry_list(cable_geometry_number)%n_parameters.ne.n_params) GOTO 9010     

      else ! FD_coax cable

        n_params=5
        if (cable_geometry_list(cable_geometry_number)%n_parameters.ne.n_params) GOTO 9010     

      end if
      
      n_conductors=2
      n_external=1
      n_shielded=1
      
      cable_geometry_list(cable_geometry_number)%n_conductors         =n_conductors      
      cable_geometry_list(cable_geometry_number)%n_shielded_conductors=n_shielded
      cable_geometry_list(cable_geometry_number)%n_external_conductors=n_external
      
      ALLOCATE( cable_geometry_list(cable_geometry_number)%external_conductor_xc(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%external_conductor_yc(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%external_conductor_radius(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%external_dielectric_radius(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%external_dielectric_permittivity(1:n_external) )
           
      cable_geometry_list(cable_geometry_number)%external_conductor_xc(1)=0d0
      cable_geometry_list(cable_geometry_number)%external_conductor_yc(1)=0d0
      
      cable_geometry_list(cable_geometry_number)%external_conductor_radius(1)=	&
      				cable_geometry_list(cable_geometry_number)%parameters(1)
				
      cable_geometry_list(cable_geometry_number)%external_dielectric_radius(1)=	&
      				cable_geometry_list(cable_geometry_number)%parameters(2)
				
      cable_geometry_list(cable_geometry_number)%external_dielectric_permittivity(1)=	&
      				cable_geometry_list(cable_geometry_number)%parameters(3)

      cable_geometry_list(cable_geometry_number)%cable_offset_radius=1.01d0*		&
         cable_geometry_list(cable_geometry_number)%external_dielectric_radius(1) 
      
      ALLOCATE( cable_geometry_list(cable_geometry_number)%shielded_conductor_xc(1:n_shielded) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%shielded_conductor_yc(1:n_shielded) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%shielded_conductor_radius(1:n_shielded) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%shielded_dielectric_radius(1:n_shielded) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%shielded_dielectric_permittivity(1:n_shielded) )
           
      cable_geometry_list(cable_geometry_number)%shielded_conductor_xc(1)=0d0
      cable_geometry_list(cable_geometry_number)%shielded_conductor_yc(1)=0d0
      
      cable_geometry_list(cable_geometry_number)%shielded_conductor_radius(1)=	&
      				cable_geometry_list(cable_geometry_number)%parameters(4)
				
      cable_geometry_list(cable_geometry_number)%shielded_dielectric_radius(1)=	&
      				cable_geometry_list(cable_geometry_number)%parameters(1)
				
      cable_geometry_list(cable_geometry_number)%shielded_dielectric_permittivity(1)=	&
      				cable_geometry_list(cable_geometry_number)%parameters(5)
				
      ALLOCATE( cable_geometry_list(cable_geometry_number)%L_internal(1:n_conductors,1:n_conductors) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%C_internal(1:n_conductors,1:n_conductors) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%R_internal(1:n_conductors,1:n_conductors) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%Tv(1:n_conductors,1:n_conductors) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%Ti(1:n_conductors,1:n_conductors) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%SC(1:n_conductors) )

      cable_geometry_list(cable_geometry_number)%L_internal(1:n_conductors,1:n_conductors)=0d0
      cable_geometry_list(cable_geometry_number)%C_internal(1:n_conductors,1:n_conductors)=0d0
      cable_geometry_list(cable_geometry_number)%R_internal(1:n_conductors,1:n_conductors)=0d0
      cable_geometry_list(cable_geometry_number)%Tv(1:n_conductors,1:n_conductors)=0
      cable_geometry_list(cable_geometry_number)%Ti(1:n_conductors,1:n_conductors)=0
      cable_geometry_list(cable_geometry_number)%SC(1:n_conductors)=0
      
      r_shield=cable_geometry_list(cable_geometry_number)%parameters(1)
      r_inner =cable_geometry_list(cable_geometry_number)%parameters(4)
      epsr    =cable_geometry_list(cable_geometry_number)%parameters(5)
    
      Lw=(mu0/(2d0*pi))*log(r_shield/r_inner)
    
      Cw=(2d0*pi*eps0*epsr)/log(r_shield/r_inner)
      
      if (cable_geometry_list(cable_geometry_number)%cable_geometry_type.EQ.cable_geometry_type_coaxial) then
        Lt    =cable_geometry_list(cable_geometry_number)%parameters(6)
        Rt    =cable_geometry_list(cable_geometry_number)%parameters(7)
      else
        Lt    =0d0   ! transfer impedance implemented as a filter function
        Rt    =0d0   
      end if 
      
      cable_geometry_list(cable_geometry_number)%L_internal(1,1)=Lw*dl/2d0
      cable_geometry_list(cable_geometry_number)%L_internal(1,2)=Lt*dl/2d0
      cable_geometry_list(cable_geometry_number)%L_internal(2,1)=Lt*dl/2d0
      
      cable_geometry_list(cable_geometry_number)%C_internal(1,1)=Cw*dl/2d0
      
      cable_geometry_list(cable_geometry_number)%R_internal(1,2)=Rt*dl/2d0
      cable_geometry_list(cable_geometry_number)%R_internal(2,1)=Rt*dl/2d0
    
      cable_geometry_list(cable_geometry_number)%Tv(1,1)= 1
      cable_geometry_list(cable_geometry_number)%Tv(1,2)=-1
      cable_geometry_list(cable_geometry_number)%Tv(2,1)= 0
      cable_geometry_list(cable_geometry_number)%Tv(2,2)= 1
    
      cable_geometry_list(cable_geometry_number)%Ti(1,1)= 1
      cable_geometry_list(cable_geometry_number)%Ti(1,2)= 0
      cable_geometry_list(cable_geometry_number)%Ti(2,1)= 1
      cable_geometry_list(cable_geometry_number)%Ti(2,2)= 1
      
      cable_geometry_list(cable_geometry_number)%SC(1)=1
      cable_geometry_list(cable_geometry_number)%SC(2)=0

      ALLOCATE( cable_geometry_list(cable_geometry_number)%filter_number(1:n_conductors,1:n_conductors) )
      
      if (cable_geometry_list(cable_geometry_number)%cable_geometry_type.EQ.cable_geometry_type_coaxial) then

        n_filters=0
        if (cable_geometry_list(cable_geometry_number)%n_filters.ne.n_filters) GOTO 9020

        cable_geometry_list(cable_geometry_number)%filter_number(1:n_conductors,1:n_conductors)=0      

      else if(cable_geometry_list(cable_geometry_number)%cable_geometry_type.EQ.cable_geometry_type_FD_coaxial) then

        n_filters=3
        if (cable_geometry_list(cable_geometry_number)%n_filters.ne.n_filters) GOTO 9020

! set inner conductor filter     
        cable_geometry_list(cable_geometry_number)%filter_number(1,1)=1
      
! set shield filter
        cable_geometry_list(cable_geometry_number)%filter_number(2,2)=2
      
! set transfer impedance filters
        cable_geometry_list(cable_geometry_number)%filter_number(1,2)=3
        cable_geometry_list(cable_geometry_number)%filter_number(2,1)=3
      
      end if
      
    else if ((cable_geometry_list(cable_geometry_number)%cable_geometry_type.EQ.cable_geometry_type_ribbon).OR.    &
             (cable_geometry_list(cable_geometry_number)%cable_geometry_type.EQ.cable_geometry_type_FD_ribbon)) then 
    
      n_params=4
      if (cable_geometry_list(cable_geometry_number)%n_parameters.ne.n_params) GOTO 9010     
    
      n_conductors=cable_geometry_list(cable_geometry_number)%n_conductors
      n_external=cable_geometry_list(cable_geometry_number)%n_conductors
      n_shielded=0
      
      cable_geometry_list(cable_geometry_number)%n_conductors         =n_conductors      
      cable_geometry_list(cable_geometry_number)%n_shielded_conductors=n_shielded
      cable_geometry_list(cable_geometry_number)%n_external_conductors=n_external

! SET UP RIBBON CABLE GEOMETRY HERE          

      r_inner =cable_geometry_list(cable_geometry_number)%parameters(1)
      r_dielectric=cable_geometry_list(cable_geometry_number)%parameters(2)
      epsr    =cable_geometry_list(cable_geometry_number)%parameters(3)
      separation=cable_geometry_list(cable_geometry_number)%parameters(4)
      
      ribbon_centre_width=(n_conductors-1)*separation
      ribbon_width=ribbon_centre_width+2d0*r_dielectric
      
      ALLOCATE( cable_geometry_list(cable_geometry_number)%external_conductor_xc(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%external_conductor_yc(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%external_conductor_radius(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%external_dielectric_radius(1:n_external) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%external_dielectric_permittivity(1:n_external) )
      
      cable_geometry_list(cable_geometry_number)%external_conductor_yc(1:n_external)=0d0     
      cable_geometry_list(cable_geometry_number)%external_conductor_radius(1:n_external)=r_inner        		      
      cable_geometry_list(cable_geometry_number)%external_dielectric_radius(1:n_external)=r_dielectric        		      
      cable_geometry_list(cable_geometry_number)%external_dielectric_permittivity(1:n_external)=epsr

      cable_geometry_list(cable_geometry_number)%external_conductor_xc(1)=-ribbon_centre_width/2d0
      do i=2,n_conductors      
        cable_geometry_list(cable_geometry_number)%external_conductor_xc(i)=			&
	   cable_geometry_list(cable_geometry_number)%external_conductor_xc(i-1)+separation	
      end do

      cable_geometry_list(cable_geometry_number)%cable_offset_radius=1.01d0*ribbon_width/2d0	 
				
      ALLOCATE( cable_geometry_list(cable_geometry_number)%L_internal(1:n_conductors,1:n_conductors) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%C_internal(1:n_conductors,1:n_conductors) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%R_internal(1:n_conductors,1:n_conductors) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%Tv(1:n_conductors,1:n_conductors) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%Ti(1:n_conductors,1:n_conductors) )
      ALLOCATE( cable_geometry_list(cable_geometry_number)%SC(1:n_conductors) )

      cable_geometry_list(cable_geometry_number)%L_internal(1:n_conductors,1:n_conductors)=0d0
      cable_geometry_list(cable_geometry_number)%C_internal(1:n_conductors,1:n_conductors)=0d0
      cable_geometry_list(cable_geometry_number)%R_internal(1:n_conductors,1:n_conductors)=0d0
      cable_geometry_list(cable_geometry_number)%Tv(1:n_conductors,1:n_conductors)=0
      cable_geometry_list(cable_geometry_number)%Ti(1:n_conductors,1:n_conductors)=0
      cable_geometry_list(cable_geometry_number)%SC(1:n_conductors)=0
      
      do row=1,n_conductors
	cable_geometry_list(cable_geometry_number)%Tv(row,row)=1
	cable_geometry_list(cable_geometry_number)%Ti(row,row)=1
      end do

      ALLOCATE( cable_geometry_list(cable_geometry_number)%filter_number(1:n_conductors,1:n_conductors) )
      
      if (cable_geometry_list(cable_geometry_number)%cable_geometry_type.EQ.cable_geometry_type_ribbon) then

        n_filters=0
        if (cable_geometry_list(cable_geometry_number)%n_filters.ne.n_filters) GOTO 9020
        cable_geometry_list(cable_geometry_number)%filter_number(1:n_conductors,1:n_conductors)=0          
      
      else if(cable_geometry_list(cable_geometry_number)%cable_geometry_type.EQ.cable_geometry_type_FD_ribbon) then
 
        n_filters=1
        if (cable_geometry_list(cable_geometry_number)%n_filters.ne.n_filters) GOTO 9020
        cable_geometry_list(cable_geometry_number)%filter_number(1:n_conductors,1:n_conductors)=0             
	do row=1,n_conductors
          cable_geometry_list(cable_geometry_number)%filter_number(row,row)=1
	end do
      
      end if
                      
    else
    
      GOTO 9000
      
    end if

     
  end do ! next cable geometry

  CALL write_line('FINISHED: create_cable_LCRG_matrices',0,output_to_screen_flag)
    
  RETURN
    
9000 CALL write_line('Error in create_cable_LCRG_matrices:',0,.TRUE.)
     CALL write_line_integer('No process for cable_geometry_type:',	&
                             cable_geometry_list(cable_geometry_number)%cable_geometry_type,0,.TRUE.)
     STOP
    
9010 CALL write_line('Error in create_cable_LCRG_matrices:',0,.TRUE.)
     CALL write_line_integer('Wrong number of parameters for cable_geometry_type:',	&
                             cable_geometry_list(cable_geometry_number)%cable_geometry_type,0,.TRUE.)
     CALL write_line_integer('Number of parameters should be:',n_params,0,.TRUE.)
     CALL write_line_integer('Number of parameters found:',	&
                              cable_geometry_list(cable_geometry_number)%n_parameters,0,.TRUE.)
     STOP
    
9020 CALL write_line('Error in create_cable_LCRG_matrices:',0,.TRUE.)
     CALL write_line_integer('Wrong number of filters for cable_geometry_type:',	&
                             cable_geometry_list(cable_geometry_number)%cable_geometry_type,0,.TRUE.)
     CALL write_line_integer('Number of filters should be:',n_filters,0,.TRUE.)
     CALL write_line_integer('Number of filters found:',	&
                              cable_geometry_list(cable_geometry_number)%n_filters,0,.TRUE.)
     STOP
 
  
END SUBROUTINE create_cable_LCRG_matrices
