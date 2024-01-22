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
! SUBROUTINE scatter
! SUBROUTINE add_shunt_admittance
! SUBROUTINE add_series_impedance
!
! NAME
!     scatter
!
! DESCRIPTION
!     SCN scatter 
!     
! COMMENTS
!     New idea is to incrementally build a Thevenin equivalent at the
!     cell centre so as to include materials and cables in a unified manner
!     Field output then also becomes straightforward if recorded at this time
!
!  Note 1: cell centre sources have a factor of 0.5 applied which gives consistency
!          between cell and face excitations but leads to discrepancies in cell output at the source
!
!  Note 2: PEC and PMC cells give the correct E and H fields at the cell centres
!          but maybe not the best way to implement it as allows some field through
!
! HISTORY
!
!     started 14/08/2012 CJS
!     revised output process 14/09/2012 CJS
!     revised excitation process 18/08/2012 CJS
!     parallel 23/11/2012 CJS
!     add current source to nodes 26/9/2013 CJS
!     allow hard and soft sources 12/2/2014 CJS
!
!
SUBROUTINE scatter

USE TLM_general
USE TLM_excitation
USE TLM_output
USE TLM_volume_materials
USE PML_module
USE mesh
USE constants

IMPLICIT NONE

! local variables

  integer cx,cy,cz
  
  real*8 	:: Vx,Vy,Vz,Ix,Iy,Iz,V_min,V_max

! SCN link line data 
  real*8 	:: Yx,Yy,Yz,Zx,Zy,Zz
  real*8	:: VEx,VEy,VEz,VHx,VHy,VHz
  
! SCN inductive and capcaitive stub data  
  real*8 	:: Ysx,Ysy,Ysz,Zsx,Zsy,Zsz
  real*8 	:: Vsx,Vsy,Vsz,Vox,Voy,Voz
  real*8 	:: Iox,Ioy,Ioz
  
! SCN electric and magnetic loss data  
  real*8 	:: Gex,Gey,Gez,Rmx,Rmy,Rmz
  
  real*8	:: Jsx,Jsy,Jsz,Msx,Msy,Msz
  
  real*8	:: Z04,Y04

! Cable currents  
  real*8	:: Iwx,Iwy,Iwz
  
! Source currents  
  real*8	:: Isx,Isy,Isz

  integer	:: cell_number
  integer	:: special_cell_count
  
! Material update parameters  
  integer	:: material_number
  integer	:: material_type
  integer	:: filter_number
   
! Output parameters  
  integer	:: output_number
  integer	:: output_cell_number
  integer	:: field_component
  real*8	:: field(6)
  
! Excitation parameters  
  integer	:: excitation_number
  integer	:: excitation_cell_number
  
! Cable parameters  
  integer	:: cable_cell_junction_number
  
  integer	:: hard_source_factor(6)
  
! PML parameters

  integer       :: PML_cell
  integer       :: PML_material_cell
  integer       :: PML_parameter
  real*8        :: sx,sy,sz
  real*8        :: ax,ay,az
  
  real*8        :: Vxt,Vyt,Vzt
  real*8        :: Vyxt,Vzxt,Vxyt,Vzyt,Vyzt,Vxzt
  
  real*8        :: Ixt,Iyt,Izt
  
  real*8        :: last_Ixt,last_Iyt,last_Izt
  real*8        :: last_Ix,last_Iy,last_Iz
  real*8        :: last_Vi(12)
  
  real*8        :: last_Vxt,last_Vyt,last_Vzt
  real*8        :: last_Vyxt,last_Vzxt,last_Vxyt,last_Vzyt,last_Vyzt,last_Vxzt

! new notation for PML in materials
  
  real*8        :: Ixt_m1,Iyt_m1,Izt_m1
  real*8        :: Ix_m1,Iy_m1,Iz_m1
  real*8        :: Vi_m1(12)
  
  real*8        :: Vxt_m1,Vyt_m1,Vzt_m1
  real*8        :: Vyxt_m1,Vzxt_m1,Vxyt_m1,Vzyt_m1,Vyzt_m1,Vxzt_m1
  
  real*8        :: Ixt_m2,Iyt_m2,Izt_m2
  real*8        :: Ix_m2,Iy_m2,Iz_m2
  real*8        :: Vi_m2(12)
  
  real*8        :: Vxt_m2,Vyt_m2,Vzt_m2
  real*8        :: Vyxt_m2,Vzxt_m2,Vxyt_m2,Vzyt_m2,Vyzt_m2,Vxzt_m2
  
  real*8        :: Vx_PML_shunt,Vy_PML_shunt,Vz_PML_shunt
  real*8        :: Vyx_PML_series,Vzx_PML_series,Vxy_PML_series,Vzy_PML_series,Vyz_PML_series,Vxz_PML_series
  
  real*8        :: Vynx_PML,Vypx_PML,Vznx_PML,Vzpx_PML
  real*8        :: Vxny_PML,Vxpy_PML,Vzny_PML,Vzpy_PML  
  real*8        :: Vxnz_PML,Vxpz_PML,Vynz_PML,Vypz_PML
  
  real*8        :: Vosx,Vosx_m1,Vosx_m2
  real*8        :: Vosy,Vosy_m1,Vosy_m2
  real*8        :: Vosz,Vosz_m1,Vosz_m2
  
  real*8        :: Vssx
  real*8        :: Vssy
  real*8        :: Vssz
  
  real*8        :: Vcx,Vcx_m1,Vcx_m2
  real*8        :: Vcy,Vcy_m1,Vcy_m2
  real*8        :: Vcz,Vcz_m1,Vcz_m2

  real*8        :: exp_x,exp_y,exp_z
  real*8        :: exp_xi,exp_yi,exp_zi
  real*8        :: exp_xr,exp_yr,exp_zr
  real*8        :: Csx,Csy,Csz
  real*8        :: Csx2,Csy2,Csz2
  real*8        :: Csy_sz,Csz_sx,Csx_sy
  
  real*8        :: Ax0,Ax1,Ax2, Ay0,Ay1,Ay2, Az0,Az1,Az2
  real*8        :: Bx0,Bx1,Bx2, By0,By1,By2, Bz0,Bz1,Bz2
  real*8        :: Cx0,Cx1,Cx2, Cy0,Cy1,Cy2, Cz0,Cz1,Cz2
  real*8        :: Dx0,Dx1,Dx2, Dy0,Dy1,Dy2, Dz0,Dz1,Dz2
  
  real*8        :: Sx0,Sx1,Sx2, Sy0,Sy1,Sy2, Sz0,Sz1,Sz2 
  real*8        :: Tx0,Tx1,Tx2, Ty0,Ty1,Ty2, Tz0,Tz1,Tz2 

  real*8        :: Uxy0,Uxy1, Uyz0,Uyz1, Uzx0,Uzx1
  real*8        :: Uyx0,Uyx1, Uzy0,Uzy1, Uxz0,Uxz1
  real*8        :: Wxy0,Wxy1, Wyz0,Wyz1, Wzx0,Wzx1
  real*8        :: Wyx0,Wyx1, Wzy0,Wzy1, Wxz0,Wxz1
  
  real*8        :: GpY2
  real*8        :: ZpR
  
  real*8        :: Px0,Px1,Px2, Py0,Py1,Py2, Pz0,Pz1,Pz2 
  real*8        :: Qx0,Qx1,Qx2, Qy0,Qy1,Qy2, Qz0,Qz1,Qz2 
  
  real*8        :: Vsx_m1,Vsx_m2
  real*8        :: Vsy_m1,Vsy_m2
  real*8        :: Vsz_m1,Vsz_m2
  
  real*8        :: Vsx_temp,Vsy_temp,Vsz_temp
  
  real*8        :: loss_dist_inc,loss_dist_ref
  
  logical :: op_flag

! START
  
  CALL write_line('CALLED: scatter',0,timestepping_output_to_screen_flag)
  
  op_flag=.FALSE.
  
  Z04=Z0*4d0
  Y04=Y0*4d0

  cell_number=0
  output_cell_number=0
  excitation_cell_number=0
  
  do cz=nz1,nz2
    do cy=1,ny
      do cx=1,nx
      
        cell_number=cell_number+1
	
	if (cell_centre_update_code(cell_number).eq.0) then   ! free space update, no output
	
          Vx=(	  V(Vx_ymin,cx,cy,cz)+V(Vx_ymax,cx,cy,cz)	&
	  	 +V(Vx_zmin,cx,cy,cz)+V(Vx_zmax,cx,cy,cz) )/ 2d0       
	      
          Vy=(	  V(Vy_xmin,cx,cy,cz)+V(Vy_xmax,cx,cy,cz)	&
		 +V(Vy_zmin,cx,cy,cz)+V(Vy_zmax,cx,cy,cz) )/ 2d0       
	      
          Vz=(	  V(Vz_xmin,cx,cy,cz)+V(Vz_xmax,cx,cy,cz)	&
		 +V(Vz_ymin,cx,cy,cz)+V(Vz_ymax,cx,cy,cz) )/ 2d0       
	          
	  Ix=( V(Vz_ymax,cx,cy,cz)-V(Vy_zmax,cx,cy,cz)	&
	      -V(Vz_ymin,cx,cy,cz)+V(Vy_zmin,cx,cy,cz) )/ (2d0*Z0)    
	
	  Iy=( V(Vx_zmax,cx,cy,cz)-V(Vz_xmax,cx,cy,cz)	&
	      -V(Vx_zmin,cx,cy,cz)+V(Vz_xmin,cx,cy,cz) )/ (2d0*Z0)    
	
	  Iz=( V(Vy_xmax,cx,cy,cz)-V(Vx_ymax,cx,cy,cz)	&
	      -V(Vy_xmin,cx,cy,cz)+V(Vx_ymin,cx,cy,cz) )/ (2d0*Z0)   
	
! Calculate Scattered voltage pulses into link lines	

          V_min=V(Vx_ymin,cx,cy,cz)
          V_max=V(Vx_ymax,cx,cy,cz)     
	  V(Vx_ymin,cx,cy,cz)=Vx-Z0*Iz-V_max
	  V(Vx_ymax,cx,cy,cz)=Vx+Z0*Iz-V_min

          V_min=V(Vx_zmin,cx,cy,cz)
          V_max=V(Vx_zmax,cx,cy,cz)     
	  V(Vx_zmin,cx,cy,cz)=Vx+Z0*Iy-V_max
	  V(Vx_zmax,cx,cy,cz)=Vx-Z0*Iy-V_min
	
          V_min=V(Vy_zmin,cx,cy,cz)
          V_max=V(Vy_zmax,cx,cy,cz)     
	  V(Vy_zmin,cx,cy,cz)=Vy-Z0*Ix-V_max
	  V(Vy_zmax,cx,cy,cz)=Vy+Z0*Ix-V_min
	
          V_min=V(Vy_xmin,cx,cy,cz)
          V_max=V(Vy_xmax,cx,cy,cz)     
	  V(Vy_xmin,cx,cy,cz)=Vy+Z0*Iz-V_max
	  V(Vy_xmax,cx,cy,cz)=Vy-Z0*Iz-V_min
	
          V_min=V(Vz_xmin,cx,cy,cz)
          V_max=V(Vz_xmax,cx,cy,cz)     
	  V(Vz_xmin,cx,cy,cz)=Vz-Z0*Iy-V_max
	  V(Vz_xmax,cx,cy,cz)=Vz+Z0*Iy-V_min
	
          V_min=V(Vz_ymin,cx,cy,cz)
          V_max=V(Vz_ymax,cx,cy,cz)     
	  V(Vz_ymin,cx,cy,cz)=Vz+Z0*Ix-V_max
	  V(Vz_ymax,cx,cy,cz)=Vz-Z0*Ix-V_min        

	else
! cell_update_code points to arrays which tell us about material number and output_number

	  special_cell_count=cell_centre_update_code(cell_number)
          
          PML_cell=cell_update_code_to_PML_data(special_cell_count)
	  
          if (PML_cell.EQ.0) then
! not PML update
          
	    material_number=cell_update_code_to_material_data(special_cell_count,1)
	    if (material_number.ne.0) then
	      material_type=volume_material_list(material_number)%type
	    else
	      material_type=0
	    end if
	    cable_cell_junction_number=cell_update_code_to_cable_cell_number(special_cell_count)
	    output_number=cell_update_code_to_output_number(special_cell_count)
	    excitation_number=cell_update_code_to_excitation_number(special_cell_count)

! Calculate internal V and I with no stubs
! Calculate equivalent circuit parameters in x, y and z directions for Voltage calculation
            VEx=(	V(Vx_ymin,cx,cy,cz)+V(Vx_ymax,cx,cy,cz)+V(Vx_zmin,cx,cy,cz)+V(Vx_zmax,cx,cy,cz) )/2d0   	      
            VEy=(	V(Vy_xmin,cx,cy,cz)+V(Vy_xmax,cx,cy,cz)+V(Vy_zmin,cx,cy,cz)+V(Vy_zmax,cx,cy,cz) )/2d0     	      
            VEz=(	V(Vz_xmin,cx,cy,cz)+V(Vz_xmax,cx,cy,cz)+V(Vz_ymin,cx,cy,cz)+V(Vz_ymax,cx,cy,cz) )/2d0    
	    Yx=Y04
	    Yy=Y04
	    Yz=Y04
	          
! Calculate equivalent circuit parameters in x, y and z directions for Current calculation
	    VHx=2d0*( V(Vz_ymax,cx,cy,cz)-V(Vy_zmax,cx,cy,cz)-V(Vz_ymin,cx,cy,cz)+V(Vy_zmin,cx,cy,cz) )	
	    VHy=2d0*( V(Vx_zmax,cx,cy,cz)-V(Vz_xmax,cx,cy,cz)-V(Vx_zmin,cx,cy,cz)+V(Vz_xmin,cx,cy,cz) )  	
	    VHz=2d0*( V(Vy_xmax,cx,cy,cz)-V(Vx_ymax,cx,cy,cz)-V(Vy_xmin,cx,cy,cz)+V(Vx_ymin,cx,cy,cz) ) 
            Zx=Z04
	    Zy=Z04
	    Zz=Z04

            if (material_type.eq.volume_material_type_PEC) then  
! Set internal voltages to zero	  

	      VEx=0d0
              VEy=0d0
              VEz=0d0
	    
            else if (material_type.eq.volume_material_type_PMC) then
! Set internal currents to zero	  
	  
              VHx=0d0
              VHy=0d0
              VHz=0d0
	  
            else if (material_type.eq.volume_material_type_DISPERSIVE) then
	  
              if (volume_material_list(material_number)%eps_filter_exists) then
	  
	        filter_number=cell_update_code_to_material_data(special_cell_count,2)

! Admittance and voltage source of capacitive (open circuit) stubs	    
                Ysx=volume_material_list(material_number)%Ys_eps_f
	        Ysy=volume_material_list(material_number)%Ys_eps_f
	        Ysz=volume_material_list(material_number)%Ys_eps_f
	    
! eps_x slow response
                call timeshift_Zfilter(volume_material_Zs_eps_filter_data(filter_number  ))
                call evaluate_Zfilter(volume_material_list(material_number  )%Zs_eps_Z,		&
	                            volume_material_Zs_eps_filter_data(filter_number),0d0)
                Vox=volume_material_Zs_eps_filter_data(filter_number  )%f
! eps_y slow response
                call timeshift_Zfilter(volume_material_Zs_eps_filter_data(filter_number+1))
                call evaluate_Zfilter(volume_material_list(material_number)%Zs_eps_Z,		&
	                            volume_material_Zs_eps_filter_data(filter_number+1),0d0)
                Voy=volume_material_Zs_eps_filter_data(filter_number+1)%f
! eps_z slow response
                call timeshift_Zfilter(volume_material_Zs_eps_filter_data(filter_number+2))
                call evaluate_Zfilter(volume_material_list(material_number)%Zs_eps_Z,		&
	                            volume_material_Zs_eps_filter_data(filter_number+2),0d0)
                Voz=volume_material_Zs_eps_filter_data(filter_number+2)%f
   
	        CALL add_shunt_admittance(VEx,Yx,Vox,Ysx,VEx,Yx)
	        CALL add_shunt_admittance(VEy,Yy,Voy,Ysy,VEy,Yy)
	        CALL add_shunt_admittance(VEz,Yz,Voz,Ysz,VEz,Yz)
	    
	      end if ! permittivity filter exists

! Electric loss 	    
              Gex=volume_material_list(material_number)%Ge
	      Gey=volume_material_list(material_number)%Ge
	      Gez=volume_material_list(material_number)%Ge
	    
	      CALL add_shunt_admittance(VEx,Yx,0d0,Gex,VEx,Yx)
	      CALL add_shunt_admittance(VEy,Yy,0d0,Gey,VEy,Yy)
	      CALL add_shunt_admittance(VEz,Yz,0d0,Gez,VEz,Yz)
	  
              if (volume_material_list(material_number)%mu_filter_exists) then
	  
	        filter_number=cell_update_code_to_material_data(special_cell_count,2)
	    
! Impedance and voltage source of inductive (short circuit) stubs	    
                Zsx=volume_material_list(material_number)%Zs_mu_f
	        Zsy=volume_material_list(material_number)%Zs_mu_f
	        Zsz=volume_material_list(material_number)%Zs_mu_f
	        
! mu_x slow response
                call timeshift_Zfilter(volume_material_Zs_mu_filter_data(filter_number  ))
                call evaluate_Zfilter(volume_material_list(material_number)%Zs_mu_Z,		&
	                            volume_material_Zs_mu_filter_data(filter_number  ),0d0)
                Vsx=volume_material_Zs_mu_filter_data(filter_number  )%f
! mu_y slow response
                call timeshift_Zfilter(volume_material_Zs_mu_filter_data(filter_number+1))
                call evaluate_Zfilter(volume_material_list(material_number)%Zs_mu_Z,		&
	                            volume_material_Zs_mu_filter_data(filter_number+1),0d0)
                Vsy=volume_material_Zs_mu_filter_data(filter_number+1)%f
! mu_z slow response
                call timeshift_Zfilter(volume_material_Zs_mu_filter_data(filter_number+2))
                call evaluate_Zfilter(volume_material_list(material_number)%Zs_mu_Z,		&
	                            volume_material_Zs_mu_filter_data(filter_number+2),0d0)
                Vsz=volume_material_Zs_mu_filter_data(filter_number+2)%f
	    
	        CALL add_series_impedance(VHx,Zx,Vsx,Zsx,VHx,Zx)
	        CALL add_series_impedance(VHy,Zy,Vsy,Zsy,VHy,Zy)
	        CALL add_series_impedance(VHz,Zz,Vsz,Zsz,VHz,Zz)
	    
	      end if ! permeability filter exists

! Magnetic loss 	    
              Rmx=volume_material_list(material_number)%Rm
	      Rmy=volume_material_list(material_number)%Rm
	      Rmz=volume_material_list(material_number)%Rm
	    
	      CALL add_series_impedance(VHx,Zx,0d0,Rmx,VHx,Zx)
	      CALL add_series_impedance(VHy,Zy,0d0,Rmy,VHy,Zy)
	      CALL add_series_impedance(VHz,Zz,0d0,Rmz,VHz,Zz)  
	    
	    end if  ! DISPERSIVE material
                    
 
! Electric and Magnetic current density source terms are added as follows (may be required...)
!	    VEx=VEx+Jsx*const1/Yx
!	    VEy=VEy+Jsy*const1/Yy
!	    VEz=VEz+Jsz*const1/Yz

!	    VHx=VHx-Msx*const2
!	    VHy=VHy-Msy*const2
!	    VHz=VHz-Msz*const2

! calculate internal Voltage and current 
            Vx=VEx
            Vy=VEy
            Vz=VEz
	  
            Ix=VHx/Zx
            Iy=VHy/Zy
            Iz=VHz/Zz

! Excitation		
            if (excitation_number.ne.0) then
	  
	      excitation_cell_number=excitation_cell_number+1
	      field(1:6)=cell_excitation_field(excitation_cell_number,1:6)
	      hard_source_factor(1:6)=cell_excitation_type(excitation_cell_number,1:6)

! work out the current source as the integral of the current density over the cell area	    
	      Isx=cell_excitation_field(excitation_cell_number,7)*dl*dl
	      Isy=cell_excitation_field(excitation_cell_number,8)*dl*dl
	      Isz=cell_excitation_field(excitation_cell_number,9)*dl*dl

! note: soft source only at the moment...
	      Vx=hard_source_factor(1)*Vx-field(Ex)*dl/2d0+Isx/Yx
	      Vy=hard_source_factor(2)*Vy-field(Ey)*dl/2d0+Isy/Yy
	      Vz=hard_source_factor(3)*Vz-field(Ez)*dl/2d0+Isz/Yz
	      Ix=hard_source_factor(4)*Ix+field(Hx)*dl/2d0
	      Iy=hard_source_factor(5)*Iy+field(Hy)*dl/2d0
	      Iz=hard_source_factor(6)*Iz+field(Hz)*dl/2d0
	    
!	      write(*,*)'finished excitation'
	    
	    end if ! excitation_number.ne.0
	  
! Cell centre cable update. Note we now have Vx, Vy, Vz, Zx, Zy, Zz to perform the cell centre cable update 	  
	    if (cable_cell_junction_number.NE.0) then
	  	  
	      CALL cell_centre_cable_junction( cable_cell_junction_number, Vx,Vy,Vz, Yx,Yy,Yz, Iwx,Iwy,Iwz )

! calculate the difference in the node voltage in the x, y and z directions due to the cable common mode current
	      Vx=Vx-Iwx/Yx
	      Vy=Vy-Iwy/Yy
	      Vz=Vz-Iwz/Yz

	    end if ! cable_number.NE.0
	  
! Output		
            if (output_number.ne.0) then
	    
	      output_cell_number=output_cell_number+1
	  
! calculate cell centre fields		     
	      field(Ex)=-Vx/dl
	      field(Ey)=-Vy/dl
	      field(Ez)=-Vz/dl
	
	      field(Hx)= Ix/dl
	      field(Hy)= Iy/dl
	      field(Hz)= Iz/dl
	  
              cell_output_field(output_cell_number,1:6)=field(1:6)
	  
	    end if ! output_number.ne.0

! If this is a dispersive material then update the filter functions
            if (material_type.eq.volume_material_type_DISPERSIVE) then
	  	  
              if (volume_material_list(material_number)%eps_filter_exists) then
	    
	        filter_number=cell_update_code_to_material_data(special_cell_count,2)

! calculate capacitor stub currents
                Iox=(Vx-Vox)*Ysx
                Ioy=(Vy-Voy)*Ysy
                Ioz=(Vz-Voz)*Ysz
	    
                call evaluate_Zfilter(volume_material_list(material_number)%Zs_eps_Z,		&
	                            volume_material_Zs_eps_filter_data(filter_number  ),Iox)
                call evaluate_Zfilter(volume_material_list(material_number)%Zs_eps_Z,		&
	                            volume_material_Zs_eps_filter_data(filter_number+1),Ioy)
                call evaluate_Zfilter(volume_material_list(material_number)%Zs_eps_Z,		&
	                            volume_material_Zs_eps_filter_data(filter_number+2),Ioz)
				  
	      end if ! permittivity filter exists
	    
              if (volume_material_list(material_number)%mu_filter_exists) then

! we already know the inductor stub currents	    
                call evaluate_Zfilter(volume_material_list(material_number)%Zs_mu_Z,		&
	                            volume_material_Zs_mu_filter_data(filter_number  ),-Ix)
                call evaluate_Zfilter(volume_material_list(material_number)%Zs_mu_Z,		&
	                            volume_material_Zs_mu_filter_data(filter_number+1),-Iy)
                call evaluate_Zfilter(volume_material_list(material_number)%Zs_mu_Z,		&
	                            volume_material_Zs_mu_filter_data(filter_number+2),-Iz)
				    
	      end if ! permeability filter exists
	      
            end if ! DISPERSIVE material
	
! Calculate Scattered voltage pulses into link lines	

            V_min=V(Vx_ymin,cx,cy,cz)
            V_max=V(Vx_ymax,cx,cy,cz)     
	    V(Vx_ymin,cx,cy,cz)=Vx-Z0*Iz-V_max
	    V(Vx_ymax,cx,cy,cz)=Vx+Z0*Iz-V_min

            V_min=V(Vx_zmin,cx,cy,cz)
            V_max=V(Vx_zmax,cx,cy,cz)     
	    V(Vx_zmin,cx,cy,cz)=Vx+Z0*Iy-V_max
	    V(Vx_zmax,cx,cy,cz)=Vx-Z0*Iy-V_min
	
            V_min=V(Vy_zmin,cx,cy,cz)
            V_max=V(Vy_zmax,cx,cy,cz)     
	    V(Vy_zmin,cx,cy,cz)=Vy-Z0*Ix-V_max
	    V(Vy_zmax,cx,cy,cz)=Vy+Z0*Ix-V_min
	
            V_min=V(Vy_xmin,cx,cy,cz)
            V_max=V(Vy_xmax,cx,cy,cz)     
	    V(Vy_xmin,cx,cy,cz)=Vy+Z0*Iz-V_max
	    V(Vy_xmax,cx,cy,cz)=Vy-Z0*Iz-V_min
	
            V_min=V(Vz_xmin,cx,cy,cz)
            V_max=V(Vz_xmax,cx,cy,cz)     
	    V(Vz_xmin,cx,cy,cz)=Vz-Z0*Iy-V_max
	    V(Vz_xmax,cx,cy,cz)=Vz+Z0*Iy-V_min
	
            V_min=V(Vz_ymin,cx,cy,cz)
            V_max=V(Vz_ymax,cx,cy,cz)     
	    V(Vz_ymin,cx,cy,cz)=Vz+Z0*Ix-V_max
  	    V(Vz_ymax,cx,cy,cz)=Vz-Z0*Ix-V_min        
          
          else

! We are in a PML cell
	    material_number=cell_update_code_to_material_data(special_cell_count,1)
	    if (material_number.ne.0) then
	      material_type=volume_material_list(material_number)%type
	    else
	      material_type=0
	    end if

            if (material_type.eq.volume_material_type_PEC) then  

              write(*,*)'ERROR: PEC volume material in PML'
	      STOP 1
	    
            else if (material_type.eq.volume_material_type_PMC) then

              write(*,*)'ERROR: PMC volume material in PML'
	      STOP 1

            else if (material_type.eq.volume_material_type_DISPERSIVE) then
	      	       
#include "scatter_PML_material.F90"
 	      
	    else
 
#include "scatter_PML.F90"

            end if      ! material or free space PML update 
  
          end if  ! PML
 
        end if ! cell_update_code      
	
      end do  ! next z cell
    end do    ! next y cell
  end do      ! next x cell
  
  CALL write_line('FINISHED: scatter',0,timestepping_output_to_screen_flag)

  RETURN

END SUBROUTINE scatter
!
! NAME
!     add_shunt_admittance
!
! DESCRIPTION
!     add_shunt_admittance
!     
! COMMENTS
!     
! HISTORY
!
!     started 5/09/2012 CJS
!
!
SUBROUTINE add_shunt_admittance(V1,Y1,V2,Y2,V3,Y3)

IMPLICIT NONE

  real*8 	:: V1,Y1,V2,Y2,V3,Y3

! local variables
  
  real*8 	:: V,Y

! START

  V=(V1*Y1+V2*Y2)/(Y1+Y2)
  Y=Y1+Y2

  V3=V
  Y3=Y

  RETURN
  
END SUBROUTINE add_shunt_admittance
!
! NAME
!     add_series_impedance
!
! DESCRIPTION
!     add_series_impedance
!     
! COMMENTS
!     
! HISTORY
!
!     started 5/09/2012 CJS
!
!
SUBROUTINE add_series_impedance(V1,Z1,V2,Z2,V3,Z3)

IMPLICIT NONE

  real*8 	:: V1,Z1,V2,Z2,V3,Z3

! local variables
  
  real*8 	:: V,Z

! START

  V=V1+V2
  Z=Z1+Z2

  V3=V
  Z3=Z

  RETURN
  
END SUBROUTINE add_series_impedance
