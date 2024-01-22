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
! scatter_PML.F90: File to be included in scatter.F90 to handle the PML cell centre update
!
!  Started 3/10/2019

 	    output_number=cell_update_code_to_output_number(special_cell_count)

! STAGE 1: get the parameters for this PML cell

            PML_parameter=PML_cell_data(PML_cell)%PML_parameter_array_pos
	    
            sx=PML_parameters(PML_parameter)%sx
            sy=PML_parameters(PML_parameter)%sy
            sz=PML_parameters(PML_parameter)%sz
   
            ax=PML_parameters(PML_parameter)%ax
            ay=PML_parameters(PML_parameter)%ay
            az=PML_parameters(PML_parameter)%az
 
! Loss factor for propagation half a cell in x, y and z        
            exp_x=PML_parameters(PML_parameter)%exp_x
            exp_y=PML_parameters(PML_parameter)%exp_y
            exp_z=PML_parameters(PML_parameter)%exp_z

! Loss factor as applied to incident fields            
            exp_xi=exp_x  
            exp_yi=exp_y  
            exp_zi=exp_z  

! Loss factor as applied to scattered fields
            exp_xr=exp_x  
            exp_yr=exp_y  
            exp_zr=exp_z  	   
          
            Csx=PML_parameters(PML_parameter)%Csx
            Csy=PML_parameters(PML_parameter)%Csy
            Csz=PML_parameters(PML_parameter)%Csz
          
            Csx2=PML_parameters(PML_parameter)%Csx2
            Csy2=PML_parameters(PML_parameter)%Csy2
            Csz2=PML_parameters(PML_parameter)%Csz2
          
            Csy_sz=PML_parameters(PML_parameter)%Csy_sz
            Csz_sx=PML_parameters(PML_parameter)%Csz_sx
            Csx_sy=PML_parameters(PML_parameter)%Csx_sy
            
!            write(*,'(A,I3,A,3I4,A,3ES12.3)')'rank=',rank,' cx, cy,cz=',cx,cy,cz,' sx, sy, sz=',sx,sy,sz
                       
! STAGE 2: Retrieve voltages and currents from the last timestep

            last_Ix=PML_cell_data(PML_cell)%Ix
            last_Iy=PML_cell_data(PML_cell)%Iy
            last_Iz=PML_cell_data(PML_cell)%Iz
            
            last_Vi(1:12)=PML_cell_data(PML_cell)%Vi(1:12)
  
            last_Vxt=PML_cell_data(PML_cell)%Vxt
            last_Vyt=PML_cell_data(PML_cell)%Vyt
            last_Vzt=PML_cell_data(PML_cell)%Vzt
            
            last_Vyxt=PML_cell_data(PML_cell)%Vyxt
            last_Vzxt=PML_cell_data(PML_cell)%Vzxt
            last_Vxyt=PML_cell_data(PML_cell)%Vxyt
            last_Vzyt=PML_cell_data(PML_cell)%Vzyt
            last_Vyzt=PML_cell_data(PML_cell)%Vyzt
            last_Vxzt=PML_cell_data(PML_cell)%Vxzt
       
! STAGE 3: Propagate incident voltages half a cell with exponential loss for dt/2 (equation 33 for half a cell propagation)
            
            V(Vynx,cx,cy,cz)=V(Vynx,cx,cy,cz)*exp_yi
            V(Vypx,cx,cy,cz)=V(Vypx,cx,cy,cz)*exp_yi
            V(Vznx,cx,cy,cz)=V(Vznx,cx,cy,cz)*exp_zi
            V(Vzpx,cx,cy,cz)=V(Vzpx,cx,cy,cz)*exp_zi
            
            V(Vxny,cx,cy,cz)=V(Vxny,cx,cy,cz)*exp_xi
            V(Vxpy,cx,cy,cz)=V(Vxpy,cx,cy,cz)*exp_xi
            V(Vzny,cx,cy,cz)=V(Vzny,cx,cy,cz)*exp_zi
            V(Vzpy,cx,cy,cz)=V(Vzpy,cx,cy,cz)*exp_zi
                                   
            V(Vxnz,cx,cy,cz)=V(Vxnz,cx,cy,cz)*exp_xi
            V(Vxpz,cx,cy,cz)=V(Vxpz,cx,cy,cz)*exp_xi
            V(Vynz,cx,cy,cz)=V(Vynz,cx,cy,cz)*exp_yi
            V(Vypz,cx,cy,cz)=V(Vypz,cx,cy,cz)*exp_yi

! STAGE 4: Calculate internal V and I with no stubs

! Voltages from shunt circuits, equation 22
            Vx=( V(Vynx,cx,cy,cz)+V(Vypx,cx,cy,cz)+V(Vznx,cx,cy,cz)+V(Vzpx,cx,cy,cz) )/2d0                 
            Vy=( V(Vxny,cx,cy,cz)+V(Vxpy,cx,cy,cz)+V(Vzny,cx,cy,cz)+V(Vzpy,cx,cy,cz) )/2d0                 
            Vz=( V(Vxnz,cx,cy,cz)+V(Vxpz,cx,cy,cz)+V(Vynz,cx,cy,cz)+V(Vypz,cx,cy,cz) )/2d0    
	          
! Currents from series circuits, equation 25
	    Ix=( V(Vypz,cx,cy,cz)-V(Vzpy,cx,cy,cz)-V(Vynz,cx,cy,cz)+V(Vzny,cx,cy,cz) )/(2d0*Z0)
	    Iy=( V(Vzpx,cx,cy,cz)-V(Vxpz,cx,cy,cz)-V(Vznx,cx,cy,cz)+V(Vxnz,cx,cy,cz) )/(2d0*Z0)        
	    Iz=( V(Vxpy,cx,cy,cz)-V(Vypx,cx,cy,cz)-V(Vxny,cx,cy,cz)+V(Vynx,cx,cy,cz) )/(2d0*Z0)

! STAGE 5: Calculate VPML_shunt terms, equation 23

! i=z j=x k=y
            Vx_PML_shunt=ax*(-(last_Vi(Vznx)+last_Vi(Vzpx)+last_Vi(Vynx)+last_Vi(Vypx))/2d0                &
                             +Csy*( (V(Vznx,cx,cy,cz)+V(Vzpx,cx,cy,cz)) + (last_Vi(Vznx)+last_Vi(Vzpx)) )    &
                             +Csz*( (V(Vynx,cx,cy,cz)+V(Vypx,cx,cy,cz)) + (last_Vi(Vynx)+last_Vi(Vypx)) )    &
                             +Csy_sz*last_Vxt )
! i=x j=y k=z           
            Vy_PML_shunt=ay*(-(last_Vi(Vxny)+last_Vi(Vxpy)+last_Vi(Vzny)+last_Vi(Vzpy))/2d0                &
                             +Csz*( (V(Vxny,cx,cy,cz)+V(Vxpy,cx,cy,cz)) + (last_Vi(Vxny)+last_Vi(Vxpy)) )    &
                             +Csx*( (V(Vzny,cx,cy,cz)+V(Vzpy,cx,cy,cz)) + (last_Vi(Vzny)+last_Vi(Vzpy)) )    &
                             +Csz_sx*last_Vyt )
! i=y j=z k=x                       
            Vz_PML_shunt=az*(-(last_Vi(Vynz)+last_Vi(Vypz)+last_Vi(Vxnz)+last_Vi(Vxpz))/2d0                &
                             +Csx*( (V(Vynz,cx,cy,cz)+V(Vypz,cx,cy,cz)) + (last_Vi(Vynz)+last_Vi(Vypz)) )    &
                             +Csy*( (V(Vxnz,cx,cy,cz)+V(Vxpz,cx,cy,cz)) + (last_Vi(Vxnz)+last_Vi(Vxpz)) )    &
                             +Csx_sy*last_Vzt )           

! STAGE 6: Stretched Vxt, Vyt, Vzt, equation 21

            Vxt=ax*Vx+Vx_PML_shunt
            
            Vyt=ay*Vy+Vy_PML_shunt
            
            Vzt=az*Vz+Vz_PML_shunt

! STAGE 7a: Calculate Ixt, Iyt, Izt terms, equation 24 ! Only required for H field output
            if (output_number.ne.0) then
  
              last_Ixt=PML_cell_data(PML_cell)%Ixt    ! Only required for H field output
              last_Iyt=PML_cell_data(PML_cell)%Iyt
              last_Izt=PML_cell_data(PML_cell)%Izt

              Ixt=ax*(Ix-last_Ix + Csx2*(Ix+last_Ix) + Csy_sz*last_Ixt )

              Iyt=ay*(Iy-last_Iy + Csy2*(Iy+last_Iy) + Csz_sx*last_Iyt )

              Izt=az*(Iz-last_Iz + Csz2*(Iz+last_Iz) + Csx_sy*last_Izt )

            end if
            
! STAGE 7: Calculate VPML_series terms, equation 28

! i=x j=y k=z
            Vxy_PML_series=az*(  Z0*( -last_Iz+Csy2*(Iz+last_Iz) )     &
                                +Csx_sy*last_Vxyt                   )
! i=y j=x k=z
            Vyx_PML_series=az*(  Z0*( -last_Iz+Csx2*(Iz+last_Iz) )     &      
                                +Csx_sy*last_Vyxt                   )
! i=z j=x k=y
            Vzx_PML_series=ay*(  Z0*( -last_Iy+Csx2*(Iy+last_Iy) )     &
                                +Csz_sx*last_Vzxt                   )
! i=x j=z k=y
            Vxz_PML_series=ay*(  Z0*( -last_Iy+Csz2*(Iy+last_Iy) )     &
                                +Csz_sx*last_Vxzt                   )
! i=y j=z k=x
            Vyz_PML_series=ax*(  Z0*( -last_Ix+Csz2*(Ix+last_Ix) )     &
                                +Csy_sz*last_Vyzt                   )
! i=z j=y k=x
            Vzy_PML_series=ax*(  Z0*( -last_Ix+Csy2*(Ix+last_Ix) )     &
                                +Csy_sz*last_Vzyt                   )

! STAGE 8: Stretched Vyxt, Vzxt, Vxyt, Vzyt, Vyzt, Vxzt, equation 27

                                         
! i=x j=y k=z
            Vxyt=az*Iz*Z0+Vxy_PML_series

! i=y j=x k=z
            Vyxt=az*Iz*Z0+Vyx_PML_series
                                         
! i=z j=x k=y
            Vzxt=ay*Iy*Z0+Vzx_PML_series
                                         
! i=x j=z k=y
            Vxzt=ay*Iy*Z0+Vxz_PML_series
                                         
! i=y j=z k=x
            Vyzt=ax*Ix*Z0+Vyz_PML_series
                                         
! i=z j=y k=x
            Vzyt=ax*Ix*Z0+Vzy_PML_series
                    
! STAGE 9:  Save voltages and currents for the next timestep

            PML_cell_data(PML_cell)%Ix=Ix
            PML_cell_data(PML_cell)%Iy=Iy
            PML_cell_data(PML_cell)%Iz=Iz
            
            PML_cell_data(PML_cell)%Vi(1:12)=V(1:12,cx,cy,cz)
  
            PML_cell_data(PML_cell)%Vxt=Vxt
            PML_cell_data(PML_cell)%Vyt=Vyt
            PML_cell_data(PML_cell)%Vzt=Vzt
  
            if (output_number.ne.0) then
              PML_cell_data(PML_cell)%Ixt=Ixt   ! Only required for output
              PML_cell_data(PML_cell)%Iyt=Iyt
              PML_cell_data(PML_cell)%Izt=Izt
            end if
            
            PML_cell_data(PML_cell)%Vxyt=Vxyt
            PML_cell_data(PML_cell)%Vyxt=Vyxt
            
            PML_cell_data(PML_cell)%Vzxt=Vzxt
            PML_cell_data(PML_cell)%Vxzt=Vxzt
            
            PML_cell_data(PML_cell)%Vzyt=Vzyt
            PML_cell_data(PML_cell)%Vyzt=Vyzt

! STAGE 10: Calculate V_PML terms from equation 34

            Vynx_PML=Vx_PML_shunt-Vyx_PML_series
            Vypx_PML=Vx_PML_shunt+Vyx_PML_series
            
            Vznx_PML=Vx_PML_shunt+Vzx_PML_series
            Vzpx_PML=Vx_PML_shunt-Vzx_PML_series
            
            Vxny_PML=Vy_PML_shunt+Vxy_PML_series
            Vxpy_PML=Vy_PML_shunt-Vxy_PML_series
            
            Vzny_PML=Vy_PML_shunt-Vzy_PML_series
            Vzpy_PML=Vy_PML_shunt+Vzy_PML_series
            
            Vynz_PML=Vz_PML_shunt+Vyz_PML_series
            Vypz_PML=Vz_PML_shunt-Vyz_PML_series
            
            Vxnz_PML=Vz_PML_shunt-Vxz_PML_series
            Vxpz_PML=Vz_PML_shunt+Vxz_PML_series


! STAGE 11: Calculate scattered voltages from equation 33

            V_max=V(Vypx,cx,cy,cz)
            V_min=V(Vynx,cx,cy,cz)    
	    V(Vynx,cx,cy,cz)=ax*Vx-az*Z0*Iz-V_max+Vynx_PML
	    V(Vypx,cx,cy,cz)=ax*Vx+az*Z0*Iz-V_min+Vypx_PML

            V_max=V(Vzpx,cx,cy,cz)
            V_min=V(Vznx,cx,cy,cz)    
	    V(Vznx,cx,cy,cz)=ax*Vx+ay*Z0*Iy-V_max+Vznx_PML
	    V(Vzpx,cx,cy,cz)=ax*Vx-ay*Z0*Iy-V_min+Vzpx_PML
	
            V_max=V(Vzpy,cx,cy,cz)
            V_min=V(Vzny,cx,cy,cz)    
	    V(Vzny,cx,cy,cz)=ay*Vy-ax*Z0*Ix-V_max+Vzny_PML
	    V(Vzpy,cx,cy,cz)=ay*Vy+ax*Z0*Ix-V_min+Vzpy_PML
	
            V_max=V(Vxpy,cx,cy,cz)
            V_min=V(Vxny,cx,cy,cz)    
	    V(Vxny,cx,cy,cz)=ay*Vy+az*Z0*Iz-V_max+Vxny_PML
	    V(Vxpy,cx,cy,cz)=ay*Vy-az*Z0*Iz-V_min+Vxpy_PML
	
            V_max=V(Vxpz,cx,cy,cz)
            V_min=V(Vxnz,cx,cy,cz)    
	    V(Vxnz,cx,cy,cz)=az*Vz-ay*Z0*Iy-V_max+Vxnz_PML
	    V(Vxpz,cx,cy,cz)=az*Vz+ay*Z0*Iy-V_min+Vxpz_PML
	
            V_max=V(Vypz,cx,cy,cz)
            V_min=V(Vynz,cx,cy,cz)    
	    V(Vynz,cx,cy,cz)=az*Vz+ax*Z0*Ix-V_max+Vynz_PML
  	    V(Vypz,cx,cy,cz)=az*Vz-ax*Z0*Ix-V_min+Vypz_PML     
         
	  
! STAGE 12: Output		
            
            if (output_number.ne.0) then
	    
	      output_cell_number=output_cell_number+1

! calculate cell centre fields
	      field(Ex)=-Vxt/dl
	      field(Ey)=-Vyt/dl
	      field(Ez)=-Vzt/dl
	
	      field(Hx)= Ixt/dl
	      field(Hy)= Iyt/dl
	      field(Hz)= Izt/dl
	  
              cell_output_field(output_cell_number,1:6)=field(1:6)
	  
	    end if ! output_number.ne.0
       
! STAGE 13: Propagate scattered voltages half a cell with exponential loss for dt/2 (equation 33 for half a cell propagation)
            
            V(Vynx,cx,cy,cz)=V(Vynx,cx,cy,cz)*exp_yr
            V(Vypx,cx,cy,cz)=V(Vypx,cx,cy,cz)*exp_yr
            V(Vznx,cx,cy,cz)=V(Vznx,cx,cy,cz)*exp_zr
            V(Vzpx,cx,cy,cz)=V(Vzpx,cx,cy,cz)*exp_zr
            
            V(Vxny,cx,cy,cz)=V(Vxny,cx,cy,cz)*exp_xr
            V(Vxpy,cx,cy,cz)=V(Vxpy,cx,cy,cz)*exp_xr
            V(Vzny,cx,cy,cz)=V(Vzny,cx,cy,cz)*exp_zr
            V(Vzpy,cx,cy,cz)=V(Vzpy,cx,cy,cz)*exp_zr
                                   
            V(Vxnz,cx,cy,cz)=V(Vxnz,cx,cy,cz)*exp_xr
            V(Vxpz,cx,cy,cz)=V(Vxpz,cx,cy,cz)*exp_xr
            V(Vynz,cx,cy,cz)=V(Vynz,cx,cy,cz)*exp_yr
            V(Vypz,cx,cy,cz)=V(Vypz,cx,cy,cz)*exp_yr
            
!9999 CONTINUE
