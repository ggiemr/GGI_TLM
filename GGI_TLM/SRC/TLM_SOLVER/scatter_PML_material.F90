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
! scatter_PML_m1aterial.F90: File to be included in scatter.F90 to handle the PML cell centre update with non-dispersive dielectric materials
!
!  Started 7/4/2022

 	    output_number=cell_update_code_to_output_number(special_cell_count)

! STAGE 1: get the parameters for this PML cell

	    Ysx=volume_material_list(material_number)%Ys_eps_f/Y0
	    Ysy=volume_material_list(material_number)%Ys_eps_f/Y0
	    Ysz=volume_material_list(material_number)%Ys_eps_f/Y0
            
	    Gex=volume_material_list(material_number)%Ge/Y0
	    Gey=volume_material_list(material_number)%Ge/Y0
	    Gez=volume_material_list(material_number)%Ge/Y0

            Zsx=volume_material_list(material_number)%Zs_mu_f
	    Zsy=volume_material_list(material_number)%Zs_mu_f
	    Zsz=volume_material_list(material_number)%Zs_mu_f
            
	    Rmx=volume_material_list(material_number)%Rm 
	    Rmy=volume_material_list(material_number)%Rm 
	    Rmz=volume_material_list(material_number)%Rm 

            PML_parameter=PML_cell_data(PML_cell)%PML_parameter_array_pos
	    
            PML_material_cell=PML_cell_data(PML_cell)%PML_material_array_pos
	    
	    if(PML_material_cell.EQ.0) then
	      write(*,*)'ERROR finding the PML material cell data'
	      STOP 1
	    end if
   
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
	    
	    GpY2=(Gex+Ysx)/2d0	    
            Ax0=      8d0+4d0*GpY2  +  2d0*dt*(sy+sz)*(1d0+GpY2)  +      dt*dt*sy*sz*GpY2	    
	    Ax1=-2d0*(8d0+4d0*GpY2)                               + 2d0*(dt*dt*sy*sz*GpY2)	    
	    Ax2=      8d0+4d0*GpY2  -  2d0*dt*(sy+sz)*(1d0+GpY2)  +      dt*dt*sy*sz*GpY2
	    
	    GpY2=(Gey+Ysy)/2d0	    
            Ay0=      8d0+4d0*GpY2  +  2d0*dt*(sz+sx)*(1d0+GpY2)  +      dt*dt*sz*sx*GpY2	    
	    Ay1=-2d0*(8d0+4d0*GpY2)                               + 2d0*(dt*dt*sz*sx*GpY2)	    
	    Ay2=      8d0+4d0*GpY2  -  2d0*dt*(sz+sx)*(1d0+GpY2)  +      dt*dt*sz*sx*GpY2
	    
	    GpY2=(Gez+Ysz)/2d0	    
            Az0=      8d0+4d0*GpY2  +  2d0*dt*(sx+sy)*(1d0+GpY2)  +      dt*dt*sx*sy*GpY2	    
	    Az1=-2d0*(8d0+4d0*GpY2)                               + 2d0*(dt*dt*sx*sy*GpY2)	    
	    Az2=      8d0+4d0*GpY2  -  2d0*dt*(sx+sy)*(1d0+GpY2)  +      dt*dt*sx*sy*GpY2
	    
            Bx0= 4d0+2d0*dt*sy
	    Bx1=-8d0
	    Bx2= 4d0-2d0*dt*sy
	    
            By0= 4d0+2d0*dt*sz
	    By1=-8d0
	    By2= 4d0-2d0*dt*sz
	    
            Bz0= 4d0+2d0*dt*sx
	    Bz1=-8d0
	    Bz2= 4d0-2d0*dt*sx

	    
            Cx0= 4d0+2d0*dt*sz
	    Cx1=-8d0
	    Cx2= 4d0-2d0*dt*sz
	    
            Cy0= 4d0+2d0*dt*sx
	    Cy1=-8d0
	    Cy2= 4d0-2d0*dt*sx
	    
            Cz0= 4d0+2d0*dt*sy
	    Cz1=-8d0
	    Cz2= 4d0-2d0*dt*sy
	    
            Dx0= 4d0*Ysx + 2d0*dt*(sy+sz)*Ysx + 1d0*(dt*dt*sy*sz*Ysx)
	    Dx1=-8d0*Ysx                      + 2d0*(dt*dt*sy*sz*Ysx)
	    Dx2= 4d0*Ysx - 2d0*dt*(sy+sz)*Ysx + 1d0*(dt*dt*sy*sz*Ysx)
	    
            Dy0= 4d0*Ysy + 2d0*dt*(sz+sx)*Ysy + 1d0*(dt*dt*sz*sx*Ysy)
	    Dy1=-8d0*Ysy                      + 2d0*(dt*dt*sz*sx*Ysy)
	    Dy2= 4d0*Ysy - 2d0*dt*(sz+sx)*Ysy + 1d0*(dt*dt*sz*sx*Ysy)
	    
            Dz0= 4d0*Ysz + 2d0*dt*(sx+sy)*Ysz + 1d0*(dt*dt*sx*sy*Ysz)
	    Dz1=-8d0*Ysz                      + 2d0*(dt*dt*sx*sy*Ysz)
	    Dz2= 4d0*Ysz - 2d0*dt*(sx+sy)*Ysz + 1d0*(dt*dt*sx*sy*Ysz)

!            write(*,'(A,I3,A,3I4,A,3ES12.3)')'rank=',rank,' cx, cy,cz=',cx,cy,cz,' sx, sy, sz=',sx,sy,sz
	    
	    ZpR=Zsx+Rmx	      
	    Sx0=16d0*Z0+4d0*Zpr  +  2d0*dt*(sy+sz)*(2d0*Z0+Zpr)  +     dt*dt*sy*sz*ZpR      
	    Sx1=-32d0*Z0-8d0*Zpr                           	 + 2d0*dt*dt*sy*sz*ZpR
	    Sx2=16d0*Z0+4d0*Zpr  -  2d0*dt*(sy+sz)*(2d0*Z0+Zpr)  +     dt*dt*sy*sz*ZpR  

	    ZpR=Zsy+Rmy	      
	    Sy0=16d0*Z0+4d0*Zpr  +  2d0*dt*(sz+sx)*(2d0*Z0+Zpr)  +     dt*dt*sz*sx*ZpR      
	    Sy1=-32d0*Z0-8d0*Zpr                           	 + 2d0*dt*dt*sz*sx*ZpR
	    Sy2=16d0*Z0+4d0*Zpr  -  2d0*dt*(sz+sx)*(2d0*Z0+Zpr)  +     dt*dt*sz*sx*ZpR

	    ZpR=Zsz+Rmz	      
	    Sz0=16d0*Z0+4d0*Zpr  +  2d0*dt*(sx+sy)*(2d0*Z0+Zpr)  +     dt*dt*sx*sy*ZpR      
	    Sz1=-32d0*Z0-8d0*Zpr                           	 + 2d0*dt*dt*sx*sy*ZpR
	    Sz2=16d0*Z0+4d0*Zpr  -  2d0*dt*(sx+sy)*(2d0*Z0+Zpr)  +     dt*dt*sx*sy*ZpR  
	      
	    Tx0= 8d0  +  4d0*dt*sx	
	    Tx1=-16d0
	    Tx2= 8d0  -  4d0*dt*sx
	      
	    Ty0= 8d0  +  4d0*dt*sy	
	    Ty1=-16d0
	    Ty2= 8d0  -  4d0*dt*sy
	      
	    Tz0= 8d0  +  4d0*dt*sz	
	    Tz1=-16d0
	    Tz2= 8d0  -  4d0*dt*sz
	    
	    Uxy0= 2d0  +  dt*sz
	    Uxy1=-2d0  +  dt*sz	    
	    Wxy0=Z0*( 2d0  +  dt*sy)
	    Wxy1=Z0*(-2d0  +  dt*sy)
	    
	    Uyx0= 2d0  +  dt*sz
	    Uyx1=-2d0  +  dt*sz	    
	    Wyx0=Z0*( 2d0  +  dt*sx)
	    Wyx1=Z0*(-2d0  +  dt*sx)
	    
	    Uyz0= 2d0  +  dt*sx
	    Uyz1=-2d0  +  dt*sx	    
	    Wyz0=Z0*( 2d0  +  dt*sz)
	    Wyz1=Z0*(-2d0  +  dt*sz)

	    Uzy0= 2d0  +  dt*sx
	    Uzy1=-2d0  +  dt*sx	    
	    Wzy0=Z0*( 2d0  +  dt*sy)
	    Wzy1=Z0*(-2d0  +  dt*sy)
	    
	    Uzx0= 2d0  +  dt*sy
	    Uzx1=-2d0  +  dt*sy  
	    Wzx0=Z0*( 2d0  +  dt*sx)
	    Wzx1=Z0*(-2d0  +  dt*sx)
	    
	    Uxz0= 2d0  +  dt*sy
	    Uxz1=-2d0  +  dt*sy  
	    Wxz0=Z0*( 2d0  +  dt*sz)
	    Wxz1=Z0*(-2d0  +  dt*sz)

! Inductive stub update stuff	    
	      
	    Px0= 4d0  +  2d0*dt*sx	
	    Px1=-8d0
	    Px2= 4d0  -  2d0*dt*sx
	      
	    Py0= 4d0  +  2d0*dt*sy	
	    Py1=-8d0
	    Py2= 4d0  -  2d0*dt*sy
	      
	    Pz0= 4d0  +  2d0*dt*sz	
	    Pz1=-8d0
	    Pz2= 4d0  -  2d0*dt*sz
	    
            Qx0= 4d0*Zsx + 2d0*dt*(sy+sz)*Zsx + 1d0*(dt*dt*sy*sz*Zsx)
	    Qx1=-8d0*Zsx                      + 2d0*(dt*dt*sy*sz*Zsx)
	    Qx2= 4d0*Zsx - 2d0*dt*(sy+sz)*Zsx + 1d0*(dt*dt*sy*sz*Zsx)
	    
            Qy0= 4d0*Zsy + 2d0*dt*(sz+sx)*Zsy + 1d0*(dt*dt*sz*sx*Zsy)
	    Qy1=-8d0*Zsy                      + 2d0*(dt*dt*sz*sx*Zsy)
	    Qy2= 4d0*Zsy - 2d0*dt*(sz+sx)*Zsy + 1d0*(dt*dt*sz*sx*Zsy)
	    
            Qz0= 4d0*Zsz + 2d0*dt*(sx+sy)*Zsz + 1d0*(dt*dt*sx*sy*Zsz)
	    Qz1=-8d0*Zsz                      + 2d0*(dt*dt*sx*sy*Zsz)
	    Qz2= 4d0*Zsz - 2d0*dt*(sx+sy)*Zsz + 1d0*(dt*dt*sx*sy*Zsz)	    
                       
! STAGE 2: Retrieve voltages and currents from the last two timesteps

! Data from PML_cell
            Ix_m1=PML_cell_data(PML_cell)%Ix
            Iy_m1=PML_cell_data(PML_cell)%Iy
            Iz_m1=PML_cell_data(PML_cell)%Iz
   
            Ixt_m1=PML_cell_data(PML_cell)%Ixt    
            Iyt_m1=PML_cell_data(PML_cell)%Iyt
            Izt_m1=PML_cell_data(PML_cell)%Izt

            Vi_m1(1:12)=PML_cell_data(PML_cell)%Vi(1:12)
  
            Vxt_m1=PML_cell_data(PML_cell)%Vxt
            Vyt_m1=PML_cell_data(PML_cell)%Vyt
            Vzt_m1=PML_cell_data(PML_cell)%Vzt	    
            
            Vyxt_m1=PML_cell_data(PML_cell)%Vyxt
            Vzxt_m1=PML_cell_data(PML_cell)%Vzxt
            Vxyt_m1=PML_cell_data(PML_cell)%Vxyt
            Vzyt_m1=PML_cell_data(PML_cell)%Vzyt
            Vyzt_m1=PML_cell_data(PML_cell)%Vyzt
            Vxzt_m1=PML_cell_data(PML_cell)%Vxzt

! Data from PML_material_cell

            Ix_m2=PML_material_cell_data(PML_material_cell)%Ix_m1
            Iy_m2=PML_material_cell_data(PML_material_cell)%Iy_m1
            Iz_m2=PML_material_cell_data(PML_material_cell)%Iz_m1
   
            Ixt_m2=PML_material_cell_data(PML_material_cell)%Ixt_m1  
            Iyt_m2=PML_material_cell_data(PML_material_cell)%Iyt_m1
            Izt_m2=PML_material_cell_data(PML_material_cell)%Izt_m1

            Vcx_m2=PML_material_cell_data(PML_material_cell)%Vcx_m1
            Vcy_m2=PML_material_cell_data(PML_material_cell)%Vcy_m1
            Vcz_m2=PML_material_cell_data(PML_material_cell)%Vcz_m1

            Vcx_m1=PML_material_cell_data(PML_material_cell)%Vcx
            Vcy_m1=PML_material_cell_data(PML_material_cell)%Vcy
            Vcz_m1=PML_material_cell_data(PML_material_cell)%Vcz
           
            Vi_m2(1:12)=PML_material_cell_data(PML_material_cell)%Vi_m1(1:12)
  
            Vxt_m2=PML_material_cell_data(PML_material_cell)%Vxt_m1
            Vyt_m2=PML_material_cell_data(PML_material_cell)%Vyt_m1
            Vzt_m2=PML_material_cell_data(PML_material_cell)%Vzt_m1
            
            Vyxt_m2=PML_material_cell_data(PML_material_cell)%Vyxt_m1
            Vzxt_m2=PML_material_cell_data(PML_material_cell)%Vzxt_m1
            Vxyt_m2=PML_material_cell_data(PML_material_cell)%Vxyt_m1
            Vzyt_m2=PML_material_cell_data(PML_material_cell)%Vzyt_m1
            Vyzt_m2=PML_material_cell_data(PML_material_cell)%Vyzt_m1
            Vxzt_m2=PML_material_cell_data(PML_material_cell)%Vxzt_m1

! capacitive stub: assumes non-dispersive capacitive stub. Vosx is the incident voltage at the current timestep = Vr at previous timestep
	   
            Vosx=PML_material_cell_data(PML_material_cell)%Vosxr    
            Vosy=PML_material_cell_data(PML_material_cell)%Vosyr    
            Vosz=PML_material_cell_data(PML_material_cell)%Voszr    

! get previous capacitive stub incident voltages
	    
	    Vosx_m1=PML_material_cell_data(PML_material_cell)%Vosx   
	    Vosx_m2=PML_material_cell_data(PML_material_cell)%Vosx_m1

	    Vosy_m1=PML_material_cell_data(PML_material_cell)%Vosy   
	    Vosy_m2=PML_material_cell_data(PML_material_cell)%Vosy_m1

	    Vosz_m1=PML_material_cell_data(PML_material_cell)%Vosz   
	    Vosz_m2=PML_material_cell_data(PML_material_cell)%Vosz_m1

! inductive stub: assumes non-dispersive inductive stub. Vssx is the incident voltage at the current timestep = -Vr at previous timestep
	   
            Vssx=-PML_material_cell_data(PML_material_cell)%Vssxr
            Vssy=-PML_material_cell_data(PML_material_cell)%Vssyr    
            Vssz=-PML_material_cell_data(PML_material_cell)%Vsszr    

! get previous inductive stub impedance voltages
	    
	    Vsx_m1=PML_material_cell_data(PML_material_cell)%Vsx   
	    Vsx_m2=PML_material_cell_data(PML_material_cell)%Vsx_m1

	    Vsy_m1=PML_material_cell_data(PML_material_cell)%Vsy   
	    Vsy_m2=PML_material_cell_data(PML_material_cell)%Vsy_m1

	    Vsz_m1=PML_material_cell_data(PML_material_cell)%Vsz   
	    Vsz_m2=PML_material_cell_data(PML_material_cell)%Vsz_m1

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


! STAGE 4: Stretched Vxt, Vyt, Vzt:
	    
            Vxt=(1d0/Ax0)*(                                                                                            &
	      Bx0*( V(Vznx,cx,cy,cz)+V(Vzpx,cx,cy,cz) ) + Cx0*( V(Vynx,cx,cy,cz)+V(Vypx,cx,cy,cz) ) + Dx0*( Vosx )     &
	     +Bx1*( Vi_m1(Vznx) +Vi_m1(Vzpx) )	        + Cx1*( Vi_m1(Vynx) +Vi_m1(Vypx) )	    + Dx1*( Vosx_m1 )  &
	     +Bx2*( Vi_m2(Vznx) +Vi_m2(Vzpx) )	        + Cx2*( Vi_m2(Vynx) +Vi_m2(Vypx) )	    + Dx2*( Vosx_m2 )  &
	     -Ax1*Vxt_m1 - Ax2*Vxt_m2)

            Vyt=(1d0/Ay0)*(                                                                                            &
	      By0*( V(Vxny,cx,cy,cz)+V(Vxpy,cx,cy,cz) ) + Cy0*( V(Vzny,cx,cy,cz)+V(Vzpy,cx,cy,cz) ) + Dy0*( Vosy )     &
	     +By1*( Vi_m1(Vxny) +Vi_m1(Vxpy) )  	+ Cy1*( Vi_m1(Vzny) +Vi_m1(Vzpy) )	    + Dy1*( Vosy_m1 )  &
	     +By2*( Vi_m2(Vxny) +Vi_m2(Vxpy) )  	+ Cy2*( Vi_m2(Vzny) +Vi_m2(Vzpy) )	    + Dy2*( Vosy_m2 )  &
	     -Ay1*Vyt_m1 - Ay2*Vyt_m2)

            Vzt=(1d0/Az0)*(                                                                                            &
	      Bz0*( V(Vynz,cx,cy,cz)+V(Vypz,cx,cy,cz) ) + Cz0*( V(Vxnz,cx,cy,cz)+V(Vxpz,cx,cy,cz) ) + Dz0*( Vosz )     &
	     +Bz1*( Vi_m1(Vynz) +Vi_m1(Vypz) )  	+ Cz1*( Vi_m1(Vxnz) +Vi_m1(Vxpz) )	    + Dz1*( Vosz_m1 )  &
	     +Bz2*( Vi_m2(Vynz) +Vi_m2(Vypz) )  	+ Cz2*( Vi_m2(Vxnz) +Vi_m2(Vxpz) )	    + Dz2*( Vosz_m2 )  &
	     -Az1*Vzt_m1 - Az2*Vzt_m2)    

! STAGE 5: Stretched Ixt, Iyt, Izt calculation
	    
	    Vcx=V(Vypz,cx,cy,cz)-V(Vzpy,cx,cy,cz)-V(Vynz,cx,cy,cz)+V(Vzny,cx,cy,cz) - Vssx
	    Vcy=V(Vzpx,cx,cy,cz)-V(Vxpz,cx,cy,cz)-V(Vznx,cx,cy,cz)+V(Vxnz,cx,cy,cz) - Vssy
	    Vcz=V(Vxpy,cx,cy,cz)-V(Vypx,cx,cy,cz)-V(Vxny,cx,cy,cz)+V(Vynx,cx,cy,cz) - Vssz
	    	    
            Ixt=(1d0/Sx0)*( Tx0*Vcx + Tx1*Vcx_m1 + Tx2*Vcx_m2 -Sx1*Ixt_m1 - Sx2*Ixt_m2)
            Iyt=(1d0/Sy0)*( Ty0*Vcy + Ty1*Vcy_m1 + Ty2*Vcy_m2 -Sy1*Iyt_m1 - Sy2*Iyt_m2)
            Izt=(1d0/Sz0)*( Tz0*Vcz + Tz1*Vcz_m1 + Tz2*Vcz_m2 -Sz1*Izt_m1 - Sz2*Izt_m2)
                                        
! i=x j=y k=z
            Vxyt=(1d0/Uxy0)*( Wxy0*Izt + Wxy1*Izt_m1 - Uxy1*Vxyt_m1 )
            Vyxt=(1d0/Uyx0)*( Wyx0*Izt + Wyx1*Izt_m1 - Uyx1*Vyxt_m1 )
	    
            Vyzt=(1d0/Uyz0)*( Wyz0*Ixt + Wyz1*Ixt_m1 - Uyz1*Vyzt_m1 )
            Vzyt=(1d0/Uzy0)*( Wzy0*Ixt + Wzy1*Ixt_m1 - Uzy1*Vzyt_m1 )
	    
            Vxzt=(1d0/Uxz0)*( Wxz0*Iyt + Wxz1*Iyt_m1 - Uxz1*Vxzt_m1 )
            Vzxt=(1d0/Uzx0)*( Wzx0*Iyt + Wzx1*Iyt_m1 - Uzx1*Vzxt_m1 )
	    
! STAGE 9:  Save voltages and currents for the next timestep

            PML_material_cell_data(PML_material_cell)%Ix_m1=PML_cell_data(PML_cell)%Ix
            PML_material_cell_data(PML_material_cell)%Iy_m1=PML_cell_data(PML_cell)%Iy
            PML_material_cell_data(PML_material_cell)%Iz_m1=PML_cell_data(PML_cell)%Iz

            PML_cell_data(PML_cell)%Ix=Ix
            PML_cell_data(PML_cell)%Iy=Iy
            PML_cell_data(PML_cell)%Iz=Iz
            
            PML_material_cell_data(PML_material_cell)%Vi_m1(1:12)=PML_cell_data(PML_cell)%Vi(1:12)            
            PML_cell_data(PML_cell)%Vi(1:12)=V(1:12,cx,cy,cz)
  
            PML_material_cell_data(PML_material_cell)%Vxt_m1=PML_cell_data(PML_cell)%Vxt
            PML_material_cell_data(PML_material_cell)%Vyt_m1=PML_cell_data(PML_cell)%Vyt
            PML_material_cell_data(PML_material_cell)%Vzt_m1=PML_cell_data(PML_cell)%Vzt
  
            PML_cell_data(PML_cell)%Vxt=Vxt
            PML_cell_data(PML_cell)%Vyt=Vyt
            PML_cell_data(PML_cell)%Vzt=Vzt

            PML_material_cell_data(PML_material_cell)%Ixt_m1=PML_cell_data(PML_cell)%Ixt
            PML_material_cell_data(PML_material_cell)%Iyt_m1=PML_cell_data(PML_cell)%Iyt
            PML_material_cell_data(PML_material_cell)%Izt_m1=PML_cell_data(PML_cell)%Izt

            PML_cell_data(PML_cell)%Ixt=Ixt
            PML_cell_data(PML_cell)%Iyt=Iyt
            PML_cell_data(PML_cell)%Izt=Izt

            PML_material_cell_data(PML_material_cell)%Vcx_m1=PML_material_cell_data(PML_material_cell)%Vcx
            PML_material_cell_data(PML_material_cell)%Vcy_m1=PML_material_cell_data(PML_material_cell)%Vcy
            PML_material_cell_data(PML_material_cell)%Vcz_m1=PML_material_cell_data(PML_material_cell)%Vcz

            PML_material_cell_data(PML_material_cell)%Vcx=Vcx
            PML_material_cell_data(PML_material_cell)%Vcy=Vcy
            PML_material_cell_data(PML_material_cell)%Vcz=Vcz
              
            PML_material_cell_data(PML_material_cell)%Vxyt_m1=PML_cell_data(PML_cell)%Vxyt
            PML_material_cell_data(PML_material_cell)%Vyxt_m1=PML_cell_data(PML_cell)%Vyxt
	    
            PML_cell_data(PML_cell)%Vxyt=Vxyt
            PML_cell_data(PML_cell)%Vyxt=Vyxt
            
	    
            PML_material_cell_data(PML_material_cell)%Vzxt_m1=PML_cell_data(PML_cell)%Vzxt
            PML_material_cell_data(PML_material_cell)%Vxzt_m1=PML_cell_data(PML_cell)%Vxzt
	    
            PML_cell_data(PML_cell)%Vzxt=Vzxt
            PML_cell_data(PML_cell)%Vxzt=Vxzt
     
            
            PML_material_cell_data(PML_material_cell)%Vzyt_m1=PML_cell_data(PML_cell)%Vzyt
            PML_material_cell_data(PML_material_cell)%Vyzt_m1=PML_cell_data(PML_cell)%Vyzt
            
            PML_cell_data(PML_cell)%Vzyt=Vzyt
            PML_cell_data(PML_cell)%Vyzt=Vyzt

! STAGE 11: Calculate scattered voltages 

            V_max=V(Vypx,cx,cy,cz)
            V_min=V(Vynx,cx,cy,cz)    
	    V(Vynx,cx,cy,cz)  =Vxt-Vyxt-V_max
	    V(Vypx,cx,cy,cz)  =Vxt+Vyxt-V_min
	    
            V_max=V(Vzpx,cx,cy,cz)
            V_min=V(Vznx,cx,cy,cz) 
	    V(Vznx,cx,cy,cz)=  Vxt+Vzxt-V_max
	    V(Vzpx,cx,cy,cz)=  Vxt-Vzxt-V_min
	
            V_max=V(Vzpy,cx,cy,cz)
            V_min=V(Vzny,cx,cy,cz) 
	    V(Vzny,cx,cy,cz)=  Vyt-Vzyt-V_max
	    V(Vzpy,cx,cy,cz)=  Vyt+Vzyt-V_min
	
            V_max=V(Vxpy,cx,cy,cz)
            V_min=V(Vxny,cx,cy,cz) 
	    V(Vxny,cx,cy,cz)=  Vyt+Vxyt-V_max
	    V(Vxpy,cx,cy,cz)=  Vyt-Vxyt-V_min
	
            V_max=V(Vxpz,cx,cy,cz)
            V_min=V(Vxnz,cx,cy,cz) 
	    V(Vxnz,cx,cy,cz)=  Vzt-Vxzt-V_max
	    V(Vxpz,cx,cy,cz)=  Vzt+Vxzt-V_min
	
            V_max=V(Vypz,cx,cy,cz)
            V_min=V(Vynz,cx,cy,cz) 
	    V(Vynz,cx,cy,cz)=  Vzt+Vyzt-V_max
  	    V(Vypz,cx,cy,cz)=  Vzt-Vyzt-V_min	  

! Scatter the capacitive stub voltage

            PML_material_cell_data(PML_material_cell)%Vosxr=Vxt-Vosx
            PML_material_cell_data(PML_material_cell)%Vosyr=Vyt-Vosy
            PML_material_cell_data(PML_material_cell)%Voszr=Vzt-Vosz

! Save the capacitive stub voltage data	    
            PML_material_cell_data(PML_material_cell)%Vosx_m1=PML_material_cell_data(PML_material_cell)%Vosx
            PML_material_cell_data(PML_material_cell)%Vosx=Vosx
	    
            PML_material_cell_data(PML_material_cell)%Vosy_m1=PML_material_cell_data(PML_material_cell)%Vosy
            PML_material_cell_data(PML_material_cell)%Vosy=Vosy
	    
            PML_material_cell_data(PML_material_cell)%Vosz_m1=PML_material_cell_data(PML_material_cell)%Vosz
            PML_material_cell_data(PML_material_cell)%Vosz=Vosz

! Scatter the inductive stub Voltage across the stub impedance only in the Thevenin equivalent cirucuit
            Vsx=(1d0/Px0)*( Qx0*Ixt + Qx1*Ixt_m1 + Qx2*Ixt_m2 -Px1*Vsx_m1 - Px2*Vsx_m2)     
	    Vsy=(1d0/Py0)*( Qy0*Iyt + Qy1*Iyt_m1 + Qy2*Iyt_m2 -Py1*Vsy_m1 - Py2*Vsy_m2)
	    Vsz=(1d0/Pz0)*( Qz0*Izt + Qz1*Izt_m1 + Qz2*Izt_m2 -Pz1*Vsz_m1 - Pz2*Vsz_m2)	    

! Note that (Vcx+Vssx) is effectively the contribution from the link line incident voltages only
	    Vsx_temp=2d0*(Vcx+Vssx)-2d0*(Vyzt+Vzyt) - 2d0*Vssx
	    Vsy_temp=2d0*(Vcy+Vssy)-2d0*(Vxzt+Vzxt) - 2d0*Vssy
	    Vsz_temp=2d0*(Vcz+Vssz)-2d0*(Vxyt+Vyxt) - 2d0*Vssz
	    
            PML_material_cell_data(PML_material_cell)%Vssxr=(Vsx+Vssx)
            PML_material_cell_data(PML_material_cell)%Vssyr=(Vsy+Vssy)
            PML_material_cell_data(PML_material_cell)%Vsszr=(Vsz+Vssz)

! Save the inductive stub impedance voltage data	    
            PML_material_cell_data(PML_material_cell)%Vsx_m1=PML_material_cell_data(PML_material_cell)%Vsx
            PML_material_cell_data(PML_material_cell)%Vsx=Vsx
	    
            PML_material_cell_data(PML_material_cell)%Vsy_m1=PML_material_cell_data(PML_material_cell)%Vsy
            PML_material_cell_data(PML_material_cell)%Vsy=Vsy
	    
            PML_material_cell_data(PML_material_cell)%Vsz_m1=PML_material_cell_data(PML_material_cell)%Vsz
            PML_material_cell_data(PML_material_cell)%Vsz=Vsz
         	  
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
            
