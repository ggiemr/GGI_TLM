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
! SUBROUTINE set_pml_volume_material_mesh
! SUBROUTINE calc_pml_dist
! SUBROUTINE calc_sigma_0
! FUNCTION PML_conductivity
!
! NAME
!     set_pml_volume_material_mesh
!
! DESCRIPTION
!     
!     Count the number of cells of each material and construct the 
!     appropriate lists for the main solver
!     
! COMMENTS
!     Allocate a local mesh array and fill with PML material data. This ensures that
!     a cell doesn't get given more than one material property.
!     In overlapping volumes the last material assignment takes precidence
!
! HISTORY
!
!     started 1/10/2019 CJS
!
!
SUBROUTINE set_pml_volume_material_mesh

USE Constants
USE TLM_general
USE geometry
USE PML_module
USE mesh
USE cell_parameters
USE file_information

IMPLICIT NONE

! local variables

  integer	:: volume_number
  integer	:: cell,number_of_cells
  integer	:: total_number_of_PML_cells_all_processes
  integer	:: n_data
  integer	:: cx,cy,cz
  integer	:: dist_xmin,dist_xmax,dist_ymin,dist_ymax,dist_zmin,dist_zmax
  integer	:: i
  integer	:: cell_number
  integer       :: pml_face
  integer       :: pml_cell
  integer       :: pml_x,pml_y,pml_z
  
  real*8        :: sigma_x,sigma_y,sigma_z
  
  logical :: write_pml_data_to_info_file

! function variables  

  real*8 :: PML_conductivity
  
! START
  
  CALL write_line('CALLED: set_pml_volume_material_mesh',0,output_to_screen_flag)

! Flag to write PML data for every PML cell. Used for testing but produces a lot of data in the info file  
!  write_pml_data_to_info_file=.TRUE.
  write_pml_data_to_info_file=.FALSE.

! INITIALISE VOLUME MATERIALS  
    
!  CALL write_line_integer('n_pml_volumes=',n_pml_volumes,0,output_to_screen_flag)
  
  write(*,*)'rank=',rank,' n_pml_volumes=',n_pml_volumes

  if (rank.eq.0) then
  
    write(info_file_unit,*)'____________________________________________________'
    write(info_file_unit,*)''
    write(info_file_unit,*)'#START OF PML DESCRIPTION'
    write(info_file_unit,*)''
    write(info_file_unit,*)'Number of PML volumes=',n_pml_volumes
    write(info_file_unit,*)''
  
  end if
  
! loop over the PML surfaces and flag the distance into the PML in x,y and z directions
! also check for intersection with other volume materials
  
  local_cell_PML(:,:,:,:)=0
  total_number_of_PML_cells=0

  do volume_number=1,n_pml_volumes
  
    pml_face=pml_volume_to_face(volume_number)

    number_of_cells=pml_volumes(volume_number)%number_of_cells

!    CALL write_line_integer('PML volume     =',volume_number  ,0,output_to_screen_flag)
!    CALL write_line_integer('Number of cells=',number_of_cells,0,output_to_screen_flag)

    write(*,*)'rank=',rank,' PML volume     =',volume_number 
    write(*,*)'rank=',rank,' Number of cells=',number_of_cells
  
    if (number_of_cells.gt.0) then
      
      total_number_of_PML_cells=total_number_of_PML_cells+number_of_cells
      
      do cell=1,number_of_cells
    
	cx=pml_volumes(volume_number)%cell_list(cell)%cell%i
	cy=pml_volumes(volume_number)%cell_list(cell)%cell%j
	cz=pml_volumes(volume_number)%cell_list(cell)%cell%k
        
        CALL calc_pml_dist(dist_xmin,dist_xmax,cx,nx,pml_txmin,pml_txmax,dl)
        CALL calc_pml_dist(dist_ymin,dist_ymax,cy,ny,pml_tymin,pml_tymax,dl)
        CALL calc_pml_dist(dist_zmin,dist_zmax,cz,nz,pml_tzmin,pml_tzmax,dl)
        
	if (rank.eq.cell_rank(cz)) then
        
          local_cell_PML(cx,cy,cz,0)=1
	  
          if (pml_face.EQ.pml_face_xmin) then
          
            local_cell_PML(cx,cy,cz,1)=dist_xmin
            
          else if (pml_face.EQ.pml_face_xmax) then
          
            local_cell_PML(cx,cy,cz,1)=dist_xmax
            
          else if (pml_face.EQ.pml_face_ymin) then
          
            local_cell_PML(cx,cy,cz,2)=dist_ymin
            
          else if (pml_face.EQ.pml_face_ymax) then
          
            local_cell_PML(cx,cy,cz,2)=dist_ymax
            
          else if (pml_face.EQ.pml_face_zmin) then
          
            local_cell_PML(cx,cy,cz,3)=dist_zmin
            
          else if (pml_face.EQ.pml_face_zmax) then
          
            local_cell_PML(cx,cy,cz,3)=dist_zmax
            
	  end if
          
	end if
	
      end do ! next cell
    
    end if ! number_of_cells.gt.0

  end do ! next volume number
  
!  CALL write_line_integer('Total number of cells (not counting overlaps)=',total_number_of_PML_cells,0,output_to_screen_flag)

  write(*,*)'rank=',rank,' Total number of cells (not counting overlaps)=',total_number_of_PML_cells
  
! count the number of cells, also the extent of the PML distances in x, y and z
! Here we can also give each PML cell a unique number which points into the PML_cell_data array

  pml_xmin=0
  pml_xmax=0
  pml_ymin=0
  pml_ymax=0
  pml_zmin=0
  pml_zmax=0

  total_number_of_PML_cells=0
  do cz=nz1,nz2
    do cy=1,ny
      do cx=1,nx
        if ( local_cell_PML(cx,cy,cz,0).NE.0 ) then
        
          total_number_of_PML_cells=total_number_of_PML_cells+1
          
          local_cell_PML(cx,cy,cz,0)=total_number_of_PML_cells
          
          pml_xmin=min( pml_xmin,local_cell_PML(cx,cy,cz,1) )
          pml_xmax=max( pml_xmax,local_cell_PML(cx,cy,cz,1) )
          pml_ymin=min( pml_ymin,local_cell_PML(cx,cy,cz,2) )
          pml_ymax=max( pml_ymax,local_cell_PML(cx,cy,cz,2) )
          pml_zmin=min( pml_zmin,local_cell_PML(cx,cy,cz,3) )
          pml_zmax=max( pml_zmax,local_cell_PML(cx,cy,cz,3) )
          
        end if
      end do
    end do
  end do
  
  npml_x=pml_xmax-pml_xmin+1
  npml_y=pml_ymax-pml_ymin+1
  npml_z=pml_zmax-pml_zmin+1  
 
!  CALL write_line_integer('Total number of PML cells =',total_number_of_PML_cells,0,output_to_screen_flag)

  write(*,*)'rank=',rank,' Total number of PML cells =',total_number_of_PML_cells

  ALLOCATE( PML_cell_data(1:total_number_of_PML_cells) )
  
! Reset PML_cell_data

  do i=1,total_number_of_PML_cells

    PML_cell_data(i)%PML_parameter_array_pos=0

    PML_cell_data(i)%PML_material_array_pos=0
  
    PML_cell_data(i)%Vxt=0d0
    PML_cell_data(i)%Vyt=0d0
    PML_cell_data(i)%Vzt=0d0
    PML_cell_data(i)%Ix=0d0
    PML_cell_data(i)%Iy=0d0
    PML_cell_data(i)%Iz=0d0
    PML_cell_data(i)%Ixt=0d0
    PML_cell_data(i)%Iyt=0d0
    PML_cell_data(i)%Izt=0d0
    PML_cell_data(i)%Vi(1:12)=0d0
    PML_cell_data(i)%Vyxt=0d0
    PML_cell_data(i)%Vzxt=0d0
    PML_cell_data(i)%Vxyt=0d0
    PML_cell_data(i)%Vzyt=0d0
    PML_cell_data(i)%Vyzt=0d0
    PML_cell_data(i)%Vxzt=0d0
     
  end do
  
  CALL write_line_integer('PML xmin=',pml_xmin,0,output_to_screen_flag)
  CALL write_line_integer('PML xmax=',pml_xmax,0,output_to_screen_flag)
  CALL write_line_integer('PML nx  =',npml_x  ,0,output_to_screen_flag)
  
  CALL write_line_integer('PML ymin=',pml_ymin,0,output_to_screen_flag)
  CALL write_line_integer('PML ymax=',pml_ymax,0,output_to_screen_flag)
  CALL write_line_integer('PML ny  =',npml_y  ,0,output_to_screen_flag)
  
  CALL write_line_integer('PML zmin=',pml_zmin,0,output_to_screen_flag)
  CALL write_line_integer('PML zmax=',pml_zmax,0,output_to_screen_flag)
  CALL write_line_integer('PML nz  =',npml_z  ,0,output_to_screen_flag)

! work out the PML conductivity profile parameters, here the electric conductivity of the first cell in the PML

  CALL calc_sigma_0(pml_s0_xmin,pml_r_xmin,pml_xmin,pml_order,dl,eps0,c0)
  CALL calc_sigma_0(pml_s0_xmax,pml_r_xmax,pml_xmax,pml_order,dl,eps0,c0)
  CALL calc_sigma_0(pml_s0_ymin,pml_r_ymin,pml_ymin,pml_order,dl,eps0,c0)
  CALL calc_sigma_0(pml_s0_ymax,pml_r_ymax,pml_ymax,pml_order,dl,eps0,c0)
  CALL calc_sigma_0(pml_s0_zmin,pml_r_zmin,pml_zmin,pml_order,dl,eps0,c0)
  CALL calc_sigma_0(pml_s0_zmax,pml_r_zmax,pml_zmax,pml_order,dl,eps0,c0)
  
! Calculate the normalised conductivity

  pml_s0_xmin=pml_s0_xmin/eps0
  pml_s0_xmax=pml_s0_xmax/eps0
  pml_s0_ymin=pml_s0_ymin/eps0
  pml_s0_ymax=pml_s0_ymax/eps0
  pml_s0_zmin=pml_s0_zmin/eps0
  pml_s0_zmax=pml_s0_zmax/eps0

! Work out the PML parameters for all the different PML cells in the mesh

  PML_array_size=npml_x*npml_y*npml_z
   
  ALLOCATE( PML_dxdydz_to_parameter_array(pml_xmin:pml_xmax,pml_ymin:pml_ymax,pml_zmin:pml_zmax) )
  PML_dxdydz_to_parameter_array(pml_xmin:pml_xmax,pml_ymin:pml_ymax,pml_zmin:pml_zmax)=0

  PML_n_parameters=0
  
  do cz=nz1,nz2
    do cy=1,ny
      do cx=1,nx
        if ( local_cell_PML(cx,cy,cz,0).NE.0 ) then
           
          pml_x=local_cell_PML(cx,cy,cz,1)
          pml_y=local_cell_PML(cx,cy,cz,2)
          pml_z=local_cell_PML(cx,cy,cz,3)
          
          if ( PML_dxdydz_to_parameter_array(pml_x,pml_y,pml_z).EQ.0 ) then ! add a new parameter to the list
            PML_n_parameters=PML_n_parameters+1           
            PML_dxdydz_to_parameter_array(pml_x,pml_y,pml_z)=PML_n_parameters
          end if
           
        end if
      end do
    end do
  end do
  
  write(*,*)'rank=',rank,' PML_n_parameters',PML_n_parameters
  
  ALLOCATE( PML_parameters(1:PML_n_parameters) )

! Fill the PML_parameters array

  PML_n_parameters=0
  
  do cz=nz1,nz2
    do cy=1,ny
      do cx=1,nx
      
        if ( local_cell_PML(cx,cy,cz,0).NE.0 ) then
           
          pml_x=local_cell_PML(cx,cy,cz,1)
          pml_y=local_cell_PML(cx,cy,cz,2)
          pml_z=local_cell_PML(cx,cy,cz,3)
                    
          PML_n_parameters=PML_dxdydz_to_parameter_array(pml_x,pml_y,pml_z)
          
          PML_parameters(PML_n_parameters)%d_x=pml_x
          PML_parameters(PML_n_parameters)%d_y=pml_y
          PML_parameters(PML_n_parameters)%d_z=pml_z
          
          sigma_x=PML_conductivity(pml_x,pml_xmin,pml_xmax,pml_s0_xmin,pml_s0_xmax,pml_order) 
          sigma_y=PML_conductivity(pml_y,pml_ymin,pml_ymax,pml_s0_ymin,pml_s0_ymax,pml_order) 
          sigma_z=PML_conductivity(pml_z,pml_zmin,pml_zmax,pml_s0_zmin,pml_s0_zmax,pml_order) 
          
          PML_parameters(PML_n_parameters)%sx=sigma_x
          PML_parameters(PML_n_parameters)%sy=sigma_y
          PML_parameters(PML_n_parameters)%sz=sigma_z
          
          PML_parameters(PML_n_parameters)%ax=4d0/(4d0+(sigma_y+sigma_z)*dt)
          PML_parameters(PML_n_parameters)%ay=4d0/(4d0+(sigma_z+sigma_x)*dt)
          PML_parameters(PML_n_parameters)%az=4d0/(4d0+(sigma_x+sigma_y)*dt)   
                            
          PML_parameters(PML_n_parameters)%exp_x=exp(-dt*sigma_x/2d0)
          PML_parameters(PML_n_parameters)%exp_y=exp(-dt*sigma_y/2d0)
          PML_parameters(PML_n_parameters)%exp_z=exp(-dt*sigma_z/2d0)
          
          PML_parameters(PML_n_parameters)%Csx=sigma_x*dt/4d0
          PML_parameters(PML_n_parameters)%Csy=sigma_y*dt/4d0
          PML_parameters(PML_n_parameters)%Csz=sigma_z*dt/4d0 
          
          PML_parameters(PML_n_parameters)%Csx2=sigma_x*dt/2d0
          PML_parameters(PML_n_parameters)%Csy2=sigma_y*dt/2d0
          PML_parameters(PML_n_parameters)%Csz2=sigma_z*dt/2d0
          
          PML_parameters(PML_n_parameters)%Csy_sz=1d0-(sigma_y+sigma_z)*dt/4d0
          PML_parameters(PML_n_parameters)%Csz_sx=1d0-(sigma_z+sigma_x)*dt/4d0
          PML_parameters(PML_n_parameters)%Csx_sy=1d0-(sigma_x+sigma_y)*dt/4d0          
           
        end if
        
      end do
    end do
  end do
    
! Link the PML_cell_data array to the PML_parameters array

  if (write_pml_data_to_info_file) then
  
    write(info_file_unit,*)''
    write(info_file_unit,*)' rank ncells cx  cy  cz  p_x p_y p_z  cell  d_x d_y d_z     sx          sy          sz'

  end if
  
  do cz=nz1,nz2
    do cy=1,ny
      do cx=1,nx
      
        if ( local_cell_PML(cx,cy,cz,0).NE.0 ) then
                  
          PML_cell=local_cell_PML(cx,cy,cz,0)

          pml_x=local_cell_PML(cx,cy,cz,1)
          pml_y=local_cell_PML(cx,cy,cz,2)
          pml_z=local_cell_PML(cx,cy,cz,3)
                              
          PML_cell_data(PML_cell)%PML_parameter_array_pos=PML_dxdydz_to_parameter_array(pml_x,pml_y,pml_z)

          if (write_pml_data_to_info_file) then
  
          i=PML_cell_data(PML_cell)%PML_parameter_array_pos
          
            write(info_file_unit,8888)rank,total_number_of_PML_cells,cx,cy,cz,pml_x,pml_y,pml_z,i, &
                                      PML_parameters(i)%d_x,PML_parameters(i)%d_y,PML_parameters(i)%d_z, &
                                      PML_parameters(i)%sx*eps0,PML_parameters(i)%sy*eps0,PML_parameters(i)%sz*eps0
                                 
8888        format(I5,I7,6I4,I6,3I4,3ES12.3)

          end if
          
        end if
      end do
    end do
  end do
   
  if (allocated( PML_dxdydz_to_parameter_array )) DEALLOCATE ( PML_dxdydz_to_parameter_array ) 

#if defined(MPI)

  n_data=1
  call MPI_REDUCE(total_number_of_PML_cells,total_number_of_PML_cells_all_processes , &
                    n_data,MPI_INTEGER, MPI_SUM, 0,MPI_COMM_WORLD,ierror)
		    
#elif defined(SEQ)

  total_number_of_PML_cells_all_processes=total_number_of_PML_cells
    
#endif			

  if (rank.eq.0) then
    CALL write_line_integer('Total number of PML cells=',	&
                            total_number_of_PML_cells_all_processes,0,output_to_screen_flag)
    write(info_file_unit,'(A,I14)')'Total number of PML cells=',total_number_of_PML_cells_all_processes
  end if
 
  if (rank.eq.0) then
    write(info_file_unit,*)''
    write(info_file_unit,*)'#END OF PML MATERIAL DESCRIPTION'
  end if
  
  CALL write_line('FINISHED: set_pml_volume_material_mesh',0,output_to_screen_flag)
    
  RETURN

END SUBROUTINE set_pml_volume_material_mesh
!
! NAME
!     
!
! DESCRIPTION
!     
!     work out the distance into the PML given the cell number, mesh extent and PML thickness
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 1/10/2019 CJS
!
!

SUBROUTINE calc_pml_dist(dist_xmin,dist_xmax,cx,nx,pml_txmin,pml_txmax,dl)

IMPLICIT NONE

  integer	:: dist_xmin,dist_xmax
  integer	:: cx
  integer	:: nx
  real*8        :: pml_txmin,pml_txmax
  real*8        :: dl

! local variables

  integer	:: d_xmin,d_xmax
  integer	:: t_xmin,t_xmax
  
! START

! Thickness of PML in cells on xmin and xmax boundaries

  t_xmin=NINT(pml_txmin/dl)
  t_xmax=NINT(pml_txmax/dl)

! distance from the xmin and xmax boundaries in cells

  d_xmin=cx 
  d_xmax=nx-cx+1 
  
! distance into PML at xmin and xmax boundaries

  dist_xmin=d_xmin-(t_xmin+1)    ! make the distance on min boundaries negative
  if (dist_xmin.GT.0) dist_xmin=0
  
  dist_xmax=t_xmax-d_xmax+1
  if (dist_xmax.LT.0) dist_xmax=0  

END SUBROUTINE calc_pml_dist
!
! NAME
!     SUBROUTINE calc_sigma_0
!
! DESCRIPTION
!     
!     work out the conductivity of the first PML cell from which the other cell conductivites are derived
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 4/10/2019 CJS
!

SUBROUTINE calc_sigma_0(pml_s0,pml_r,N,o,dl,eps0,c0)

IMPLICIT NONE

  real*8        :: pml_s0,pml_r
  integer	:: N
  integer	:: o
  real*8        :: dl,eps0,c0

! local variables
  
! START

  pml_s0=-(eps0*c0/2d0)*log(pml_r)/( dl*abs(N)**(o+1) )
  
END SUBROUTINE calc_sigma_0
!
! NAME
!     FUNCTION PML_conductivity
!
! DESCRIPTION
!     
!     work out the electrical conductivity given the distance into the PML
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 2/10/2019 CJS
!
!

FUNCTION PML_conductivity(pml_x,pml_xmin,pml_xmax,pml_s0_min,pml_s0_max,pml_order) RESULT(sigma)

real*8 sigma

integer :: pml_x,pml_xmin,pml_xmax
real*8  :: pml_s0_min,pml_s0_max
integer :: pml_order

! local variables

integer :: L    ! distance into PML
integer :: N    ! thickness of PML

! START

! work out if we are on the min or max side of the PML

if (pml_x.EQ.0) then

  sigma=0d0
  RETURN
  
else if (pml_x.LT.0) then

! min boundary

  L=abs(pml_x)-1     ! L is in the range 0 to N-1 where N is the number of layers
  N=pml_xmin
  if (N.EQ.0) then
    write(*,*)'ERROR in PML_conductivity: PML thickness is zero'
    STOP 1
  end if
  
  sigma=pml_s0_min*dble( (L+1)**(pml_order+1)-(L)**(pml_order+1) )

else 

! max boundary

  L=abs(pml_x)-1    ! L is in the range 0 to N-1 where N is the number of layers
  N=pml_xmax
  if (N.EQ.0) then
    write(*,*)'ERROR in PML_conductivity: PML thickness is zero'
    STOP 1
  end if
  
  sigma=pml_s0_max*dble( (L+1)**(pml_order+1)-(L)**(pml_order+1) )

end if

END FUNCTION PML_conductivity
!
! NAME
!     SUBROUTINE allocate_PML_material_data()
!
! DESCRIPTION
!     
!     Allocate the additional data required when material regions intersect PML layers
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 20/4/2020 CJS
!
SUBROUTINE allocate_PML_material_data()

USE Constants
USE TLM_general
USE geometry
USE PML_module
USE TLM_volume_materials
USE mesh
USE cell_parameters
USE file_information

IMPLICIT NONE

! local variables

  integer	:: cx,cy,cz
  integer	:: i
  integer       :: cell_number
  integer       :: special_cell_count
  integer       :: material_number
  integer       :: material_type
  integer       :: pml_cell
  integer       :: pml_material_cell

! START

  pml_material_cell=0
  
  do cz=nz1,nz2
    do cy=1,ny
      do cx=1,nx
      
        cell_number=cell_number+1
	
	if (cell_centre_update_code(cell_number).NE.0) then   ! special cell update
	
	  special_cell_count=cell_centre_update_code(cell_number)
          
          PML_cell=cell_update_code_to_PML_data(special_cell_count)
	  
          if (PML_cell.NE.0) then
	  
! This is a PML cell so check if this is a material cell too

	    material_number=cell_update_code_to_material_data(special_cell_count,1) 
	    if (material_number.ne.0) then
	      material_type=volume_material_list(material_number)%type
	    else
	      material_type=0
	    end if
	    
	    if (material_type.eq.volume_material_type_DISPERSIVE) then
            
	      pml_material_cell=pml_material_cell+1
                              
              PML_cell_data(PML_cell)%PML_material_array_pos=pml_material_cell
	    
	    else
	  
              PML_cell_data(PML_cell)%PML_material_array_pos=0
	  
	    end if
	    
	  end if
          
        end if
	
      end do
    end do
  end do
  
! Check that we have counted correctly...
  write(*,*)'rank=',rank
  write(*,*)'pml_material_cell                 =',pml_material_cell
  write(*,*)'total_number_of_PML_material_cells=',total_number_of_PML_material_cells

  if (pml_material_cell.NE.total_number_of_PML_material_cells) then
    write(*,*)'PML material cell count discrepancy: pml_material_cell.NE.total_number_of_PML_material_cells'
    STOP 1
  end if
  
  if (total_number_of_PML_material_cells.Eq.0) RETURN

  ALLOCATE( PML_material_cell_data(1:total_number_of_PML_material_cells) )
  
! Reset PML_cell_data

  do i=1,total_number_of_PML_material_cells
  
    PML_material_cell_data(i)%Vxt_m1=0d0
    PML_material_cell_data(i)%Vyt_m1=0d0
    PML_material_cell_data(i)%Vzt_m1=0d0
    PML_material_cell_data(i)%Ix_m1=0d0
    PML_material_cell_data(i)%Iy_m1=0d0
    PML_material_cell_data(i)%Iz_m1=0d0
    PML_material_cell_data(i)%Ixt_m1=0d0
    PML_material_cell_data(i)%Iyt_m1=0d0
    PML_material_cell_data(i)%Izt_m1=0d0
    PML_material_cell_data(i)%Vi_m1(1:12)=0d0
    PML_material_cell_data(i)%Vyxt_m1=0d0
    PML_material_cell_data(i)%Vzxt_m1=0d0
    PML_material_cell_data(i)%Vxyt_m1=0d0
    PML_material_cell_data(i)%Vzyt_m1=0d0
    PML_material_cell_data(i)%Vyzt_m1=0d0
    PML_material_cell_data(i)%Vxzt_m1=0d0
    
    PML_material_cell_data(i)%Vosxr=0d0
    PML_material_cell_data(i)%Vosyr=0d0
    PML_material_cell_data(i)%Voszr=0d0
     
    PML_material_cell_data(i)%Vosx=0d0
    PML_material_cell_data(i)%Vosy=0d0
    PML_material_cell_data(i)%Vosz=0d0
       
    PML_material_cell_data(i)%Vosx_m1=0d0
    PML_material_cell_data(i)%Vosy_m1=0d0
    PML_material_cell_data(i)%Vosz_m1=0d0
     
    PML_material_cell_data(i)%Vsx=0d0
    PML_material_cell_data(i)%Vsy=0d0
    PML_material_cell_data(i)%Vsz=0d0
     
    PML_material_cell_data(i)%Vsx_m1=0d0
    PML_material_cell_data(i)%Vsy_m1=0d0
    PML_material_cell_data(i)%Vsz_m1=0d0
     
    PML_material_cell_data(i)%Vssxr=0d0
    PML_material_cell_data(i)%Vssyr=0d0
    PML_material_cell_data(i)%Vsszr=0d0
     
    PML_material_cell_data(i)%Vcx=0d0
    PML_material_cell_data(i)%Vcy=0d0
    PML_material_cell_data(i)%Vcz=0d0
       
    PML_material_cell_data(i)%Vcx_m1=0d0
    PML_material_cell_data(i)%Vcy_m1=0d0
    PML_material_cell_data(i)%Vcz_m1=0d0
   
  end do

  RETURN

END SUBROUTINE allocate_PML_material_data
