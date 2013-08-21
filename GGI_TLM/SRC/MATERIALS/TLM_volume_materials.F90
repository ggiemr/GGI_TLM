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
! SUBROUTINE set_volume_material_mesh
! SUBROUTINE calculate_volume_material_filter_coefficients
! SUBROUTINE allocate_volume_material_filter_data
! SUBROUTINE volume_material_update
!
! NAME
!     set_volume_material_mesh
!
! DESCRIPTION
!     
!     Count the number of cells of each material and construct the 
!     appropriate lists for the main solver
!     
! COMMENTS
!     Allocate a local mesh array and fill with material data. This ensures that
!     a cell doesn't get given more than one material property.
!     In overlapping volumes the last material assignment takes precidence
!
! HISTORY
!
!     started 5/09/2012 CJS
!     frequency warping included 12/02/2013 CJS
!
!
SUBROUTINE set_volume_material_mesh

USE TLM_general
USE geometry
USE TLM_volume_materials
USE mesh
USE cell_parameters
USE file_information

IMPLICIT NONE

! local variables

  integer	:: material_number
  integer	:: volume_number
  integer	:: cell,number_of_cells
  integer	:: total_number_of_material_cells
  integer	:: total_number_of_material_cells_all_procsesses
  integer	:: n_data
  integer	:: cx,cy,cz
  
  integer	:: i
  integer	:: cell_number

! START
  
  CALL write_line('CALLED: set_volume_material_mesh',0,output_to_screen_flag)

! INITIALISE VOLUME MATERIALS  
    
  CALL write_line_integer('n_volume_materials=',n_volume_materials,0,output_to_screen_flag)

  if (rank.eq.0) then
  
    write(info_file_unit,*)'____________________________________________________'
    write(info_file_unit,*)''
    write(info_file_unit,*)'Number of volume materials=',n_volume_materials
    write(info_file_unit,*)''
    write(info_file_unit,*)'Material_number  Total_number_of_material_cells'
  
  end if

  do material_number=1,n_volume_materials
  
    total_number_of_material_cells=0
    
    CALL write_line_integer('volume material number=',material_number,0,output_to_screen_flag)
    CALL write_line_integer('Number of geometric volumes=',	&
                             volume_material_list(material_number)%n_volumes,0,output_to_screen_flag)

! loop over the geometric volumes with this material type and initialise the TLM cells
    do i=1,volume_material_list(material_number)%n_volumes
    
      volume_number=volume_material_list(material_number)%volume_list(i)

      CALL write_line_integer('Geometric volume material number=',volume_number,0,output_to_screen_flag)
      
      problem_volumes(volume_number)%volume_material_number=material_number
      
      number_of_cells=problem_volumes(volume_number)%number_of_cells

      CALL write_line_integer('Number of cells=',number_of_cells,0,output_to_screen_flag)
      
      total_number_of_material_cells=total_number_of_material_cells+number_of_cells
      
      do cell=1,number_of_cells
    
	cx=problem_volumes(volume_number)%cell_list(cell)%cell%i
	cy=problem_volumes(volume_number)%cell_list(cell)%cell%j
	cz=problem_volumes(volume_number)%cell_list(cell)%cell%k
	
	if (rank.eq.cell_rank(cz)) then
	
          local_cell_material(cx,cy,cz)=material_number
	  
	end if
	
      end do ! next cell
    
    end do ! next volume number

#if defined(MPI)

    n_data=1
    call MPI_REDUCE(total_number_of_material_cells,total_number_of_material_cells_all_procsesses , &
                    n_data,MPI_INTEGER, MPI_SUM, 0,MPI_COMM_WORLD,ierror)
		    
#elif defined(SEQ)

    total_number_of_material_cells_all_procsesses=total_number_of_material_cells
    
#endif			

    if (rank.eq.0) then
      CALL write_line_integer('Total number of material cells=',	&
                              total_number_of_material_cells_all_procsesses,0,output_to_screen_flag)
      write(info_file_unit,'(2I14)')material_number,total_number_of_material_cells_all_procsesses
    end if

  end do ! next volume material number
  
  CALL write_line('FINISHED: set_volume_material_mesh',0,output_to_screen_flag)

  RETURN

END SUBROUTINE set_volume_material_mesh
!
! Name calculate_volume_material_filter_coefficients
!     
!
! Description:
!      Once the timestep is found we can determine the volume material filter coefficients
!      
!
! Comments:
!      
!      
!
! History
!
!     started 04/09/12 CJS
!

SUBROUTINE calculate_volume_material_filter_coefficients()

USE TLM_general
USE TLM_volume_materials
USE filter_types
USE filter_operators
USE filter_functions
USE constants

IMPLICIT NONE

! variables passed to subroutine

! local_variables

  integer material_number
  
  type(Sfilter) :: Zs_eps_Sfilter
  type(Zfilter) :: Zs_eps_Zfilter1,Zs_eps_Zfilter2
  real*8	:: Zs_eps_f
  real*8	:: Cstub_constant
  
  type(Sfilter) :: Zs_mu_Sfilter
  type(Zfilter) :: Zs_mu_Zfilter1,Zs_mu_Zfilter2
  real*8	:: Zs_mu_f
  real*8	:: Lstub_constant

! function_types

! START

  write(*,*)'CALLED: calculate_volume_material_filter_coefficients'
  
  do material_number=1,n_volume_materials
    
    if (volume_material_list(material_number)%type.EQ.volume_material_type_DISPERSIVE) then
    
! Calculate Filter coefficients for the frequency dependent impedance of electric susceptibility
      
! Electric loss admittance
      volume_material_list(material_number)%Ge=volume_material_list(material_number)%sigma_e*dl
      if (volume_material_list(material_number)%Ge.lt.0d0) then
! if we think this is a result of tolerance issuse then set to zero, else cause an error
        if (volume_material_list(material_number)%Ge.lt.-small) then
            GOTO 9010
	else
	  volume_material_list(material_number)%Ge=0d0
        end if
      end if

      write(*,*)'Ge=',volume_material_list(material_number)%Ge

      write(*,*)'eps_S filter in'
      call write_Sfilter(volume_material_list(material_number)%eps_S,0)

! see if the filter exists or is epsr=1d0 i.e. free space
      if ( (volume_material_list(material_number)%eps_S%a%order.EQ.0) .AND.	&
           (volume_material_list(material_number)%eps_S%a%coeff(0).EQ.1d0) ) then ! free space properties
	   
        volume_material_list(material_number)%eps_filter_exists=.FALSE.
        volume_material_list(material_number)%Ys_eps_f=0d0
	
      else
      
        volume_material_list(material_number)%eps_filter_exists=.TRUE.     
	
      end if

      if (volume_material_list(material_number)%eps_filter_exists) then

! Stage 1. Calculate susceptibility impedance= 1/( jw*eps0*dl*(epsr(jw)-1) )

        Cstub_constant=eps0*dl
          
        CALL calculate_electric_susceptibility_impedance_Sfilter(		&
                     volume_material_list(material_number)%eps_S,Cstub_constant,Zs_eps_Sfilter ) 
        write(*,*)'Zs_eps_S filter'
        call write_Sfilter(Zs_eps_Sfilter,0)


! Stage 2. calculate Z domain susceptibility admittance filter coefficients by bilinear transformation
!        Zs_eps_Zfilter1=s_to_z(Zs_eps_Sfilter,dt) 
        Zs_eps_Zfilter1=s_to_z_warp(Zs_eps_Sfilter,dt,bicubic_warp_flag,frequency_scale) 
          
        write(*,*)'Zs_eps_Z filter1'
        call write_Zfilter(Zs_eps_Zfilter1,0)
    
! Stage 3. calculate fast permittivity and slow permittivity filter

        call Z_fast_slow_docomposition( Zs_eps_Zfilter1 ,Zs_eps_f  , Zs_eps_Zfilter2 )
          
        write(*,*)'Zs_eps_f',Zs_eps_f
        write(*,*)'Ys_eps_f',1d0/Zs_eps_f
        write(*,*)'Zs_eps_Z filter2'
        call write_Zfilter(Zs_eps_Zfilter2,0)
    
        volume_material_list(material_number)%Zs_eps_Z=Zs_eps_Zfilter2
        volume_material_list(material_number)%Zs_eps_f=Zs_eps_f
        volume_material_list(material_number)%Ys_eps_f=1d0/Zs_eps_f
 
      end if ! eps filter exists i.e. epsr.NE.1
      
! Magnetic loss impedance
      volume_material_list(material_number)%Rm=volume_material_list(material_number)%sigma_m*dl
      if (volume_material_list(material_number)%Rm.lt.0d0) then
! if we think this is a result of tolerance issuse then set to zero, else cause an error
        if (volume_material_list(material_number)%Rm.lt.-small) then
          GOTO 9030
	else
	  volume_material_list(material_number)%Rm=0d0
        end if
      end if

      write(*,*)'Rm=',volume_material_list(material_number)%Rm
    
! Calculate Filter coefficients for the frequency dependent impedance of magnetic susceptibility

      write(*,*)'mu_S filter in'
      call write_Sfilter(volume_material_list(material_number)%mu_S,0)
      
! see if the filter exists or is epsr=1d0 i.e. free space
      if ( (volume_material_list(material_number)%mu_S%a%order.EQ.0) .AND.	&
           (volume_material_list(material_number)%mu_S%a%coeff(0).EQ.1d0) ) then ! free space properties
	   
        volume_material_list(material_number)%mu_filter_exists=.FALSE.
        volume_material_list(material_number)%Zs_mu_f=0d0
	
      else
      
        volume_material_list(material_number)%mu_filter_exists=.TRUE.     
	
      end if

      if (volume_material_list(material_number)%mu_filter_exists) then

! Stage 1. Calculate susceptibility impedance= jw*mu0*dl*(mur(jw)-1)

        Lstub_constant=mu0*dl
      
        CALL calculate_magnetic_susceptibility_impedance_Sfilter(		&
                   volume_material_list(material_number)%mu_S,Lstub_constant,Zs_mu_Sfilter)
		   
        write(*,*)'Zs_mu_S filter'
        call write_Sfilter(Zs_mu_Sfilter,0)


! Stage 2. calculate Z domain susceptibility admittance filter coefficients by bilinear transformation
!        Zs_mu_Zfilter1=s_to_z(Zs_mu_Sfilter,dt) 
        Zs_mu_Zfilter1=s_to_z_warp(Zs_mu_Sfilter,dt,bicubic_warp_flag,frequency_scale) 
          
        write(*,*)'Zs_mu_Z filter1'
        call write_Zfilter(Zs_mu_Zfilter1,0)
    
! Stage 3. calculate fast permeability and slow permeability filter

        call Z_fast_slow_docomposition( Zs_mu_Zfilter1 ,Zs_mu_f  , Zs_mu_Zfilter2 )
          
        write(*,*)'Zs_mu_f',Zs_mu_f
        write(*,*)'Zs_mu_Z filter2'
        call write_Zfilter(Zs_mu_Zfilter2,0)
    
        volume_material_list(material_number)%Zs_mu_Z=Zs_mu_Zfilter2
        volume_material_list(material_number)%Zs_mu_f=Zs_mu_f
		
      end if ! permeability filter exists i.e. mur.NE.1
		     
    end if ! material type is dispersive
    
  end do ! next material to set

  write(*,*)'FINISHED: calculate_volume_material_filter_coefficients'
  
  RETURN
  
9000 CALL write_line('Error in :calculate_volume_material_filter_coefficients',0,.TRUE.)
     CALL write_line_integer('Negative Capacitive stub, material number',material_number,0,.TRUE.)
     write(*,*)volume_material_list(material_number)%Ys_eps_f
     STOP
  
9010 CALL write_line('Error in :calculate_volume_material_filter_coefficients',0,.TRUE.)
     CALL write_line_integer('Negative Electric loss, material number',material_number,0,.TRUE.)
     write(*,*)volume_material_list(material_number)%Ge
     STOP
  
9020 CALL write_line('Error in :calculate_volume_material_filter_coefficients',0,.TRUE.)
     CALL write_line_integer('Negative Inductive stub, material number',material_number,0,.TRUE.)
     write(*,*)volume_material_list(material_number)%Zs_mu_f
     STOP
  
9030 CALL write_line('Error in :calculate_volume_material_filter_coefficients',0,.TRUE.)
     CALL write_line_integer('Negative Magnetic loss, material number',material_number,0,.TRUE.)
     write(*,*)volume_material_list(material_number)%Rm
     STOP

  
END SUBROUTINE calculate_volume_material_filter_coefficients
!
! NAME
!     allocate_volume_material_filter_data
!
! DESCRIPTION
!     allocate the memory required for volume material filter data
!     
! COMMENTS
!
!
! HISTORY
!
!     started 4/09/2012 CJS
!
!
SUBROUTINE allocate_volume_material_filter_data

USE TLM_general
USE TLM_volume_materials
USE filter_types
USE filter_operators
USE filter_functions
USE mesh

IMPLICIT NONE

! local variables
  
  integer	:: material_number
  integer	:: cell
  integer	:: volume_filter_number
  integer       :: epsaorder,epsborder
  integer       :: muaorder,muborder

! START
  
  CALL write_line('CALLED: allocate_volume_material_filter_data',0,output_to_screen_flag)
  
  n_volume_material_cells=0
  volume_material_storage=0

! loop over special cells
  DO cell=1,n_special_cells
  
    if (cell_update_code_to_material_data(cell,1).ne.0) then
! this is a material cell add to the memory allocation list if required

      material_number=cell_update_code_to_material_data(cell,1)
      
      if ( (volume_material_list(material_number)%type.eq.volume_material_type_DISPERSIVE) ) then
	  	      
        n_volume_material_cells=n_volume_material_cells+1

        cell_update_code_to_material_data(cell,2)=volume_material_storage+1
        volume_material_storage=volume_material_storage+3   ! note 3 filters per cell, one for each field component    

      end if ! volume material at this cell cell

    end if ! is this a material cell? 

  end do ! next special cell

! allocate volume material filter data storage for each filter required
  if (volume_material_storage.ne.0) then
  
    allocate (volume_material_Zs_eps_filter_data(1:volume_material_storage)) 
    allocate (volume_material_Zs_mu_filter_data(1:volume_material_storage)) 
    
  end if
 
! loop over filters allocating memory for the particular material order

! loop over special cells
  DO cell=1,n_special_cells
  
    if (cell_update_code_to_material_data(cell,1).ne.0) then
! this is a material cell add to the memory allocation list if required

      material_number=cell_update_code_to_material_data(cell,1)
      
      if ( (volume_material_list(material_number)%type.eq.volume_material_type_DISPERSIVE) ) then	  	      
      
        if (volume_material_list(material_number)%eps_filter_exists) then

	  epsaorder=volume_material_list(material_number)%Zs_eps_Z%a%order
	  epsborder=volume_material_list(material_number)%Zs_eps_Z%b%order

! x polarisation 
          volume_filter_number=cell_update_code_to_material_data(cell,2)	  
          volume_material_Zs_eps_filter_data(volume_filter_number)=allocate_Zfilter_response(epsaorder,epsborder)	

! y polarisation 
          volume_filter_number=cell_update_code_to_material_data(cell,2)+1	  
          volume_material_Zs_eps_filter_data(volume_filter_number)=allocate_Zfilter_response(epsaorder,epsborder)	

! z polarisation 
          volume_filter_number=cell_update_code_to_material_data(cell,2)+2	  
          volume_material_Zs_eps_filter_data(volume_filter_number)=allocate_Zfilter_response(epsaorder,epsborder)	
      
        end if ! permittivity filter exists
      
        if (volume_material_list(material_number)%mu_filter_exists) then
     
	  muaorder =volume_material_list(material_number)%Zs_mu_Z%a%order
	  muborder =volume_material_list(material_number)%Zs_mu_Z%b%order

! x polarisation 
          volume_filter_number=cell_update_code_to_material_data(cell,2)	  
          volume_material_Zs_mu_filter_data(volume_filter_number) =allocate_Zfilter_response(muaorder,muborder)	

! y polarisation 
          volume_filter_number=cell_update_code_to_material_data(cell,2)+1	  
          volume_material_Zs_mu_filter_data(volume_filter_number) =allocate_Zfilter_response(muaorder,muborder)	

! z polarisation 
          volume_filter_number=cell_update_code_to_material_data(cell,2)+2	  
          volume_material_Zs_mu_filter_data(volume_filter_number)=allocate_Zfilter_response(muaorder,muborder)	
       
        end if ! permeability filter exists
      
      end if ! volume material at this cell

    end if ! is this a material cell? 

  end do ! next special cell
 
 
  CALL write_line('FINISHED: allocate_volume_material_filter_data',0,output_to_screen_flag)


  RETURN

END SUBROUTINE allocate_volume_material_filter_data
