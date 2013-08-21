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
! SUBROUTINE Output_cable_geometry()
! SUBROUTINE Output_cable_route()
! SUBROUTINE Output_bundle_geometry()
! SUBROUTINE Output_junction_specification()
! SUBROUTINE Output_face_junction_specification()


!
!SUBROUTINE Output_cable_geometry
!
! NAME
!     SUBROUTINE Output_cable_geometry
!
! DESCRIPTION
!      
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/11/2012 CJS
!
!
SUBROUTINE Output_cable_geometry()

USE TLM_general
USE Cables
USE File_information

IMPLICIT NONE

! local variables

  integer cable_geometry
  integer first_cable_geometry
  integer last_cable_geometry
  integer n_rows,n_cols
  integer row,col

! START

  CALL write_line('CALLED: Output_cable_geometry',0,output_to_screen_flag)

  write(*,*)'Number of cable geometries=',n_cable_geometries
  
10 continue

  write(*,*)
  write(*,*)'Enter the number of the cable geometry to output or 0 to output all of then'
  read(*,*)cable_geometry
  
  if (cable_geometry.eq.0) then
    first_cable_geometry=1
    last_cable_geometry=n_cable_geometries
  else 
  
    if ( (cable_geometry.lt.0).OR.(cable_geometry.gt.n_cable_geometries) ) then
      write(*,*)'Cable geometry is outside the available range'
      RETURN
    end if
    
    first_cable_geometry=cable_geometry
    last_cable_geometry=cable_geometry
  end if
  
  do cable_geometry=first_cable_geometry,last_cable_geometry
  
    write(*,*)
    write(*,*)'Cable geometry number		=',cable_geometry
    write(*,*)
    write(*,*)'Geometry type			=',trim(cable_geometry_list(cable_geometry)%cable_geometry_type_string)
    write(*,*)'Geometry type number		=',cable_geometry_list(cable_geometry)%cable_geometry_type
    write(*,*)'Number of conductors		=',cable_geometry_list(cable_geometry)%n_conductors
    write(*,*)'Number of Shielded conductors	=',cable_geometry_list(cable_geometry)%n_shielded_conductors
    write(*,*)'External_conductor_radius	=',cable_geometry_list(cable_geometry)%external_conductor_radius
    write(*,*)'External_dielectric_radius	=',cable_geometry_list(cable_geometry)%external_dielectric_radius
    write(*,*)'External_dielectric_permittivity	=',cable_geometry_list(cable_geometry)%external_dielectric_permittivity

! write parameters
    write(*,*)'Number of parameters		=',cable_geometry_list(cable_geometry)%n_parameters
    
    n_cols=cable_geometry_list(cable_geometry)%n_parameters
    write(*,*)'Parameter list:'
    write(*,*)(cable_geometry_list(cable_geometry)%parameters(col),	&
                                   col=1,cable_geometry_list(cable_geometry)%n_parameters)
				   
    n_rows=cable_geometry_list(cable_geometry)%n_conductors
    n_cols=n_rows
    
! write Sc
    write(*,*)'Shielded conductor flag vector (Sc):'
    write(*,8000)(cable_geometry_list(cable_geometry)%Sc(row),row=1,n_rows)
    
! write Tv
    write(*,*)'Voltage reference matrix (Tv):'
    do row=1,n_rows	   
      write(*,8000)(cable_geometry_list(cable_geometry)%Tv(row,col),col=1,n_cols)
      write(*,*)
    end do
    
! write Ti
    write(*,*)'Current reference matrix (Ti):'
    do row=1,n_rows	   
      write(*,8000)(cable_geometry_list(cable_geometry)%Ti(row,col),col=1,n_cols)
      write(*,*)
    end do
    
! write L_internal
    write(*,*)'Internal inductance matrix (L_internal):'
    do row=1,n_rows	   
      write(*,8010)(cable_geometry_list(cable_geometry)%L_internal(row,col),col=1,n_cols)
      write(*,*)
    end do
    
! write C_internal
    write(*,*)'Internal capacitance matrix (C_internal):'
    do row=1,n_rows	   
      write(*,8010)(cable_geometry_list(cable_geometry)%C_internal(row,col),col=1,n_cols)
      write(*,*)
    end do
    
    write(*,*)'________________________________________________________'
    write(*,*)

8000 format(100I4)
8010 format(100E15.6) 
  end do ! next cable geometry
  


  CALL write_line('FINISHED: Output_cable_geometry',0,output_to_screen_flag)

    
  RETURN
    
END SUBROUTINE Output_cable_geometry


!
!SUBROUTINE Output_cable_route
!
! NAME
!     SUBROUTINE Output_cable_route
!
! DESCRIPTION
!      
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/11/2012 CJS
!
!
SUBROUTINE Output_cable_route()

USE TLM_general
USE geometry_types
USE Cables
USE Mesh
USE File_information

IMPLICIT NONE

! local variables

  integer cable
  integer first_cable
  integer last_cable
  integer segment
  integer i
  integer number_of_cell_segments

! START

  CALL write_line('CALLED: Output_cable_route',0,output_to_screen_flag)

  write(*,*)'Number of cables=',n_cables
  
10 continue

  write(*,*)
  write(*,*)'Enter the number of the cable to output or 0 to output all of then'
  read(*,*)cable
  
  if (cable.eq.0) then
    first_cable=1
    last_cable=n_cables
  else 
  
    if ( (cable.lt.0).OR.(cable.gt.n_cables) ) then
      write(*,*)'Cable is outside the available range'
      RETURN
    end if
    
    first_cable=cable
    last_cable=cable
  end if
  
  do cable=first_cable,last_cable
  
    write(*,*)
    write(*,*)'Cable number			=',cable
    write(*,*)
    write(*,*)'Cable geometry number		=',cable_list(cable)%cable_geometry_number
    write(*,*)'End 1 junction number		=',cable_list(cable)%junction_1  
    write(*,*)'End 2 junction number		=',cable_list(cable)%junction_2  
    write(*,*)'Number of lines on cable route	=',cable_list(cable)%n_lines  
    
    write(*,*)'Cable line list:'
    write(*,8000)(cable_list(cable)%line_list(i),i=1,cable_list(cable)%n_lines) 
    
    
    write(*,*)'Number of segments on cable route	=',cable_list(cable)%number_of_cable_segments  
    write(*,*)'Cable segment list:'
    write(*,*)'     Segment End 1             Segment End 2'
    write(*,*)'  i    j    k   point       i    j    k   point'
    
    do segment=1,cable_list(cable)%number_of_cable_segments
      write(*,8010)cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%i,	&
                   cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%j,	&
                   cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%k,' ',	&
                   face_string(cable_list(cable)%cable_segment_list(segment)%segment_point(1)%point),	&
		   '    ',	&
                   cable_list(cable)%cable_segment_list(segment)%segment_point(2)%cell%i,	&
                   cable_list(cable)%cable_segment_list(segment)%segment_point(2)%cell%j,	&
                   cable_list(cable)%cable_segment_list(segment)%segment_point(2)%cell%k,' ',	&
                   face_string(cable_list(cable)%cable_segment_list(segment)%segment_point(2)%point)
                   
    end do ! next segment%i

8000 format(100I8)
8010 format(3I5,A,A,A,3I5,A,A)


    number_of_cell_segments=cable_list(cable)%number_of_cable_segments

    if (number_of_cell_segments.gt.0) then
    
! open and write line to vtk format file
      CALL open_vtk_file(cable_mesh_file_unit,cable_mesh_file_extension,cable) 
      
      CALL write_line_mesh_list_vtk(cable_mesh_file_unit,number_of_cell_segments,	&
                        cable_list(cable)%cable_segment_list)
      
      CALL close_vtk_file(cable_mesh_file_unit) 
    
    end if ! number_of_cell_segments.gt.0

  
  end do ! next cable


  CALL write_line('FINISHED: Output_cable_route',0,output_to_screen_flag)

    
  RETURN
    
END SUBROUTINE Output_cable_route

!
!SUBROUTINE Output_bundle_geometry
!
! NAME
!     SUBROUTINE Output_bundle_geometry
!
! DESCRIPTION
!      
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/11/2012 CJS
!
!
SUBROUTINE Output_bundle_geometry()

USE TLM_general
USE Cables
USE File_information

IMPLICIT NONE

! local variables

  integer bundle_geometry
  integer row
  integer col,n_cols
  integer i
  
  real*8 x,y,r
  
  character*256	:: gnuplot_filename
  character*256	:: command
  
  character*256	:: base_title
  character*256	:: bundle_title
  character*256	:: bundle_filename
  
  character	:: plot_option

! START

  CALL write_line('CALLED: Output_bundle_segment_geometry',0,output_to_screen_flag)
  
  write(*,*)'Number of cable bundle geometries=',n_bundle_segment_geometries

  write(*,*)
  write(*,*)'Enter the number of the bundle geometry to output or 0 to output a summary '
  read(*,*)bundle_geometry
  
  if (bundle_geometry.eq.0) then
    
    write(*,*)'geometry  n_cables  n_conductors	   Cable list'
    write(*,*)' number   '
      
    do bundle_geometry=1,n_bundle_segment_geometries
       
      n_cols= bundle_segment_geometry_list(bundle_geometry)%n_cables
      write(*,8000)bundle_geometry,bundle_segment_geometry_list(bundle_geometry)%n_cables,	&
                   bundle_segment_geometry_list(bundle_geometry)%n_conductors,'           ',	&
		  (bundle_segment_geometry_list(bundle_geometry)%cable_list(i),i=1,n_cols)
       
    end do
    
    RETURN    
    
  end if
  
  if ( (bundle_geometry.lt.1).OR.(bundle_geometry.gt.n_bundle_segment_geometries) ) then
  
    write(*,*)'Bundle geometry number should greater than 0 and less than ',n_bundle_segment_geometries
    RETURN
    
  end if

  write(*,*)'geometry  n_cables  n_conductors	   Cable list'
  write(*,*)' number   '
  n_cols= bundle_segment_geometry_list(bundle_geometry)%n_cables
  write(*,8000)bundle_geometry,bundle_segment_geometry_list(bundle_geometry)%n_cables,      &
  	       bundle_segment_geometry_list(bundle_geometry)%n_conductors,'	      ',    &
    	      (bundle_segment_geometry_list(bundle_geometry)%cable_list(i),i=1,n_cols)

  write(*,*)''
  write(*,*)'Number of conductors  ',bundle_segment_geometry_list(bundle_geometry)%n_conductors
  write(*,*)''
  write(*,*)' conductor     xc            yc            rc            ri'
  write(*,*)'  number    '

  do row=1,bundle_segment_geometry_list(bundle_geometry)%n_conductors
  
    write(*,8030)row,bundle_segment_geometry_list(bundle_geometry)%xc(row),	&
                     bundle_segment_geometry_list(bundle_geometry)%yc(row),	&
                     bundle_segment_geometry_list(bundle_geometry)%rc(row),	&
                     bundle_segment_geometry_list(bundle_geometry)%ri(row)
    
  end do ! next row
  
  write(*,*)''
  write(*,*)'Cable bundle radius    : ',bundle_segment_geometry_list(bundle_geometry)%cable_bundle_radius
  write(*,*)'TLM_reference_radius_rL: ',bundle_segment_geometry_list(bundle_geometry)%TLM_reference_radius_rL
  write(*,*)'TLM_reference_radius_rC: ',bundle_segment_geometry_list(bundle_geometry)%TLM_reference_radius_rC

! Write cross section data to file(s)

! write conductors
  OPEN(unit=local_file_unit,file='bundle_geometry_conductors.dat')
  
  do row=1,bundle_segment_geometry_list(bundle_geometry)%n_conductors
  
    x=bundle_segment_geometry_list(bundle_geometry)%xc(row)
    y=bundle_segment_geometry_list(bundle_geometry)%yc(row)
    r=bundle_segment_geometry_list(bundle_geometry)%rc(row)
    CALL write_circle(x,y,r,local_file_unit)
    
  end do ! next row
  
  CLOSE (unit=local_file_unit)

! write insulators
  OPEN(unit=local_file_unit,file='bundle_geometry_insulators.dat') 
  
  do row=1,bundle_segment_geometry_list(bundle_geometry)%n_conductors
  
    x=bundle_segment_geometry_list(bundle_geometry)%xc(row)
    y=bundle_segment_geometry_list(bundle_geometry)%yc(row)
    r=bundle_segment_geometry_list(bundle_geometry)%ri(row)
    CALL write_circle(x,y,r,local_file_unit)
    
  end do ! next row
  
  CLOSE (unit=local_file_unit)

! write TLM_cell equivalent radius
  OPEN(unit=local_file_unit,file='bundle_geometry_TLM_cell.dat')

  x=0d0
  y=0d0
  r=bundle_segment_geometry_list(bundle_geometry)%TLM_cell_equivalent_radius
  CALL write_circle_dash(x,y,r,local_file_unit)
  
  CLOSE (unit=local_file_unit)

! write Inductance_reference
  OPEN(unit=local_file_unit,file='bundle_geometry_Inductance_reference.dat')

  x=0d0
  y=0d0
  r=bundle_segment_geometry_list(bundle_geometry)%TLM_reference_radius_rL
  CALL write_circle(x,y,r,local_file_unit)
  
  CLOSE (unit=local_file_unit)

! write Capacitance_reference
  OPEN(unit=local_file_unit,file='bundle_geometry_Capacitance_reference.dat')

  x=0d0
  y=0d0
  r=bundle_segment_geometry_list(bundle_geometry)%TLM_reference_radius_rC
  CALL write_circle(x,y,r,local_file_unit)
  
  CLOSE (unit=local_file_unit)

! Write gnuplot file to plot this data    

  gnuplot_filename='cable_plot.plt'
    
  write(*,*)'Enter the output required s(creen) g(if) j(peg)'
  read(*,'(A)')plot_option
  
  CALL convert_to_lower_case(plot_option,1)
  
  base_title="Bundle_geometry_number"

  CALL add_integer_to_filename(base_title,bundle_geometry,bundle_title)
  
  OPEN(unit=local_file_unit,file=gnuplot_filename)
  
  if (plot_option.eq.'g') then
    bundle_filename=trim(bundle_title)//'.gif'
    write(local_file_unit,'(A)')'set term gif'
    write(local_file_unit,'(A,A,A)')'set output "',trim(bundle_filename),'"' 
  else if (plot_option.eq.'j') then
    bundle_filename=trim(bundle_title)//'.jpeg'
    write(local_file_unit,'(A)')'set term jpeg'
    write(local_file_unit,'(A,A,A)')'set output "',trim(bundle_filename),'"' 
  end if
 
  write(local_file_unit,'(A)')'set size square'
  write(local_file_unit,'(A)')'set autoscale x'
  write(local_file_unit,'(A)')'set autoscale y'
  write(local_file_unit,'(A,I8,A)')'set title "Bundle geometry number',bundle_geometry,'"'
  write(local_file_unit,'(A)')'set xlabel "x"'
  write(local_file_unit,'(A)')'set ylabel "y"'
  
  write(local_file_unit,'(A)')'plot "bundle_geometry_conductors.dat" u 1:2 title "Conductors" w l,\'
  write(local_file_unit,'(A)')'     "bundle_geometry_insulators.dat" u 1:2 title "Insulators" w l,\'
  write(local_file_unit,'(A)')'     "bundle_geometry_TLM_cell.dat" u 1:2 title "TLM cell" w l,\'
  write(local_file_unit,'(A)')'     "bundle_geometry_Inductance_reference.dat" u 1:2 title "Inductance_reference" w l,\'
  write(local_file_unit,'(A)')'     "bundle_geometry_Capacitance_reference.dat" u 1:2 title "Capacitance_reference" w l'
  
  if (plot_option.eq.'s') then
    write(local_file_unit,'(A)')'pause 100'
  end if
  
  CLOSE(unit=local_file_unit)
    
  command='gnuplot cable_plot.plt & '
  CALL system(command)

! Write L,C materices to screen

  n_cols=bundle_segment_geometry_list(bundle_geometry)%n_conductors
  write(*,*)'L'    
  do row=1,n_cols
    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%L(row,col),col=1,n_cols)
  end do ! next row
    
  write(*,*)'C'    
  do row=1,n_cols
    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%C(row,col),col=1,n_cols)
  end do ! next row
    
  write(*,*)'R'    
  do row=1,n_cols
    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%R(row,col),col=1,n_cols)	 
  end do ! next row
    
  write(*,*)'Zlink'    
  do row=1,n_cols
    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%Zlink(row,col),col=1,n_cols)	 
  end do ! next row
    
  write(*,*)'Ylink'    
  do row=1,n_cols
    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%Ylink(row,col),col=1,n_cols)	 
  end do ! next row
    
  write(*,*)'ZLstub'    
  do row=1,n_cols
    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%ZLstub(row,col),col=1,n_cols)	 
  end do ! next row
    
  write(*,*)'Yf'    
  do row=1,n_cols
    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%Yf(row,col),col=1,n_cols)	 
  end do ! next row
    
  write(*,*)'Tv'    
  do row=1,n_cols
    write(*,8020)(bundle_segment_geometry_list(bundle_geometry)%Tv(row,col),col=1,n_cols)     
  end do ! next row
    
  write(*,*)'Ti'    
  do row=1,n_cols
    write(*,8020)(bundle_segment_geometry_list(bundle_geometry)%Ti(row,col),col=1,n_cols)     
  end do ! next row
    
  write(*,*)'Sc'    
  do row=1,n_cols
    write(*,*)bundle_segment_geometry_list(bundle_geometry)%SC(row)
  end do ! next row
  
8000  format(I7,I9,I12,A,100I5)
8010  format(100E14.4)
8020  format(100I2)
8030  format(I7,4E14.4)

  CALL write_line('FINISHED: Output_bundle_geometry',0,output_to_screen_flag)

    
  RETURN
    
END SUBROUTINE Output_bundle_geometry

!
!SUBROUTINE Output_junction_specification
!
! NAME
!     SUBROUTINE Output_junction_specification
!
! DESCRIPTION
!      
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/11/2012 CJS
!
!
SUBROUTINE Output_junction_specification()

USE TLM_general
USE cell_parameters
USE Cables
USE File_information

IMPLICIT NONE

! local variables

  integer bundle_junction
  integer ix,iy,iz,face
  integer segment
  integer n_internal,n_external
  integer row,col
  integer i
  integer filter

! START

  CALL write_line('CALLED: Output_junction_specification',0,output_to_screen_flag)

  write(*,*)'Number of bundle junctions=',n_cell_centre_junctions
  
10 continue

  write(*,*)
  write(*,*)'Enter the number of the cell centre junction to output or 0 to output a summary '
  read(*,*)bundle_junction
  
  if (bundle_junction.eq.0) then
  
    write(*,*)
    write(*,*)'junction   ix   iy   iz n_internal segment segment segment segment segment segment'
    write(*,*)'                                     xmin    ymin    zmin    xmax    ymax    zmax '
    
    do bundle_junction=1,n_cell_centre_junctions
    
     ix=cell_centre_junction_list(bundle_junction)%cell_point%cell%i
     iy=cell_centre_junction_list(bundle_junction)%cell_point%cell%j
     iz=cell_centre_junction_list(bundle_junction)%cell_point%cell%k
     n_internal=cell_centre_junction_list(bundle_junction)%n_internal_connection_nodes
    
      write(*,8000)bundle_junction,ix,iy,iz,n_internal,	&
      cell_centre_junction_list(bundle_junction)%segment_list(1),	&
      cell_centre_junction_list(bundle_junction)%segment_list(2),	&
      cell_centre_junction_list(bundle_junction)%segment_list(3),	&
      cell_centre_junction_list(bundle_junction)%segment_list(4),	&
      cell_centre_junction_list(bundle_junction)%segment_list(5),	&
      cell_centre_junction_list(bundle_junction)%segment_list(6)
8000  format(I7,I7,2I5,I11,6I8)
   
    end do
  end if
  
  if ( (bundle_junction.lt.1).OR.(bundle_junction.gt.n_cell_centre_junctions) ) then
    write(*,*)'Bundle junction number should greater than 0 and less than ',n_cell_centre_junctions
    RETURN
  end if

  ix=cell_centre_junction_list(bundle_junction)%cell_point%cell%i
  iy=cell_centre_junction_list(bundle_junction)%cell_point%cell%j
  iz=cell_centre_junction_list(bundle_junction)%cell_point%cell%k

  n_internal=cell_centre_junction_list(bundle_junction)%n_internal_connection_nodes

  write(*,*)' Bundle junction number',bundle_junction
  write(*,*)
  write(*,*)'ix=',ix,' iy=',iy,' iz=',iz
  write(*,*)
  write(*,*)'Number of internal connection nodes=',n_internal
  write(*,*)
  do face=1,6
  
    segment=cell_centre_junction_list(bundle_junction)%segment_list(face)

    write(*,*)'Face:',face_string(face),' Segment number=',segment
    if (segment.ne.0) then
    
      n_external=cell_centre_junction_list(bundle_junction)%n_external_conductors(face)
    
      write(*,*)'Number of external cables=',bundle_segment_list(segment)%n_cables
      write(*,*)'Number of external conductors=',bundle_segment_list(segment)%n_conductors,' =',n_external
      
      write(*,*)'Cable list:'
      
      do i=1,bundle_segment_list(segment)%n_cables
        write(*,*)bundle_segment_list(segment)%cable_list(i)
      end do
      
      write(*,*)'P matrix (n_internal rows x n_external columns):'
      write(*,*)
      
      do row=1,n_internal
      
        write(*,8010)(cell_centre_junction_list(bundle_junction)%P_matrix_list(face)%P(row,col),col=1,n_external)
        write(*,*)
8010	format(1000I2)
	
      end do ! next internal connection node
      
    end if ! there is a cable bundle on this face

  end do ! next face
  
  face=7
  n_external=cell_centre_junction_list(bundle_junction)%n_external_conductors(face)
     
  if (n_external.ne.0) then
    write(*,*)'Internal impedance data'	    
    write(*,*)'Number of impedance filters=',n_external	    
    write(*,*)'P matrix (n_internal rows x n_external columns):'
    write(*,*)
    do row=1,n_internal
    
      write(*,8010)(cell_centre_junction_list(bundle_junction)%P_matrix_list(face)%P(row,col),col=1,n_external)
      write(*,*)

    end do ! next internal connection node
    do filter=1,cell_centre_junction_list(bundle_junction)%n_internal_impedance_filters

      write(*,*)'! Sfilter number',filter	    
      CALL write_Sfilter(cell_centre_junction_list(bundle_junction)%Sfilter(filter),0)
      write(*,*)'! Z_f',filter	  
      write(*,*)cell_centre_junction_list(bundle_junction)%Z_f(filter)
      write(*,*)'! Zfilter number',filter	 
      CALL write_Zfilter(cell_centre_junction_list(bundle_junction)%Zfilter(filter),0)
    
    end do ! next filter
  end if
  
  CALL write_line('FINISHED: Output_junction_specification',0,output_to_screen_flag)

    
  RETURN
    
END SUBROUTINE Output_junction_specification

!
!SUBROUTINE Output_face_junction_specification
!
! NAME
!     SUBROUTINE Output_face_junction_specification
!
! DESCRIPTION
!      
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/11/2012 CJS
!
!
SUBROUTINE Output_face_junction_specification()

USE TLM_general
USE Cables
USE File_information

IMPLICIT NONE

! local variables

  integer face_junction
  integer ix,iy,iz,face
  integer segment
  integer n_internal,n_external
  integer row,col
  integer i

! START

  CALL write_line('CALLED: Output_face_junction_specification',0,output_to_screen_flag)

  write(*,*)'Number of face junctions=',n_face_junctions
  
10 continue

  write(*,*)
  write(*,*)'Enter the number of the face junction to output or 0 to output a summary '
  read(*,*)face_junction
  
  if (face_junction.eq.0) then
  
    write(*,*)
    write(*,*)'junction   ix   iy   iz n_internal segment1 segment2 '
    
    do face_junction=1,n_face_junctions
    
     ix=face_junction_list(face_junction)%cell_point%cell%i
     iy=face_junction_list(face_junction)%cell_point%cell%j
     iz=face_junction_list(face_junction)%cell_point%cell%k
     n_internal=face_junction_list(face_junction)%n_internal_connection_nodes
    
      write(*,8000)face_junction,ix,iy,iz,n_internal,	&
      face_junction_list(face_junction)%segment_list(1),	&
      face_junction_list(face_junction)%segment_list(2)
8000  format(I7,3I5,I11,2I9)
   
    end do
  end if
  
  if ( (face_junction.lt.1).OR.(face_junction.gt.n_face_junctions) ) then
    write(*,*)'Bundle junction number should greater than 0 and less than ',n_face_junctions
    RETURN
  end if

  ix=face_junction_list(face_junction)%cell_point%cell%i
  iy=face_junction_list(face_junction)%cell_point%cell%j
  iz=face_junction_list(face_junction)%cell_point%cell%k

  n_internal=face_junction_list(face_junction)%n_internal_connection_nodes

  write(*,*)' Bundle junction number',face_junction
  write(*,*)
  write(*,*)'ix=',ix,' iy=',iy,' iz=',iz
  write(*,*)
  write(*,*)'Number of internal connection nodes=',n_internal
  write(*,*)
  do face=1,2
  
    segment=face_junction_list(face_junction)%segment_list(face)

    write(*,*)'Face:',face,' Segment number=',segment
    if (segment.ne.0) then
    
      n_external=face_junction_list(face_junction)%n_external_conductors(face)
    
      write(*,*)'Number of external cables=',bundle_segment_list(segment)%n_cables
      write(*,*)'Number of external conductors=',bundle_segment_list(segment)%n_conductors,' =',n_external
      
      write(*,*)'Cable list:'
      
      do i=1,bundle_segment_list(segment)%n_cables
        write(*,*)bundle_segment_list(segment)%cable_list(i)
      end do
      
      write(*,*)'P matrix (n_internal rows x n_external columns):'
      write(*,*)
      
      do row=1,n_internal
      
        write(*,8010)(face_junction_list(face_junction)%P_matrix_list(face)%P(row,col),col=1,n_external)
        write(*,*)
8010	format(1000I2)
	
      end do ! next internal connection node
      
      write(*,*)'BC vector (n_internal rows ):'
      write(*,*)
      
      do row=1,n_internal
      
        write(*,8020)face_junction_list(face_junction)%BC(row)
        write(*,*)
8020	format(I4)
	
      end do ! next internal connection node
      
    end if ! there is a cable bundle on this face

  end do ! next face

  CALL write_line('FINISHED: Output_face_junction_specification',0,output_to_screen_flag)

    
  RETURN
    
END SUBROUTINE Output_face_junction_specification

!
!SUBROUTINE write_circle
!
! NAME
!     SUBROUTINE write_circle
!
! DESCRIPTION
!     write a circle with specified x,y centre and radius to file for plotting with gnuplot
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 10/05/2013 CJS
!
!
  SUBROUTINE write_circle(x,y,r,unit)
   
USE constants	

IMPLICIT NONE

  real*8 x,y,r
  integer unit
  
! local variables  

  real*8 t
  real*8 xp,yp

  integer tloop
  
! START

! loop over theta    

  do tloop=0,50
  
    t=2d0*pi*dble(tloop)/50d0

    xp=x+r*cos(t)
    yp=y+r*sin(t)

    write(unit,8000)xp,yp
8000 format (4E14.6)
  
  end do
  
  write(unit,*)
  write(unit,*)
   
  return
  
  end SUBROUTINE write_circle

!
!SUBROUTINE write_circle_dash
!
! NAME
!     SUBROUTINE write_circle_dash
!
! DESCRIPTION
!     write a dashed circle with specified x,y centre and radius to file for plotting with gnuplot
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 10/05/2013 CJS
!
!
  SUBROUTINE write_circle_dash(x,y,r,unit)
   
USE constants	

IMPLICIT NONE

  real*8 x,y,r
  integer unit
  
! local variables  

  real*8 t,t2
  real*8 xp,yp

  integer tloop
  
! START

! loop over theta    

  do tloop=0,50
  
    t=2d0*pi*dble(tloop)/50d0

    xp=x+r*cos(t)
    yp=y+r*sin(t)

    write(unit,8000)xp,yp
8000 format (4E14.6)

    t2=2d0*pi*(dble(tloop)+0.5d0)/50d0
    
    xp=x+r*cos(t2)
    yp=y+r*sin(t2)

    write(unit,8000)xp,yp
	
    write(unit,*)
    write(unit,*)
    
  end do
  
    
  return
  
  end SUBROUTINE write_circle_dash
