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
!
! NAME
!     SUBROUTINE set_terminal_connection_cells
!     SUBROUTINE mesh_to_stl
!     SUBROUTINE get_TLM_cell_centre_coordinate
!     SUBROUTINE get_TLM_cell_from_coordinate
!
! DESCRIPTION
!     Create a stl file for the geometry of all the surfaces set in the mesh
!     
! COMMENTS
!     
!
!
!
! HISTORY
!
!     started 14/3/19 CJS
!     
!

SUBROUTINE set_terminal_connection_cells(x1,y1,z1,x2,y2,z2)

USE PCB_simulation

IMPLICIT NONE

! local variables

real*8  :: x1,y1,z1,x2,y2,z2

integer :: ix1,iy1,iz1,ix2,iy2,iz2

integer :: ix,iy,iz,dx,dy,dz,nx,ny,nz

integer :: i

! START

! get the cell coordinates of the end points

  CALL get_TLM_cell_from_coordinate(x1,y1,z1,ix1,iy1,iz1)
  CALL get_TLM_cell_from_coordinate(x2,y2,z2,ix2,iy2,iz2)

  dx=0
  dy=0
  dz=0

  if (ix2.GT.ix1) then
    dx=1
  else if (ix2.LT.ix1) then
    dx=-1
  end if
  nx=abs(ix2-ix1)
  
  if (iy2.GT.iy1) then
    dy=1
  else if (iy2.LT.iy1) then
    dy=-1
  end if
  ny=abs(iy2-iy1)
  
  if (iz2.GT.iz1) then
    dz=1
  else if (iz2.LT.iz1) then
    dz=-1
  end if
  nz=abs(iz2-iz1)
  write(*,*)'ix1=',ix1,' iy1=',iy1,' iz1=',iz1
  write(*,*)'ix2=',ix2,' iy2=',iy2,' iz2=',iz2
  write(*,*)'dx =',dx, ' dy =',dy ,' dz =',dz
  
! set the first cell
  ix=ix1
  iy=iy1
  iz=iz1
  write(*,*)'Setting cell:',ix,iy,iz
  material_mesh(centre,ix,iy,iz)=1
  
! Move in the z direction first

  do i=1,nz
    iz=iz+dz
    write(*,*)'Setting cell:',ix,iy,iz
    material_mesh(centre,ix,iy,iz)=1
  end do

! move in y
  do i=1,ny
    iy=iy+dy
    write(*,*)'Setting cell:',ix,iy,iz
    material_mesh(centre,ix,iy,iz)=1
  end do

! move in x
  do i=1,nx
    ix=ix+dx
    write(*,*)'Setting cell:',ix,iy,iz
    material_mesh(centre,ix,iy,iz)=1
  end do
  
! set the last cell
    write(*,*)'Setting cell:',ix2,iy2,iz2
  material_mesh(centre,ix2,iy2,iz2)=1

RETURN  
  
END SUBROUTINE set_terminal_connection_cells
!
! DESCRIPTION
!     Create a stl file for the geometry of all the surfaces set in the mesh
!     
! COMMENTS
!     
!
!
!
! HISTORY
!
!     started 14/3/19 CJS
!     
!

SUBROUTINE mesh_to_stl

USE PCB_simulation

IMPLICIT NONE

! local variables

integer :: i,ix,iy,iz,loop

real*8  :: cx,cy,cz

integer :: n_faces,face
integer :: n_triangles
integer :: n_points
real*8,allocatable :: node_list(:,:)

! START

write(*,*)'CALLED: mesh_to_stl'

open(unit=40,file=terminal_connection_geometry_filename)

write(40,'(A)')'solid ascii'

! faces normal to x

do loop=1,2

  n_faces=0
  n_points=0
  
  do ix=1,mesh_nx
    do iy=1,mesh_ny
      do iz=1,mesh_nz

! Look for a discontinuity in material in x and place a surface there    
        if (material_mesh(centre,ix,iy,iz).NE.material_mesh(centre,ix+1,iy,iz)) then
                  
          if (loop.EQ.2) then
          
            n_faces=n_faces+1
            
            CALL get_TLM_cell_centre_coordinate(ix,iy,iz,cx,cy,cz)
          
            n_points=n_points+1
            node_list(n_points,1)=cx+dl/2d0
            node_list(n_points,2)=cy-dl/2d0
            node_list(n_points,3)=cz-dl/2d0
            
            n_points=n_points+1
            node_list(n_points,1)=cx+dl/2d0
            node_list(n_points,2)=cy+dl/2d0
            node_list(n_points,3)=cz-dl/2d0
            
            n_points=n_points+1
            node_list(n_points,1)=cx+dl/2d0
            node_list(n_points,2)=cy+dl/2d0
            node_list(n_points,3)=cz+dl/2d0
            
            n_points=n_points+1
            node_list(n_points,1)=cx+dl/2d0
            node_list(n_points,2)=cy-dl/2d0
            node_list(n_points,3)=cz+dl/2d0
            
          else
          
            n_faces=n_faces+1
            n_points=n_points+4
            
          end if
          
        end if

      end do ! next x cell
    end do ! next y cell
  end do ! next x cell

  if (loop.EQ.1) then
    ALLOCATE( node_list(n_points,3) )
    write(*,*)'Number of faces normal to x:',n_faces
  end if

end do ! next loop

n_points=0
do face=1,n_faces
  
! write this face as two triangles, normal to z

! traingle 1, points 1, 2, 3 on face
  write(40,'(A)')'facet normal 1.0 0.0 0.0'
  write(40,'(A)')' outer loop'
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+1,1),node_list(n_points+1,2),node_list(n_points+1,3)
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+2,1),node_list(n_points+2,2),node_list(n_points+2,3)
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+3,1),node_list(n_points+3,2),node_list(n_points+3,3)
  write(40,'(A)')'  endloop'
  write(40,'(A)')'endfacet'

! traingle 2, points 1, 4, 4 on face
  write(40,'(A)')'facet normal 1.0 0.0 0.0'
  write(40,'(A)')' outer loop'
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+1,1),node_list(n_points+1,2),node_list(n_points+1,3)
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+3,1),node_list(n_points+3,2),node_list(n_points+3,3)
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+4,1),node_list(n_points+4,2),node_list(n_points+4,3)
  write(40,'(A)')'  endloop'
  write(40,'(A)')'endfacet'
  
  n_points=n_points+4
  
end do ! next face

DEALLOCATE( node_list )

! faces normal to y

do loop=1,2

  n_faces=0
  n_points=0

  do ix=1,mesh_nx
    do iy=1,mesh_ny
      do iz=1,mesh_nz

! Look for a discontinuity in material in z and place a surface there    
        if (material_mesh(centre,ix,iy,iz).NE.material_mesh(centre,ix,iy+1,iz)) then
                  
          if (loop.EQ.2) then
          
            n_faces=n_faces+1
            
            CALL get_TLM_cell_centre_coordinate(ix,iy,iz,cx,cy,cz)
          
            n_points=n_points+1
            node_list(n_points,1)=cx-dl/2d0
            node_list(n_points,2)=cy+dl/2d0
            node_list(n_points,3)=cz-dl/2d0
            
            n_points=n_points+1
            node_list(n_points,1)=cx-dl/2d0
            node_list(n_points,2)=cy+dl/2d0
            node_list(n_points,3)=cz+dl/2d0
            
            n_points=n_points+1
            node_list(n_points,1)=cx+dl/2d0
            node_list(n_points,2)=cy+dl/2d0
            node_list(n_points,3)=cz+dl/2d0
            
            n_points=n_points+1
            node_list(n_points,1)=cx+dl/2d0
            node_list(n_points,2)=cy+dl/2d0
            node_list(n_points,3)=cz-dl/2d0
            
          else
          
            n_faces=n_faces+1
            n_points=n_points+4
            
          end if
          
        end if

      end do ! next x cell
    end do ! next y cell
  end do ! next x cell

  if (loop.EQ.1) then
    ALLOCATE( node_list(n_points,3) )
    write(*,*)'Number of faces normal to y:',n_faces
  end if

end do ! next loop

n_points=0
do face=1,n_faces
  
! write this face as two triangles, normal to z

! traingle 1, points 1, 2, 3 on face
  write(40,'(A)')'facet normal 0.0 1.0 0.0'
  write(40,'(A)')' outer loop'
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+1,1),node_list(n_points+1,2),node_list(n_points+1,3)
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+2,1),node_list(n_points+2,2),node_list(n_points+2,3)
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+3,1),node_list(n_points+3,2),node_list(n_points+3,3)
  write(40,'(A)')'  endloop'
  write(40,'(A)')'endfacet'

! traingle 2, points 1, 4, 4 on face
  write(40,'(A)')'facet normal 0.0 1.0 0.0'
  write(40,'(A)')' outer loop'
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+1,1),node_list(n_points+1,2),node_list(n_points+1,3)
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+3,1),node_list(n_points+3,2),node_list(n_points+3,3)
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+4,1),node_list(n_points+4,2),node_list(n_points+4,3)
  write(40,'(A)')'  endloop'
  write(40,'(A)')'endfacet'
  
  n_points=n_points+4
  
end do ! next face

DEALLOCATE( node_list )

! faces normal to z

do loop=1,2

  n_faces=0
  n_points=0

  do ix=1,mesh_nx
    do iy=1,mesh_ny
      do iz=1,mesh_nz

! Look for a discontinuity in material in z and place a surface there    
        if (material_mesh(centre,ix,iy,iz).NE.material_mesh(centre,ix,iy,iz+1)) then
                  
          if (loop.EQ.2) then
          
            n_faces=n_faces+1
            
            CALL get_TLM_cell_centre_coordinate(ix,iy,iz,cx,cy,cz)
          
            n_points=n_points+1
            node_list(n_points,1)=cx-dl/2d0
            node_list(n_points,2)=cy-dl/2d0
            node_list(n_points,3)=cz+dl/2d0
            
            n_points=n_points+1
            node_list(n_points,1)=cx+dl/2d0
            node_list(n_points,2)=cy-dl/2d0
            node_list(n_points,3)=cz+dl/2d0
            
            n_points=n_points+1
            node_list(n_points,1)=cx+dl/2d0
            node_list(n_points,2)=cy+dl/2d0
            node_list(n_points,3)=cz+dl/2d0
            
            n_points=n_points+1
            node_list(n_points,1)=cx-dl/2d0
            node_list(n_points,2)=cy+dl/2d0
            node_list(n_points,3)=cz+dl/2d0
            
          else
          
            n_faces=n_faces+1
            n_points=n_points+4
            
          end if
          
        end if

      end do ! next x cell
    end do ! next y cell
  end do ! next x cell

  if (loop.EQ.1) then
    ALLOCATE( node_list(n_points,3) )
    write(*,*)'Number of faces normal to z:',n_faces
  end if

end do ! next loop

n_points=0
do face=1,n_faces
  
! write this face as two triangles, normal to z

! triangle 1, points 1, 2, 3 on face
  write(40,'(A)')'facet normal 0.0 0.0 1.0'
  write(40,'(A)')' outer loop'
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+1,1),node_list(n_points+1,2),node_list(n_points+1,3)
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+2,1),node_list(n_points+2,2),node_list(n_points+2,3)
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+3,1),node_list(n_points+3,2),node_list(n_points+3,3)
  write(40,'(A)')'  endloop'
  write(40,'(A)')'endfacet'

! traingle 2, points 1, 4, 4 on face
  write(40,'(A)')'facet normal 0.0 0.0 1.0'
  write(40,'(A)')' outer loop'
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+1,1),node_list(n_points+1,2),node_list(n_points+1,3)
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+3,1),node_list(n_points+3,2),node_list(n_points+3,3)
  write(40,'(A13,3ES16.6)')'  vertex ',node_list(n_points+4,1),node_list(n_points+4,2),node_list(n_points+4,3)
  write(40,'(A)')'  endloop'
  write(40,'(A)')'endfacet'
  
  n_points=n_points+4
  
end do ! next face

DEALLOCATE( node_list )

write(40,'(A)')'endsolid'

close(unit=40) 

write(*,*)'FINISHED: mesh_to_stl'

RETURN  
  
END SUBROUTINE mesh_to_stl
!
! ______________________________________________________________________
!
!
SUBROUTINE get_TLM_cell_centre_coordinate(ix,iy,iz,cx,cy,cz)

USE PCB_simulation

IMPLICIT NONE

integer :: ix,iy,iz
real*8  :: cx,cy,cz

! START

cx=xmin + ix*dl - dl/2d0
cy=ymin + iy*dl - dl/2d0
cz=zmin + iz*dl - dl/2d0

END SUBROUTINE get_TLM_cell_centre_coordinate
!
! ______________________________________________________________________
!
!
SUBROUTINE get_TLM_cell_from_coordinate(cx,cy,cz,ix,iy,iz)

USE PCB_simulation

IMPLICIT NONE

real*8  :: cx,cy,cz
integer :: ix,iy,iz

! START

ix=INT((cx-xmin)/dl)+1
iy=INT((cy-ymin)/dl)+1
iz=INT((cz-zmin)/dl)+1

END SUBROUTINE get_TLM_cell_from_coordinate
