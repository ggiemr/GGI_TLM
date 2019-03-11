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
!       SUBROUTINE  triangulate_surface
!     
!
! DESCRIPTION
!     
!     
! COMMENTS
!
!
! HISTORY
!
!     started 10/12/19 CJS
!     
!
SUBROUTINE  triangulate_surface()

USE gerber

IMPLICIT NONE

integer :: ix,iy,ixmax,iymax,ix2,iy2
real*8  :: x1,y1,x2,y2
logical :: max_x,max_y
integer :: loop

! START

  write(*,*)'CALLED triangulate_surface'
  
  do loop=1,2
  
    n_nodes=0
    n_triangles=0

! step through the pixel array looking for 1's i.e. un-triangulated areas 

    do ix=1,nx  
      do iy=1,ny
    
        if (p(ix,iy).EQ.1) then
        
! find a rectangle with corners (ix,iy),(ixmax,iymax) which is completely ones and then mesh this with two triangles
          max_x=.FALSE.
          max_y=.FALSE.
          ixmax=ix
          iymax=iy
          
          do while( (.NOT.max_x).OR.(.NOT.max_y) )
          
            if (.NOT.max_x) then

! try to increase the x extent of the rectangle
              ixmax=ixmax+1
! check the new row of cells to see if they are all set to 1

              do iy2=iy,iymax
              
                if (p(ixmax,iy2).NE.1) then
! this is one row beyond the x limit of the rectangle so set the x limit to ixmax-1 and flag that we have reached max_x
                  ixmax=ixmax-1
                  max_x=.TRUE.
                  EXIT
                end if
                
              end do
              
            end if  !.NOT.max_x
          
            if (.NOT.max_y) then

! try to increase the y extent of the rectangle
              iymax=iymax+1
! check the new row of cells to see if they are all set to 1

              do ix2=ix,ixmax
              
                if (p(ix2,iymax).NE.1) then
! this is one row beyond the y limit of the rectangle so set the y limit to iymax-1 and flag that we have reached max_y
                  iymax=iymax-1
                  max_y=.TRUE.
                  EXIT
                end if
                
              end do
              
            end if  !.NOT.max_y
          
          end do   ! (.NOT.max_x).OR.(.NOT.max_y) 
      
          if (loop.EQ.2) then
! add the two triangles defining this rectangle to the node_list and the triangle_to_node_list
            x1=xmin+ix*dl-dl
            y1=ymin+iy*dl-dl
            x2=xmin+ixmax*dl
            y2=ymin+iymax*dl
            
            n_nodes=n_nodes+1
            node_list(n_nodes,1)=x1
            node_list(n_nodes,2)=y1
            n_nodes=n_nodes+1
            node_list(n_nodes,1)=x2
            node_list(n_nodes,2)=y1
            n_nodes=n_nodes+1
            node_list(n_nodes,1)=x2
            node_list(n_nodes,2)=y2
            n_nodes=n_nodes+1
            node_list(n_nodes,1)=x1
            node_list(n_nodes,2)=y2
            
            n_triangles=n_triangles+1
            triangle_to_node_list(n_triangles,1)=n_nodes-3
            triangle_to_node_list(n_triangles,2)=n_nodes-2
            triangle_to_node_list(n_triangles,3)=n_nodes-1
            n_triangles=n_triangles+1
            triangle_to_node_list(n_triangles,1)=n_nodes-3
            triangle_to_node_list(n_triangles,2)=n_nodes-1
            triangle_to_node_list(n_triangles,3)=n_nodes

          else
! just increase the number of nodes and number of triangles
            n_nodes=n_nodes+4
            n_triangles=n_triangles+2
            
          end if

! update the rectangle of pixels         
          do ix2=ix,ixmax
            do iy2=iy,iymax
              p(ix2,iy2)=2
            end do
          end do
      
        end if ! p(ix,iy).EQ.1
      
      end do
    end do
  
    if (loop.EQ.1) then
! Allocate memory for the triangle to node list

      write(*,*)'Number of triangles=',n_triangles
      write(*,*)'Number of nodes=',n_nodes

      ALLOCATE( node_list(1:n_nodes,1:2) )
      ALLOCATE( triangle_to_node_list(1:n_triangles,1:3) )
    
! put the renumbered pixel values back as they were before
      do ix2=1,nx
        do iy2=1,ny
          if (p(ix2,iy2).EQ.1) then
            write(*,*)'Error in triangulate surface: unfilled cell',ix2,iy2
          end if
          if (p(ix2,iy2).EQ.2) p(ix2,iy2)=1
        end do
      end do
  
    end if
  
  end do

  write(*,*)'Done: triangulate_surface'

END SUBROUTINE triangulate_surface
