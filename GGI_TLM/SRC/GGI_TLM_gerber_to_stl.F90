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
!     PROGRAM GGI_TLM_gerber_to_stl
!
! DESCRIPTION
!     Convert a gerber format file into an stl file suitable for import into GGI_TLM as a surface
!	
!     
! COMMENTS
!     This process does not include some of the more complex aspects of the gerber format
!     so may require further development to convert some gerber files
!     correctly. The 'image' defined by the gerber file is built on a cartesian grid of specified
!     resolution then this cartesian pixelated image is converted to stl format. The resolution of
!     the cartesian grid should be higher than the discretisation of the subsequent GGI_TLM simulation
!     to try and prevent problems with re-meshing. The stl format output can be viewed with paraview.
! 
!
!
! HISTORY
!
!     started 10/12/19 CJS
!     
!

PROGRAM GGI_gerber_to_stl

USE gerber

IMPLICIT NONE

character(LEN=line_length) :: line_in
character(LEN=line_length) :: line
character(LEN=line_length) :: temp_line
character :: ch
character :: endchar
integer   :: first_search_char
integer :: ll,le,ll2

real*8 :: x,y,r,t,xi,yj     ! coordinates
real*8 :: xm,ym
logical :: new_coordinate

integer :: pos
integer :: rloop
integer :: count

integer ix,iy
integer iax,iay
integer :: node
integer :: triangle

! coordinate format setting default values

integer :: xfni=2
integer :: xfnd=6
integer :: yfni=2
integer :: yfnd=6

! Command line arguments
integer :: i
character(len=256) :: ipfilename
character(len=256) :: arg
logical :: DL_SET

! output filename
character(len=256) :: opfilename

! START

  CALL write_progress('STARTED: GGI_TLM_gerber_to_stl')

  i = 1
  CALL get_command_argument(i, ipfilename)
  if (len_trim(ipfilename).eq.0) then
    write(*,*)'Please provide the gerber file name as the first command line argument'
    STOP
  end if

! Open gerber file

  OPEN(unit=10,file=trim(ipfilename),status='OLD',err=9000)

  i = 2
  CALL get_command_argument(i, arg)
  if (len_trim(arg).eq.0) then
    write(*,*)'No dl has been set. This will be now be set automatically...'
    dl=0d0
    DL_SET=.FALSE.
  else
  
    read(arg,*)dl
    if (dl.LE.0d0) then
      write(*,*)'Error: dl<= 0'
      STOP
    end if
    write(*,*)'Setting dl=',dl
    DL_SET=.TRUE.
    
  end if
  
  xmin= 1D30
  xmax=-1D30
  ymin= 1D30
  ymax=-1D30
  
  do rloop=1,2
  
! no Dcodes specified initially
    nDcodes=0
    Dcode(1:maxDcodes)=0
    Atype(1:maxDcodes)=''
    Aparams(1:maxDcodes,1:maxAparams)=0d0
    Asize(1:maxDcodes)=0d0
    A_AMnumber(1:maxDcodes)=0
    A_nvars(1:maxDcodes)=0
    Avars(1:maxDcodes,1:maxAparams)=''
  
! no aperture macros specified initially
    n_AM=0
    AMname(1:maxAM)=''
    total_AM_primitives=0
    AM_n_primitives(1:maxAM)=0
    AM_to_primitive_list(1:maxAM,1:maxAM_primitives)=0
    AM_primitive_number(1:maxAM_primitives)=0
    n_AM_modifiers(1:maxAM_primitives)=0
    AM_modifiers(1:maxAM_primitives,1:maxAM_modifiers)=''
  
! initial coordinate
    x=0d0
    y=0d0

! no aperture specified initially
    aperture=0
  
    interpolation_mode=linear
    
    quadrant_mode=single
 
    polarity=dark
  
! not within a region statement
    region=.FALSE.
! 
    line_in=''
    
    DO
    
! Assemble the next string to process - this may cover multiple lines...

      line=''
      new_coordinate=.TRUE.
      
! read the next block or extended command
3 CONTINUE

!      write(*,*)'Started process with line_in=',trim(line_in),':',len_trim(line_in)

      if (len_trim(line_in).EQ.0) then
        read(10,'(A)',END=1000)line_in
      end if
      ch=line_in(1:1)
      pos=1
      
! is this an extended command or a block?
      if (ch.EQ.'%') then
! this is an extended command so read until the end of the command i.e. '%'  
        endchar='%'    
      else
! this is a block so read until the end of the block i.e. '*'  
        endchar='*'         
      end if
      first_search_char=2   ! when we look for the end char, start looking from character 2
      
5     CONTINUE
      
      ll=len_trim(line_in)
      
!      write(*,*)line_in(2:ll),' looking for:',endchar,' index=',index(line_in(2:ll),endchar)  
          
      if(index(line_in(first_search_char:ll),endchar).EQ.0) then
      
! end character not found so add this line to the string and read the next line 
        line(pos:pos+ll-1)=line_in(1:ll)
        pos=pos+ll
        read(10,'(A)',END=1000)line_in
        first_search_char=1   ! when we look for the end char, start looking from character 1 as we have read a new line
        GOTO 5
        
      else
! the end character is in this line so break the line there      
        le=index(line_in(2:ll),endchar)+1
! add the part of the line up to the end character to the string
        line(pos:pos+le-1)=line_in(1:le)
! save the rest of the string in line_in
        ll2=ll-le
        if (ll2.GT.0) then
          temp_line=line_in(le+1:ll)
          line_in=temp_line
        else
          line_in=''
        end if
        
!        write(*,*)'exit with line_in=',trim(line_in),':',len_trim(line_in)
        
      end if      
      
! At this point line should contain either the complete extended command or a single block
      
      write(*,*)'Processing line:',trim(line)
      
! Look for specific extended commands
      if (line(1:3).EQ.'%FS') then
        
        CALL read_format(line,xfni,xfnd,yfni,yfnd)
    
      else if (line(1:3).EQ.'%MO') then
      
        CALL read_scale(line,s)
    
      else if (line(1:3).EQ.'%AD') then
      
        CALL read_Dcode(line)
    
      else if (line(1:6).EQ.'%LPD*%') then
       
        polarity=dark
        write(*,*)'Setting polarity=dark'
    
      else if (line(1:6).EQ.'%LPC*%') then
      
        polarity=clear
        write(*,*)'Setting polarity=clear'        
    
      else if (line(1:3).EQ.'%AM') then
      
        CALL read_Aperture_Macro(line,n_AM)
    
      else if (line(1:3).EQ.'%AB') then
      
        ABflag=.TRUE.
    
      else if (line(1:3).EQ.'%LM') then
      
        LMflag=.TRUE.
    
      else if (line(1:3).EQ.'%LR') then
      
        LRflag=.TRUE.
    
      else if (line(1:3).EQ.'%LS') then
      
        LSflag=.TRUE.
    
      else if (line(1:3).EQ.'%SR') then
      
        SRflag=.TRUE.
    
      else if (line(1:3).EQ.'%TF') then
      
        TFflag=.TRUE.
    
      else if (line(1:3).EQ.'%TA') then
      
        TAflag=.TRUE.
    
      else if (line(1:3).EQ.'%TO') then
      
        TOflag=.TRUE.
    
      else if (line(1:3).EQ.'%TD') then
      
        TDflag=.TRUE.

      else
! process the line character by character
      
        pos=1
      
10      CONTINUE

        write(*,*)'Processing line from pos=',pos
    
        if (line(pos:pos+2).EQ.'G01') then
        
          pos=pos+3
          interpolation_mode=linear
          write(*,*)'Setting interpolation mode=linear'
              
        else if (line(pos:pos+2).EQ.'G02') then
        
          pos=pos+3
          interpolation_mode=clockwise
          write(*,*)'Setting interpolation mode=clockwise'
    
        else if (line(pos:pos+2).EQ.'G03') then
        
          pos=pos+3
          interpolation_mode=counterclockwise
          write(*,*)'Setting interpolation mode=linear'
    
        else if (line(pos:pos+2).EQ.'G36') then
        
          pos=pos+3
          region=.TRUE.
          write(*,*)'Setting interpolation region=TRUE'
          reg(1:nx,1:ny)=0       ! reset the region filling array
    
        else if (line(pos:pos+2).EQ.'G37') then
        
          if (region) then
! If a contour has been specified then go through the fill procedure and then transfer the region array to the pixel array
            CALL set_region(reg,p,nx,ny,polarity)
          end if
       
          pos=pos+3
          region=.FALSE.
          write(*,*)'Setting interpolation region=FALSE'
    
        else if (line(pos:pos+2).EQ.'G74') then
        
          pos=pos+3
          quadrant_mode=single
          write(*,*)'Setting interpolation quadrant_mode=single'
    
        else if (line(pos:pos+2).EQ.'G75') then
 
          pos=pos+3
          quadrant_mode=multi
          write(*,*)'Setting interpolation quadrant_mode=multi'
    
        else if (line(pos:pos+2).EQ.'D01') then
 
! Create line segment      
  
          pos=pos+3
          write(*,*)'Sweep from: ',xm,ym
          write(*,*)'       to : ',x,y
          write(*,*)'  xi, yj  : ',xi,yj
 
          if (rloop.EQ.2) then      
            if ( (.NOT.region).AND.(aperture.EQ.0) ) GOTO 9010  
            if (region) then
! place the line in the region filling array
              CALL sweep_line(reg,nx,ny,ap,anx,any,x,y,xm,ym,xi,yj,xmin,ymin,dl,interpolation_mode,quadrant_mode,polarity,region)
            else
              CALL sweep_line(p,nx,ny,ap,anx,any,x,y,xm,ym,xi,yj,xmin,ymin,dl,interpolation_mode,quadrant_mode,polarity,region)           
            end if
          end if
        
        else if (line(pos:pos+2).EQ.'D02') then
        
          if (region) then
! If a contour has been specified then go through the fill procedure and then transfer the region array to the pixel array
          
          end if
                
          pos=pos+3
          if (rloop.EQ.1) then        
            xmin=min(xmin,x)
            xmax=max(xmax,x)
            ymin=min(ymin,y)
            ymax=max(ymax,y)        
          end if
          
        else if (line(pos:pos+2).EQ.'D03') then
        
 ! Flash the aperture

          pos=pos+3
          write(*,*)'Flash aperture'  
          ipx=NINT((x-xmin)/dl)
          ipy=NINT((y-ymin)/dl)
          write(*,*)'Centre',ipx,ipy,' of ',nx,ny
          write(*,*)'size of Aperture:',anx,any

          if (rloop.EQ.2) then        
            if(.NOT.region) then
              if (aperture.EQ.0) GOTO 9010  
              CALL flash(p,nx,ny,ap,anx,any,x,y,xmin,ymin,dl,polarity)       
            end if
          end if
          
        else if (line(pos:pos).EQ.'D') then
      
          pos=pos+1
          
          CALL set_aperture(line,pos,nDcodes,maxDcodes,Dcode,aperture)
          
          if (Atype(aperture).NE.'M') then
            CALL set_aperture_pattern(rloop)
          else
            CALL set_aperture_macro_pattern(rloop)
          end if
        
        else if (line(pos:pos).EQ.'X') then

          pos=pos+1
          if (new_coordinate) then
            xm=x
            ym=y
            xi=0d0
            yj=0d0
            new_coordinate=.FALSE.
          end if
          CALL read_formatted_coordinate_from_line(line,pos,x,xfni,xfnd,s)
          if (rloop.EQ.1) then        
            xmin=min(xmin,x)
            xmax=max(xmax,x)
          end if
        
        else if (line(pos:pos).EQ.'Y') then

          pos=pos+1
          if (new_coordinate) then
            xm=x
            ym=y
            xi=0d0
            yj=0d0
            new_coordinate=.FALSE.
          end if
          CALL read_formatted_coordinate_from_line(line,pos,y,yfni,yfnd,s)
          if (rloop.EQ.1) then        
            ymin=min(ymin,y)
            ymax=max(ymax,y)
          end if
        
        else if (line(pos:pos).EQ.'I') then

          pos=pos+1
          if (new_coordinate) then
            write(*,*)'Error: reading centre offset before curve end point'
            STOP
          end if
          CALL read_formatted_coordinate_from_line(line,pos,xi,xfni,xfnd,s)
        
        else if (line(pos:pos).EQ.'J') then

          pos=pos+1
          if (new_coordinate) then
            write(*,*)'Error: reading centre offset before curve end point'
            STOP
          end if
          CALL read_formatted_coordinate_from_line(line,pos,yj,yfni,yfnd,s)
    
        else if (line(pos:pos+3).EQ.'M02*') then
        
          EXIT ! end of file
          
        else if (line(pos:pos).EQ.'*') then

          write(*,*)'Found end of block'
! end of the block so go to the next line
          GOTO 100
          
        else if (line(pos:pos).NE.'*') then

! Don't know how to interpret the next character so ignore this line and read the next one
          GOTO 100
                
        end if

! we get here if we have completed one of the actions and we have not reached '*' in the line        
        GOTO 10
        
      end if
      
100   CONTINUE
   
    END DO  ! go to the next line

1000 CONTINUE

    if (rloop.EQ.1) then
    
      max_Asize=0d0
      do i=1,nDcodes
        max_Asize=max(Asize(i),max_Asize)
      end do
  
      write(*,*)'Xmin=',xmin,' Xmax=',xmax
      write(*,*)'Ymin=',ymin,' Ymax=',ymax    
      write(*,*)'Max aperture size=',max_Asize
      dframe=max_Asize+abs(xmax-xmin)*0.05+abs(ymax-ymin)*0.05
      xmin=xmin-dframe
      xmax=xmax+dframe
      ymin=ymin-dframe
      ymax=ymax+dframe
      write(*,*)'Add a frame...'
      write(*,*)'Xmin=',xmin,' Xmax=',xmax
      write(*,*)'Ymin=',ymin,' Ymax=',ymax    
      
      if (.NOT.DL_SET) then
! set dl related to the average edge length of the whole PCB design
        dl=(abs(xmax-xmin)+abs(ymax-ymin))/2000d0
      end if
      
      nx=NINT((xmax-xmin)/dl)
      ny=NINT((ymax-ymin)/dl)
      
      write(*,*)'Nx=',nx,' Ny=',ny
      
      ALLOCATE( p(1:nx,1:ny) )
      p(1:nx,1:ny)=0
      ALLOCATE( reg(1:nx,1:ny) )
      reg(1:nx,1:ny)=0
  
      rewind(unit=10)
      
    end if

  end do ! next reading loop

  CLOSE(unit=10)

! We now have all the filled pixels set to '1' in the p(ix,iy) array
! Now turn this into a triangular mesh

  CALL triangulate_surface()
  
! write stl format data (can be imported into GGI_TLM)
  opfilename=trim(ipfilename)//'.stl'
  OPEN(unit=20,file=trim(opfilename))
 
  write(20,'(A)')'solid ascii'
  
  do triangle=1,n_triangles
    
      ! loop over the nodes of this triangle, going back to the first to close the triangle
    write(20,'(A)')'facet normal 0.0 0.0 1.0'
    write(20,'(A)')'  outer loop'
    do i=1,3
    
      node=triangle_to_node_list(triangle,i)      
      write(20,'(A13,3ES16.6)')'    vertex   ',node_list(node,1),node_list(node,2),0d0
      
    end do ! next node

    write(20,'(A)')'  endloop'
    write(20,'(A)')'endfacet'

  end do ! next triangle
  write(20,'(A)')'endsolid'
  
  CLOSE(unit=20)
  
  DEALLOCATE( p )
  DEALLOCATE( reg )
  if (allocated( ap )) Deallocate( ap )
  
  if (allocated(node_list)) DEALLOCATE( node_list )
  if (allocated(triangle_to_node_list)) DEALLOCATE( triangle_to_node_list )
  
! Write warnings for unprocessed commands:

  write(*,*)
  write(*,*)'Warnings:'
  write(*,*)

  if (AMflag) then
  
    write(*,*)'WARNING: AM command not processed (Aperture Macro)'
    if (AMflag_outline) write(*,*)'Outline primitive not processes'
    if (AMflag_Moire) write(*,*)'Outline Moire not processes'
    if (AMflag_thermal) write(*,*)'Outline thermal not processes'

  end if

  if (ABflag) then
    write(*,*)'WARNING: AB command not processed (Aperture Block)'
  end if

  if (LMflag) then
    write(*,*)'WARNING: LM command not processed (Load Mirror)'
  end if

  if (LRflag) then
    write(*,*)'WARNING: LR command not processed (Load Reflection)'
  end if

  if (LSflag) then
    write(*,*)'WARNING: LS command not processed (Load Scale)'
  end if

  if (SRflag) then
    write(*,*)'WARNING: SR command not processed (Step and Repeat)'
  end if

  if (TFflag) then
    write(*,*)'WARNING: TF command not processed (Attribute File)'
  end if

  if (TAflag) then
    write(*,*)'WARNING: TA command not processed (Attribute Aperture)'
  end if

  if (TOflag) then
    write(*,*)'WARNING: TO command not processed (Attribute Object)'
  end if

  if (TDflag) then
    write(*,*)'WARNING: TD command not processed (Attribute Delete)'
  end if
  
  write(*,*)'_______________________________________________________________'
  write(*,*)''
  write(*,*)'stl filename is:',trim(opfilename)
  write(*,*)''
  write(*,*)'stl file generated with dl set to ',dl,'m'
  write(*,*)''
  write(*,*)'If the resolution of the stl file is not sufficient then r-run'    
  write(*,*)'with dl set as the second command line argument e.g.'    
  write(*,*)'GGI_TLM_gerber_to_stl ',trim(ipfilename),' ',dl/2d0
  write(*,*)'_______________________________________________________________'
  write(*,*)''
  
  CALL write_progress('FINISHED: GGI_TLM_gerber_to_stl')

  STOP
  
9000 write(*,*)'ERROR: Cannot open the file:',trim(arg)
  STOP
  
9010 write(*,*)'ERROR: Aperture is not defined'
  STOP

END PROGRAM GGI_gerber_to_stl
