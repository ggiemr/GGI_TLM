! SUBROUTINE equivalent_source_model()
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
! SUBROUTINE 
!
! NAME equivalent_source_model
!    
!
! DESCRIPTION
!    Read complex frequency domain near field scan data on a uniform grid and calculate an equivalent source
!    model consisting of electric and magnetic dipoles on a defined source plane grid
!    
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/3/2024 CJS
!
SUBROUTINE equivalent_source_model()

USE constants
USE post_process
USE file_information

IMPLICIT NONE

integer	:: function_number

! local variables

  character(len=256)	:: Hx_filename
  character(len=256)	:: Hy_filename
    
  character(len=256)	:: command
  
  logical		:: file_exists
  
! Hx Near field scan plane data
  
  integer                :: Hx_n_lines
  complex*16,allocatable :: Hx_measured(:)
  
! Hy Near field scan plane data
  
  integer                :: Hy_n_lines
  complex*16,allocatable :: Hy_measured(:)

! combined data
  integer               :: n_lines
  integer               :: nx_measured,ny_measured,n_measured

  real*8,allocatable    :: x_measured(:)
  real*8,allocatable    :: y_measured(:)
  real*8,allocatable    :: z_measured(:)
  
  real*8  :: xmin_measured,xmax_measured
  real*8  :: ymin_measured,ymax_measured
  
! equivalent source plane specification

  real*8  :: xmin_source,xmax_source,dx_source
  real*8  :: ymin_source,ymax_source,dy_source
  real*8  :: z_source

  real*8,allocatable    :: x_source(:)
  real*8,allocatable    :: y_source(:)
  
  integer :: nx_source,ny_source,n_source
  
! Electric dipoles

  logical :: use_electric_dipoles

! Magnetic dipoles

  logical :: use_magnetic_dipoles

! Tangential dipoles on source plane 

  logical :: use_tangential_dipoles_only
  integer :: n_source_components

! PEC plane

  logical :: use_PEC_plane
  real*8  :: z_PEC
  
  integer :: n_unknows_per_source_point

! matrix equation 

  integer :: matdim 
  complex*16,allocatable :: LHS(:),LHS2(:),RHS(:),G(:,:),GI(:,:)

! Source calculation stuff

  real*8    :: k0,f,w
  
  complex*16:: Hx_Px,Hx_Py,Hx_Pz,Hy_Px,Hy_Py,Hy_Pz
  
  complex*16:: Hx_Mx,Hx_My,Hx_Mz,Hy_Mx,Hy_My,Hy_Mz
  
  complex*16 :: Mx,My,Mz,Px,Py,Pz
  
  complex*16 :: Hx,Hy
  
! near to far field transformation data

  logical :: calc_FF
  
  integer :: theta_loop,n_theta
  real*8  :: Theta_min,Theta_max,Theta_step
  
  integer :: phi_loop,n_phi
  real*8 :: Phi_min,Phi_max,Phi_step
  
  integer :: n_data,count
  complex*16,allocatable :: Etheta(:)
  complex*16,allocatable :: Ephi(:)
  
  real*8  	:: rp_cos_psi
  real*8  	:: r(3),rp(3)
  real*8  	:: r0,theta,phi
  
  real*8 	:: ct,cp,st,sp
  
  real*8 	:: dA
  
  complex*16 	:: Ntheta,Nphi,Ltheta,Lphi
  complex*16 	:: Nz2D,Lz2D
  complex*16 	:: prop,prop2
  complex*16 	:: Jsx,Jsy,Jsz,Msx,Msy,Msz

! general data

  real*8 :: hzm   ! height (z value) of near field scan measurement plane
  real*8 :: hzs   ! height (z value) of near field source plane

  integer :: line
  real*8  :: x,y,z,re,im
  real*8  :: xs0,ys0,zs0
  
  integer :: ix,iy
  integer :: row,col
  integer :: col_offset
  
  integer :: i_source,i_measured
  
  integer :: n_LHS,n_RHS
  integer :: Gn_rows,Gn_cols
  
  character :: ch


! START

  write(*,*)
  write(*,*)'Output files:'
  command='ls -ltr '
  CALL system(command)

! get the filename for the Hx data

  write(*,*)'Enter the filename for the Hx data on the scan plane'
  read(*,'(A256)')Hx_filename
  inquire(file=trim(Hx_filename),exist=file_exists)
  if (.NOT.file_exists) then
    GOTO 9000
  end if
  write(record_user_inputs_unit,'(A)')trim(Hx_filename)

  OPEN(unit=local_file_unit,file=Hx_filename)

  CALL write_file_format_information(local_file_unit,Hx_n_lines)
  
  rewind(unit=local_file_unit)


! get the filename for the Hy data

  write(*,*)'Enter the filename for the Hy data on the scan plane'
  read(*,'(A256)')Hy_filename
  inquire(file=trim(Hy_filename),exist=file_exists)
  if (.NOT.file_exists) then
    GOTO 9020
  end if
  write(record_user_inputs_unit,'(A)')trim(Hy_filename)

  OPEN(unit=local_file_unit2,file=Hy_filename)

  CALL write_file_format_information(local_file_unit2,Hy_n_lines)
  
  rewind(unit=local_file_unit2)

! Initial checks:
 
  if (Hx_n_lines.NE.Hy_n_lines) then
    CALL write_line('Error: different number of points in Hx and Hy data files ',0,.TRUE.)
    write(*,*)'Hx_n_lines=',Hx_n_lines
    write(*,*)'Hy_n_lines=',Hy_n_lines
    STOP 1
  end if 

! Read the Hx nad Hy scan data

  n_lines=Hx_n_lines

  ALLOCATE( x_measured(n_lines) )
  ALLOCATE( y_measured(n_lines) )
  ALLOCATE( z_measured(n_lines) )
  ALLOCATE( Hx_measured(n_lines) )
  ALLOCATE( Hy_measured(n_lines) )
  
  do line=1,n_lines

! Read Hx data  
    read(local_file_unit,*,ERR=9010,END=9010)x,y,z,re,im
    
    x_measured(line)=x
    y_measured(line)=y
    z_measured(line)=z
    Hx_measured(line)=dcmplx(re,im)

! check data is on a z plane    
    if (line.Eq.1) then
      hzm=z
    else
      if (z_measured(line).NE.hzm) then
        write(*,*)'Error: Hx measured data is not on a plane normal to z. z.NE.hzm on line',line
        STOP 1
      end if
    end if

! Read Hy data  
    read(local_file_unit2,*,ERR=9030,END=9030)x,y,z,re,im
    
    Hy_measured(line)=dcmplx(re,im)

! check Hx and Hy are at the same point
    if (    (x.NE.x_measured(line))       &
        .OR.(y.NE.y_measured(line))       &
        .OR.(z.NE.z_measured(line)) ) then
    
      write(*,*)'Error: x,y,z position discrepancy on line',line
      write(*,*)'Hx position:',real(x_measured(line)),real(y_measured(line)),real(z_measured(line))
      write(*,*)'Hy position:',real(x),real(y),real(z)
      STOP 1   
    end if
    
  end do ! next line to read

  write(*,*)n_lines,' lines of data read successfully'
  
  close(unit=local_file_unit)
  close(unit=local_file_unit2)

  n_measured=n_lines
  
! Work out the dimensions of the measured scan plane
  xmin_measured= 1e30
  xmax_measured=-1e30
  ymin_measured= 1e30
  ymax_measured=-1e30
  
  do line=1,n_lines
    xmin_measured=min(xmin_measured,x_measured(line))
    xmax_measured=max(xmax_measured,x_measured(line))
    ymin_measured=min(ymin_measured,y_measured(line))
    ymax_measured=max(ymax_measured,y_measured(line))
  end do
  
  write(*,*)'xmin_measured=',xmin_measured
  write(*,*)'xmax_measured=',xmax_measured
  write(*,*)'ymin_measured=',ymin_measured
  write(*,*)'ymax_measured=',ymax_measured
  write(*,*)'z_measured=',hzm

! read the analysis frequency
  write(*,*)
  write(*,*)'Enter the frequency at which the H fields have been measured (Hz)'
  read(*,*)f
  write(record_user_inputs_unit,*)real(f),'  # measurement frequency'

! Read the equivalent source grid specification

  write(*,*)
  write(*,*)'Enter the source plane xmin and xmax values'
  read(*,*)xmin_source,xmax_source
  write(record_user_inputs_unit,*)real(xmin_source),real(xmax_source),'  # source plane xmin and xmax values'

  write(*,*)'Enter the source plane number of values in x'
  read(*,*)nx_source
  write(record_user_inputs_unit,*)nx_source,'  # nx_source '

  write(*,*)'Enter the source plane ymin and ymax values'
  read(*,*)ymin_source,ymax_source
  write(record_user_inputs_unit,*)real(xmin_source),real(xmax_source),'  # source plane ymin and ymax values'

  write(*,*)'Enter the source plane number of values in y'
  read(*,*)ny_source
  write(record_user_inputs_unit,*)ny_source,'  # ny_source '

  write(*,*)'Enter the source plane z value'
  read(*,*)z_source
  write(record_user_inputs_unit,*)real(z_source),'  # source plane z value'

  n_source=nx_source*ny_source

! set the x and y positions on the source plane grid

  ALLOCATE( x_source(n_source) )
  ALLOCATE( y_source(n_source) )

  if (nx_source.Eq.1) then
    dx_source=0d0
  else
    dx_source=(xmax_source-xmin_source)/dble(nx_source-1)
  end if

  if (ny_source.Eq.1) then
    dy_source=0d0
  else
    dy_source=(ymax_source-ymin_source)/dble(ny_source-1)
  end if
  
  count=0
  
  do ix=1,nx_source
  
    do iy=1,ny_source
    
      count=count+1
    
      if (nx_source.EQ.1) then
        x_source(count)=(xmin_source+xmax_source)/2d0
      else
        x_source(count)=xmin_source+dx_source*dble(ix-1)
      end if
    
      if (ny_source.EQ.1) then
        y_source(count)=(ymin_source+ymax_source)/2d0
      else
        y_source(count)=ymin_source+dy_source*dble(iy-1)
      end if  
    
    end do
    
  end do
  
  write(*,*)''
  write(*,*)'Total number of source plane points = ',n_source

  write(*,*)'Do you want to use electric dipoles on the source plane (y/n)'
  read(*,'(A1)')ch
  if ( (ch.Eq.'y').OR.(ch.Eq.'Y') ) then
    use_electric_dipoles=.TRUE.
  else
    use_electric_dipoles=.FALSE.
  end if  
  write(record_user_inputs_unit,'(A1,A)')ch,'   #use electric dipoles on the source plane (y/n)'

  
  write(*,*)'Do you want to use magnetic dipoles on the source plane (y/n)'
  read(*,'(A1)')ch
  if ( (ch.Eq.'y').OR.(ch.Eq.'Y') ) then
    use_magnetic_dipoles=.TRUE.
  else
    use_magnetic_dipoles=.FALSE.
  end if  
  write(record_user_inputs_unit,'(A1,A)')ch,'   #use magnetic dipoles on the source plane (y/n)'
  
  if ( (.NOT.use_electric_dipoles).AND.(.NOT.use_magnetic_dipoles) ) then
    write(*,*)'Error: No electric or magnetic dipoles on equivalent source plane...'
    STOP 1   
  end if
    
  write(*,*)'Do you want to use only tangential dipole source components on the source plane (y/n)'
  read(*,'(A1)')ch
  if ( (ch.Eq.'y').OR.(ch.Eq.'Y') ) then
    use_tangential_dipoles_only=.TRUE.
    n_source_components=2
  else
    use_tangential_dipoles_only=.FALSE.
    n_source_components=3
  end if  
  write(record_user_inputs_unit,'(A1,A)')ch,'   #use only tangential dipole source components on the source plane (y/n)'

  write(*,*)'Do you want to include image source model of a conducting ground plane (y/n)'
  read(*,'(A1)')ch
  if ( (ch.Eq.'y').OR.(ch.Eq.'Y') ) then
    use_PEC_plane=.TRUE.
  else
    use_PEC_plane=.FALSE.
  end if  
  write(record_user_inputs_unit,'(A1,A)')ch,'   #use image source model of a conducting ground plane (y/n)'
 
  if (use_PEC_plane) then
    write(*,*)'Enter the z coordinate of the conducting ground plane'
    read(*,*)z_PEC
    write(record_user_inputs_unit,*)real(z_PEC),'  # z_PEC '
  end if

! Work out the dimensions of the matrix system

  n_LHS=2*n_measured              ! factor of 2 for the two field components
  n_unknows_per_source_point=0
  
  if (use_electric_dipoles) then
    n_RHS=n_source_components*n_source              ! 
    n_unknows_per_source_point=n_source_components
  else
    n_RHS=0
  end if
  
  if (use_magnetic_dipoles) then
    n_RHS=n_RHS+n_source_components*n_source        ! 
    n_unknows_per_source_point=n_unknows_per_source_point+n_source_components
  end if
  
  write(*,*)
  if (use_electric_dipoles) write(*,*)'Use electric Dipoles'
  if (use_magnetic_dipoles) write(*,*)'Use magnetic Dipoles'
  if (use_tangential_dipoles_only) then
    write(*,*)'Use tangential Dipole components only'
  else
    write(*,*)'Use all three Dipole components '
  end if
  if (use_PEC_plane) write(*,*)'Use PEC plane, z_PEC=',real(z_PEC)
  
  write(*,*)'Number of known field components on LHS          =',n_LHS 
  write(*,*)'Number of unknown source dipole components on RHS=',n_RHS 
  write(*,*)'Number of unknowns on each source point=',n_unknows_per_source_point
  
  if (n_RHS.GT.n_LHS) then
    write(*,*)'ERROR: number of unknowns is greater than the number of measured field samples'
    STOP 1
  end if
  
! Allocate and reset matrix equation data

  matdim=n_LHS

  ALLOCATE( LHS(matdim) )
  ALLOCATE( LHS2(matdim) )
  ALLOCATE( RHS(matdim) )
  ALLOCATE( G(matdim,matdim) )
  ALLOCATE( GI(matdim,matdim) )
  
  LHS(:)=(0d0,0d0)
  LHS2(:)=(0d0,0d0)
  RHS(:)=(0d0,0d0)
  G(:,:)=(0d0,0d0)
  GI(:,:)=(0d0,0d0)
  
! loop over the measured data points setting the LHS vector 
  row=0
  do i_measured=1,n_measured
    LHS(row+1)=Hx_measured(i_measured)
    LHS(row+2)=Hy_measured(i_measured)
    row=row+2
  end do ! next measured point
  
! loop over the source plane setting the G matrix elements

  w=2d0*pi*f
  k0=w/c0

  col=0
  do i_source=1,n_source

! Coordinates of source dipole point  
    xs0=x_source(i_source)
    ys0=y_source(i_source)
    zs0=z_source
    
    row=0
    do i_measured=1,n_measured

! Coordinates of measured data point  
      x=x_measured(i_measured)
      y=y_measured(i_measured)
      z=hzm
      
      col_offset=0
      
      if (use_electric_dipoles) then
      
        CALL electric_dipole_source(xs0,ys0,zs0,x,y,z,k0,Hx_Px,Hx_Py,Hx_Pz,Hy_Px,Hy_Py,Hy_Pz,use_PEC_plane,z_PEC)
      
        col_offset=col_offset+1          
! contribution to Hx from Px  
        G(row+1,col+col_offset)=Hx_Px
      
! contribution to Hy from Px  
        G(row+2,col+col_offset)=Hy_Px
      
        col_offset=col_offset+1        
! contribution to Hx from Py  
        G(row+1,col+col_offset)=Hx_Py
      
! contribution to Hy from Py 
        G(row+2,col+col_offset)=Hy_Py
        
        if (.NOT.use_tangential_dipoles_only) then
        
          col_offset=col_offset+1        
! contribution to Hx from Pz  
          G(row+1,col+col_offset)=Hx_Pz
      
! contribution to Hy from Pz 
          G(row+2,col+col_offset)=Hy_Pz
        
        end if  ! use z component of source dipoles
              
      end if ! contributions from electric dipoles
            
      if (use_magnetic_dipoles) then
      
        CALL magnetic_dipole_source(xs0,ys0,zs0,x,y,z,k0,Hx_Mx,Hx_My,Hx_Mz,Hy_Mx,Hy_My,Hy_Mz,use_PEC_plane,z_PEC)
      
        col_offset=col_offset+1          
! contribution to Hx from Px  
        G(row+1,col+col_offset)=Hx_Mx
      
! contribution to Hy from Px  
        G(row+2,col+col_offset)=Hy_Mx
      
        col_offset=col_offset+1        
! contribution to Hx from Py  
        G(row+1,col+col_offset)=Hx_My
      
! contribution to Hy from Py 
        G(row+2,col+col_offset)=Hy_My
        
        if (.NOT.use_tangential_dipoles_only) then
        
          col_offset=col_offset+1        
! contribution to Hx from Pz  
          G(row+1,col+col_offset)=Hx_Mz
      
! contribution to Hy from Pz 
          G(row+2,col+col_offset)=Hy_Mz
        
        end if  ! use z component of source dipoles
              
      end if ! contributions from electric dipoles
      
      row=row+2
      
    end do ! next measured point
    
    col=col+n_unknows_per_source_point
    
  end do ! next source dipole point

!! Write g matrix to screen  
!  do ix=1,n_LHS
!    do iy=1,n_RHS
!      write(*,*)ix,iy,cmplx(G(ix,iy))
!    end do
!  end do
  
! Solve the system of equations using Moore_Penrose inverse

  CALL cinvert_Moore_Penrose(G,n_LHS,n_RHS,GI,matdim)
  
  CALL cmatvmul(GI,n_RHS,n_LHS,LHS,n_LHS,RHS,matdim)
  
! Now we have the RHS vector of source dipoles, calcuulate LHS2=[G](RHS)
! to determine the field on the measurement plane prediicted by these sources

  CALL cmatvmul(G,n_LHS,n_RHS,RHS,n_RHS,LHS2,matdim)
  
! Write the equivalent source field on measurement plane to file

  open(unit=local_file_unit , file='Hx_from_equiv_sources.dat')
  open(unit=local_file_unit2, file='Hy_from_equiv_sources.dat')
  
  row=0
  do i_measured=1,n_measured

! Coordinates of measured data point  
    x=x_measured(i_measured)
    y=y_measured(i_measured)
    z=hzm
    
    Hx=LHS2(row+1)
    Hy=LHS2(row+2)
    
    write(local_file_unit ,'(6ES12.4)')x,y,z,real(Hx),imag(Hx),abs(Hx)
    write(local_file_unit2,'(6ES12.4)')x,y,z,real(Hy),imag(Hy),abs(Hy)    
    
    row=row+2

  end do ! next measured data point
  
  close(unit=local_file_unit)
  close(unit=local_file_unit2)

! Extract the equivalent sources and write to file

  if (use_electric_dipoles) then
    write(*,*)'Writing electric dipole data to file: Electric_dipoles.dat'
    open(unit=local_file_unit , file='Electric_dipoles.dat')
  end if
  
  if (use_magnetic_dipoles) then
    write(*,*)'Writing magnetic dipole data to file: Magnetic_dipoles.dat'
    open(unit=local_file_unit2, file='Magnetic_dipoles.dat')
  end if
  
  col=0
  do i_source=1,n_source
  
! Coordinates of source dipole point  
    xs0=x_source(i_source)
    ys0=y_source(i_source)
    zs0=z_source
      
    col_offset=0
      
    if (use_electric_dipoles) then

      col_offset=col_offset+1          
      Px=RHS(col+col_offset)
      
      col_offset=col_offset+1        
      Py=RHS(col+col_offset)
              
      if (.NOT.use_tangential_dipoles_only) then       
        col_offset=col_offset+1        
        Pz=RHS(col+col_offset)      
      else
        Pz=(0d0,0d0)        
      end if  ! use z component of source dipoles
       
!      write(*,*)'Electric Dipoles for Source point ',i_source
!      write(*,*)'Px=',cmplx(Px)
!      write(*,*)'Py=',cmplx(Py)
!      write(*,*)'Pz=',cmplx(Pz)

       write(local_file_unit,'(6ES12.4)')xs0,ys0,zs0,abs(Px),abs(Py),abs(Pz)
             
    end if ! contributions from electric dipoles
    
    if (use_magnetic_dipoles) then

      col_offset=col_offset+1          
      Mx=RHS(col+col_offset)
      
      col_offset=col_offset+1        
      My=RHS(col+col_offset)
              
      if (.NOT.use_tangential_dipoles_only) then        
        col_offset=col_offset+1        
        Mz=RHS(col+col_offset)
      else
        Mz=(0d0,0d0)                
      end if  ! use z component of source dipoles
      
!      write(*,*)'Magnetic Dipoles for Source point ',i_source
!      write(*,*)'Mx=',cmplx(Mx)
!      write(*,*)'My=',cmplx(My)
!      write(*,*)'Mz=',cmplx(Mz)

      write(local_file_unit2,'(6ES12.4)')xs0,ys0,zs0,abs(Mx),abs(My),abs(Mz)
      
    end if ! contributions from magnetic dipoles
    
    col=col+n_unknows_per_source_point

  end do ! next source point
  
  if (use_electric_dipoles) close(unit=local_file_unit)
  
  if (use_magnetic_dipoles) close(unit=local_file_unit2)

! Here we could write the equivalent source fields on the measurement plane or other output plane for comparison

! Here we could write the far field pattern if required

  write(*,*)'Do you want to calculate the far field pattern (y/n)'
  read(*,'(A1)')ch
  if ( (ch.Eq.'y').OR.(ch.Eq.'Y') ) then
    calc_FF=.TRUE.
  else
    calc_FF=.FALSE.
  end if  
  write(record_user_inputs_unit,'(A1,A)')ch,'   # calculate far field pattern (y/n)'

  if (calc_FF) then
  
    write(*,*)'Enter far field pattern Theta_min Theta_max and Theta_step (degrees)'
    read(*,*)Theta_min,Theta_max,Theta_step
    write(record_user_inputs_unit,'(3ES12.4,A)')Theta_min,Theta_max,Theta_step,' # Theta_min Theta_max and Theta_step (degrees)'

! convert to radians   
    Theta_min =Theta_min*pi/180d0
    Theta_max =Theta_max*pi/180d0
    Theta_step=Theta_step*pi/180d0   
   
    write(*,*)'Enter far field pattern Phi_min Phi_max and Phi_step (degrees)'
    read(*,*)Phi_min,Phi_max,Phi_step
    write(record_user_inputs_unit,'(3ES12.4,A)')Phi_min,Phi_max,Phi_step,' # Phi_min Phi_max and Phi_step (degrees)'
   
! convert to radians   
    Phi_min =Phi_min*pi/180d0
    Phi_max =Phi_max*pi/180d0
    Phi_step=Phi_step*pi/180d0   
  
    n_theta=INT( (theta_max-theta_min)/theta_step )+1
    n_phi  =INT( (phi_max-phi_min)/phi_step )+1
    
    n_data=n_theta*n_phi
    
    ALLOCATE( Etheta(1:n_data) )
    ALLOCATE( Ephi(1:n_data) )
    
    Etheta(1:n_data)      =(0d0,0d0)
    Ephi(1:n_data)        =(0d0,0d0)

! Calculate far field data
    r0=1d0
    prop2=j*k0*exp(j*k0*r0)/(4d0*pi*r0)

! loop over theta and phi
     
    count=0
 
    do theta_loop=1,n_theta
      theta= theta_min+(theta_loop-1)*theta_step

      do phi_loop=1,n_phi
        phi= phi_min+(phi_loop-1)*phi_step

! Observation point       
        r(1)=r0*sin(theta)*cos(phi)
        r(2)=r0*sin(theta)*sin(phi)
        r(3)=r0*cos(theta)
    
    	Ntheta=(0d0,0d0)
        Nphi  =(0d0,0d0)
        Ltheta=(0d0,0d0)
        Lphi  =(0d0,0d0)

        ct=cos(theta)
        st=sin(theta)
        cp=cos(phi)
        sp=sin(phi)  

! Loop over sources, adding contributions to the far field  
        col=0
        do i_source=1,n_source
  
! Coordinates of source dipole point  
          xs0=x_source(i_source)
          ys0=y_source(i_source)
          zs0=z_source
       
          rp(1)=xs0
          rp(2)=ys0
          rp(3)=zs0
          
          rp_cos_psi=r(1)*rp(1)+r(2)*rp(2)+r(3)*rp(3)
       
          col_offset=0
      
          if (use_electric_dipoles) then

            col_offset=col_offset+1          
            Px=RHS(col+col_offset)
      
            col_offset=col_offset+1        
            Py=RHS(col+col_offset)
              
            if (.NOT.use_tangential_dipoles_only) then       
              col_offset=col_offset+1        
              Pz=RHS(col+col_offset)      
            else
              Pz=(0d0,0d0)        
            end if  ! use z component of source dipoles
            
          else
            Px=(0d0,0d0)
            Py=(0d0,0d0)
            Pz=(0d0,0d0)  
          end if ! contributions from electric dipoles
    
          if (use_magnetic_dipoles) then

            col_offset=col_offset+1          
            Mx=RHS(col+col_offset)
      
            col_offset=col_offset+1        
            My=RHS(col+col_offset)
              
            if (.NOT.use_tangential_dipoles_only) then        
              col_offset=col_offset+1        
              Mz=RHS(col+col_offset)
            else
              Mz=(0d0,0d0)                
            end if  ! use z component of source dipoles
          else
            Mx=(0d0,0d0)
            My=(0d0,0d0)
            Mz=(0d0,0d0)              
          end if ! contributions from magnetic dipoles
    
    	  prop=exp(j*k0*rp_cos_psi)
    
    	  Ntheta=Ntheta+(Px*ct*cp+Py*ct*sp-Pz*st)*prop
    	  Nphi=Nphi+(-Px*sp+Py*cp)*prop
        
    	  Ltheta=Ltheta+(Mx*ct*cp+My*ct*sp-Mz*st)*prop
    	  Lphi=Lphi+(-Mx*sp+My*cp)*prop
          
! add image sources if present
          if (use_PEC_plane) then
  
! Coordinates of image source dipole point  
            xs0=x_source(i_source)
            ys0=y_source(i_source)
            zs0=2d0*z_PEC-z_source
       
            rp(1)=xs0
            rp(2)=ys0
            rp(3)=zs0
          
            rp_cos_psi=r(1)*rp(1)+r(2)*rp(2)+r(3)*rp(3)
              
! Calculate image sources
            Px=-Px
            Py=-Py
            Mz=-Mz
    
    	    prop=exp(j*k0*rp_cos_psi)
    
    	    Ntheta=Ntheta+(Px*ct*cp+Py*ct*sp-Pz*st)*prop
    	    Nphi=Nphi+(-Px*sp+Py*cp)*prop
        
    	    Ltheta=Ltheta+(Mx*ct*cp+My*ct*sp-Mz*st)*prop
    	    Lphi=Lphi+(-Mx*sp+My*cp)*prop
          
          end if ! use_PEC_plane
          
          col=col+n_unknows_per_source_point

        end do ! next source point
      
        count=count+1
                                    
        Etheta(count)=-prop2*(Lphi+Z0*Ntheta)
        Ephi(count)  = prop2*(Ltheta-Z0*Nphi)
        
      end do ! next phi

    end do ! next theta


! Write far field data to file
 
    OPEN(unit=far_field_output_unit,file='Far_field_from_equivalent_source.fout')
        
    write(far_field_output_unit,8010)'# ',n_theta,n_phi
8010  format(A2,2I6)

! loop over theta and phi
     
    count=0
 
    do theta_loop=1,n_theta
      theta= theta_min+(theta_loop-1)*theta_step

      do phi_loop=1,n_phi
        phi= phi_min+(phi_loop-1)*phi_step
      
        count=count+1
                                    
        write(far_field_output_unit,8020)theta*180.0/pi,phi*180.0/pi,        &
                                       abs(Etheta(count)),abs(Ephi(count))
8020    format(4E14.6)

      end do ! next phi

    end do ! next theta

    
    DEALLOCATE( Etheta )
    DEALLOCATE( Ephi )
 
  end if ! calc_FF


  DEALLOCATE( x_measured )
  DEALLOCATE( y_measured )
  DEALLOCATE( z_measured )
  DEALLOCATE( Hx_measured )
  DEALLOCATE( Hy_measured )

  DEALLOCATE( x_source )
  DEALLOCATE( y_source )

  DEALLOCATE( LHS )
  DEALLOCATE( LHS2 )
  DEALLOCATE( RHS )
  DEALLOCATE( G )
  DEALLOCATE( GI )

  RETURN

     
9000 CALL write_line('Error opening Hx data file',0,.TRUE.)
     CALL write_line('filename:',0,.TRUE.)
     write(*,*)trim(Hx_filename)
     STOP 1
      
9010 CALL write_line('Error reading Hx data file',0,.TRUE.)
     CALL write_line('filename:',0,.TRUE.)
     write(*,*)trim(Hx_filename)
     STOP
     
9020 CALL write_line('Error opening Hy data file',0,.TRUE.)
     CALL write_line('filename:',0,.TRUE.)
     write(*,*)trim(Hy_filename)
     STOP 1
      
9030 CALL write_line('Error reading Hy data file',0,.TRUE.)
     CALL write_line('filename:',0,.TRUE.)
     write(*,*)trim(Hy_filename)
     STOP


END SUBROUTINE equivalent_source_model
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
! SUBROUTINE 
!
! NAME electric_dipole_source
!    
!
! DESCRIPTION
!    Calculate the contributions to the radiated field components Hx and Hy at x,y,z due to an 
!    electric dipole with components Px,Py,Pz at x0,y0,z0
!    Contributions from an image source in a PEC plane normal to z can be added if requred
!    
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/3/2024 CJS
!
SUBROUTINE electric_dipole_source(xs0,ys0,zs0,x,y,z,k0,Hx_Px,Hx_Py,Hx_Pz,Hy_Px,Hy_Py,Hy_Pz,add_image,z_PEC)

USE constants

IMPLICIT NONE

real*8 ::xs0,ys0,zs0,x,y,z,k0

complex*16 :: Hx_Px,Hx_Py,Hx_Pz,Hy_Px,Hy_Py,Hy_Pz

logical :: add_image

real*8  :: z_PEC

! local variables

complex*16 :: Hz_Px,Hz_Py,Hz_Pz

real*8     :: r
complex*16 :: C1

real*8 ::xs0p,ys0p,zs0p          ! Image coordinates

! START

  r=sqrt( (x-xs0)**2+(y-ys0)**2+(z-zs0)**2 )
  
  C1=(j*k0*exp(-j*k0*r)/(4d0*pi*r*r) )*(1d0+1d0/(j*k0*r))

! x->y y->z z->x
  Hx_Px=(0d0,0d0)  
  Hy_Px=-C1*(z-zs0)
  Hz_Px= C1*(y-ys0)

! x->z y->x z->y
  Hx_Py= C1*(z-zs0)
  Hy_Py=(0d0,0d0)  
  Hz_Py=-C1*(x-xs0)

! x->x y->y z->z          ! ORIGINAL, not rotated
  Hx_Pz=-C1*(y-ys0)
  Hy_Pz= C1*(x-xs0)
  Hz_Pz=(0d0,0d0)
  
  if (add_image) then

! Position of image source 
    xs0p=xs0
    ys0p=ys0
    zs0p=2d0*z_PEC-zs0 
  
    r=sqrt( (x-xs0p)**2+(y-ys0p)**2+(z-zs0p)**2 )
  
    C1=(j*k0*exp(-j*k0*r)/(4d0*pi*r*r) )*(1d0+1d0/(j*k0*r))

! Reverse Px, Py

! x->y y->z z->x
    Hx_Px=(0d0,0d0)  
    Hy_Px= C1*(z-zs0p)
    Hz_Px=-C1*(y-ys0p)

! x->z y->x z->y
    Hx_Py=-C1*(z-zs0p)
    Hy_Py=(0d0,0d0)  
    Hz_Py= C1*(x-xs0p)

! x->x y->y z->z          ! ORIGINAL, not rotated
    Hx_Pz=-C1*(y-ys0p)
    Hy_Pz= C1*(x-xs0p)
    Hz_Pz=(0d0,0d0)
  
  end if

  RETURN
  
end SUBROUTINE electric_dipole_source

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
! SUBROUTINE 
!
! NAME magnetic_dipole_source
!    
!
! DESCRIPTION
!    Calculate the contributions to the radiated field components Hx and Hy at x,y,z due to an 
!    magnetic dipole with components Mx,My,Mz at x0,y0,z0
!    Contributions from an image source in a PEC plane normal to z can be added if requred
!    
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/3/2024 CJS
!
SUBROUTINE magnetic_dipole_source(xs0,ys0,zs0,x,y,z,k0,Hx_Mx,Hx_My,Hx_Mz,Hy_Mx,Hy_My,Hy_Mz,add_image,z_PEC)

USE constants

IMPLICIT NONE

real*8 ::xs0,ys0,zs0,x,y,z,k0

complex*16 :: Hx_Mx,Hx_My,Hx_Mz,Hy_Mx,Hy_My,Hy_Mz

logical :: add_image

real*8  :: z_PEC

! local variables

complex*16 :: Hz_Mx,Hz_My,Hz_Mz

real*8     :: r
complex*16 :: C1,C2,C3a,C3b

real*8 ::xs0p,ys0p,zs0p          ! Image coordinates

! START

  r=sqrt( (x-xs0)**2+(y-ys0)**2+(z-zs0)**2 )
  
  C1=(j*k0*exp(-j*k0*r)/(4d0*pi*(r**4)) )*(j*k0*r+3d0+3d0/(j*k0*r))
  
  C2=C1
  
  C3a=(j*k0*k0*exp(-j*k0*r)/(4d0*pi*(r**3)) )*(j+3d0/(k0*r)+3d0/(j*k0*k0*r*r))
  C3b=(j*k0*k0*exp(-j*k0*r)/(4d0*pi*r) )*(j+1d0/(k0*r)+1d0/(j*k0*k0*r*r))

! x->y y->z z->x
  Hy_Mx=C1*(y-ys0)*(z-zs0)
  Hz_Mx=C2*(z-zs0)*(x-xs0)
  Hx_Mx=C3a*(x-xs0)*(x-xs0)-C3b

! x->z y->x z->y
  Hz_My=C1*(z-zs0)*(x-xs0)
  Hx_My=C2*(x-xs0)*(y-ys0)
  Hy_My=C3a*(y-ys0)*(y-ys0)-C3b

! x->x y->y z->z          ! ORIGINAL, not rotated
  Hx_Mz=C1*(x-xs0)*(y-ys0)
  Hy_Mz=C2*(y-ys0)*(z-zs0)
   Hz_Mz=C3a*(z-zs0)*(z-zs0)-C3b
  
  if (add_image) then

! Position of image source 
    xs0p=xs0
    ys0p=ys0
    zs0p=2d0*z_PEC-zs0 
  
    r=sqrt( (x-xs0p)**2+(y-ys0p)**2+(z-zs0p)**2 )
  
    C1=(j*k0*exp(-j*k0*r)/(4d0*pi*(r**4)) )*(j*k0*r+3d0+3d0/(j*k0*r))
  
    C2=C1
  
    C3a=(j*k0*k0*exp(-j*k0*r)/(4d0*pi*(r**3)) )*(j+3d0/(k0*r)+3d0/(j*k0*k0*r*r))
    C3b=(j*k0*k0*exp(-j*k0*r)/(4d0*pi*r) )*(j+1d0/(k0*r)+1d0/(j*k0*k0*r*r))

! Reverse Mz only

! x->y y->z z->x
    Hy_Mx=C1*(y-ys0p)*(z-zs0p)
    Hz_Mx=C2*(z-zs0p)*(x-xs0p)
    Hx_Mx=C3a*(x-xs0p)*(x-xs0p)-C3b

! x->z y->x z->y
    Hz_My=C1*(z-zs0p)*(x-xs0p)
    Hx_My=C2*(x-xs0p)*(y-ys0p)
    Hy_My=C3a*(y-ys0p)*(y-ys0p)-C3b

! x->x y->y z->z          ! ORIGINAL, not rotated
    Hx_Mz=-C1*(x-xs0p)*(y-ys0p)
    Hy_Mz=-C2*(y-ys0p)*(z-zs0p)
    Hz_Mz=-(C3a*(z-zs0p)*(z-zs0p)-C3b)
  
  end if

  RETURN
  
end SUBROUTINE magnetic_dipole_source

