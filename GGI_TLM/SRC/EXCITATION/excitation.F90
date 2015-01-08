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
! SUBROUTINE set_excitations_in_mesh
! SUBROUTINE initialise_excitations
! SUBROUTINE initialise_excitation_functions
! SUBROUTINE excitation
! SUBROUTINE get_interpolated_excitation_value
!!
! NAME
!     set_excitations_in_mesh
!
! DESCRIPTION
!     loop over all the required excitations and flag all the excitation cells/ faces 
!     in the arrays local_cell_excitation(i,j,k) or local_surface_excitation(i,j,k,face)
!     as required
!     
! COMMENTS
!     Excitation_surfaces could be made more efficient in parallel - every process gets the whole excitation plane...
!  
!
! HISTORY
!
!     started 17/09/2012 CJS
!     Parallel 23/11/2012 CJS
!     Huygens surface 4/12/2012 CJS
!
!
SUBROUTINE set_excitations_in_mesh

USE TLM_general
USE geometry
USE cell_parameters
USE mesh
USE TLM_excitation
USE file_information

IMPLICIT NONE

! local variables

! START

  
  CALL write_line('CALLED: set_excitations_in_mesh',0,output_to_screen_flag)
  
  CALL set_excitation_points_in_mesh()
  
  CALL set_excitation_surfaces_in_mesh()
  
  CALL set_huygens_surface_in_mesh()

  CALL set_mode_excitations_in_mesh()
  
  CALL write_line('FINISHED: set_excitations_in_mesh',0,output_to_screen_flag)

  RETURN

END SUBROUTINE set_excitations_in_mesh
!
! NAME
!     initialise_excitations
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!     Huygens surface 4/12/2012 CJS
!
!
SUBROUTINE initialise_excitations

USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE mesh
USE TLM_excitation
USE file_information
USE constants

IMPLICIT NONE

! local variables

  integer	:: excitation_surface
  integer	:: excitation_point
  integer	:: cx,cy,cz,face
  integer	:: number_of_faces
  integer	:: excitation_face
  integer	:: huygens_face

  integer 	:: face_number
  integer 	:: number_of_points
  integer 	:: point_number
  type(xyz)	:: point1,point2,point3,point4
    
! START

  
  CALL write_line('CALLED: initialise_excitations',0,output_to_screen_flag)

  CALL initialise_excitation_points()

  CALL initialise_excitation_surfaces()

  CALL initialise_huygens_surface()

  CALL initialise_mode_excitations()
  
  CALL write_line('FINISHED: initialise_excitations',0,output_to_screen_flag)

  RETURN

END SUBROUTINE initialise_excitations
!
! NAME
!     initialise_excitation_functions
!
! DESCRIPTION
!     
!     Set up the sampled excitation functions in arrays 1:n_timesteps
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!
!
SUBROUTINE initialise_excitation_functions

USE TLM_general
USE TLM_excitation
USE mesh
USE file_information
USE output_formats
USE constants

IMPLICIT NONE

! local variables

  integer	:: excitation_number
  
  real*8	:: amplitude
  real*8	:: delay,width
  real*8	:: frequency,phase
  real*8	:: gaussian_function
  real*8	:: alpha
  real*8	:: beta
  real*8	:: tpeak
  real*8	:: fpeak
  real*8	:: r_random
  real*8	:: t1,t2,t3,t4,pulse_amplitude
  
  real*8 	:: interpolate_function

! START

  if (n_excitation_functions.LE.0) RETURN
  
  CALL write_line('CALLED: initialise_excitation_functions',0,output_to_screen_flag)
  
  OPEN(unit=excitation_output_unit,file=trim(problem_name)//excitation_output_extn)
  
  CALL write_time_domain_header_data(excitation_output_unit,n_excitation_functions,n_timesteps)

  do excitation_number=1,n_excitation_functions
  
    ALLOCATE( excitation_functions(excitation_number)%value(1:n_timesteps) )
    ALLOCATE( excitation_functions(excitation_number)%value_face(1:n_timesteps) )
    
    do timestep=1,n_timesteps
    
      time=(timestep-1)*dt

      if (excitation_functions(excitation_number)%type.EQ.excitation_function_type_impulse) then

        amplitude=excitation_functions(excitation_number)%parameters(1)
        if (timestep.eq.1) then
	  excitation_functions(excitation_number)%value(timestep)=amplitude
        else
	  excitation_functions(excitation_number)%value(timestep)=0d0
	end if
        if (timestep.eq.1) then
	  excitation_functions(excitation_number)%value_face(timestep)=amplitude
        else
	  excitation_functions(excitation_number)%value_face(timestep)=0d0
	end if
	
      else if (excitation_functions(excitation_number)%type.EQ.excitation_function_type_gaussian) then
      
        amplitude=excitation_functions(excitation_number)%parameters(1)
        width=excitation_functions(excitation_number)%parameters(2)
        delay=excitation_functions(excitation_number)%parameters(3)

! evaluate gaussian function, setting result to zero if we are too far into the tails of the function 	
	if (abs(time-delay).lt.width*5d0) then
	  excitation_functions(excitation_number)%value(timestep)=	&
		                               amplitude*( exp(-((time-delay)/width)**2) )
	else
	  excitation_functions(excitation_number)%value(timestep)=0D0
	end if
	
	if (abs(time+dt/2d0-delay).lt.width*5) then
	  excitation_functions(excitation_number)%value_face(timestep)=	&
		                               amplitude*( exp(-((time+dt/2d0-delay)/width)**2) )
	else
	  excitation_functions(excitation_number)%value_face(timestep)=0D0
	end if
	
      else if (excitation_functions(excitation_number)%type.EQ.excitation_function_type_gaussian_step) then
      
        amplitude=excitation_functions(excitation_number)%parameters(1)
        width=excitation_functions(excitation_number)%parameters(2)
        delay=excitation_functions(excitation_number)%parameters(3)

! evaluate gaussian function, setting result to zero if we are too far into the tails of the function 	
	if (time.LT.delay) then
	
	  if (abs(time-delay).lt.width*5) then
	    excitation_functions(excitation_number)%value(timestep)=	&
		                               amplitude*( exp(-((time-delay)/width)**2) )
	  else
	    excitation_functions(excitation_number)%value(timestep)=0D0
	  end if
	  
	else
	  excitation_functions(excitation_number)%value(timestep)=amplitude
	end if
	
	if (time-dt/2d0.LT.delay) then
	
	  if (abs(time-dt/2d0-delay).lt.width*5) then
	    excitation_functions(excitation_number)%value_face(timestep)=	&
		                               amplitude*( exp(-((time-dt/2d0-delay)/width)**2) )
	  else
	    excitation_functions(excitation_number)%value_face(timestep)=0D0
	  end if
	  
	else
	  excitation_functions(excitation_number)%value_face(timestep)=amplitude
	end if
	
      else if (excitation_functions(excitation_number)%type.EQ.excitation_function_type_step) then
      
        amplitude=excitation_functions(excitation_number)%parameters(1)
        delay=excitation_functions(excitation_number)%parameters(2)

! evaluate gaussian function, setting result to zero if we are too far into the tails of the function 	
	if (time.GE.delay) then
	  excitation_functions(excitation_number)%value(timestep)=amplitude
	else
	  excitation_functions(excitation_number)%value(timestep)=0D0
	end if
	
	if (time+dt/2d0.GE.delay) then
	  excitation_functions(excitation_number)%value_face(timestep)=amplitude
	else
	  excitation_functions(excitation_number)%value_face(timestep)=0D0
	end if
	
      else if (excitation_functions(excitation_number)%type.EQ.excitation_function_type_sinusoid) then
      
        amplitude=excitation_functions(excitation_number)%parameters(1)
        frequency=excitation_functions(excitation_number)%parameters(2)
        phase=excitation_functions(excitation_number)%parameters(3)*pi/180d0

	excitation_functions(excitation_number)%value(timestep)=	&
                             amplitude*sin(2d0*pi*frequency*time-phase)

	excitation_functions(excitation_number)%value_face(timestep)=	&
                             amplitude*sin(2d0*pi*frequency*(time+dt/2d0)-phase)
	
      else if (excitation_functions(excitation_number)%type.EQ.excitation_function_type_gaussian_sinusoid) then
      
        amplitude=excitation_functions(excitation_number)%parameters(1)
        width=excitation_functions(excitation_number)%parameters(2)
        delay=excitation_functions(excitation_number)%parameters(3)
        frequency=excitation_functions(excitation_number)%parameters(4)
        phase=excitation_functions(excitation_number)%parameters(5)*pi/180d0
	
! evaluate gaussian function, setting result to zero if we are too far into the tails of the function 	
	if (abs(time-delay).lt.width*5d0) then
	  gaussian_function=( exp(-((time-delay)/width)**2) )
	else
	  gaussian_function=0D0
	end if

	excitation_functions(excitation_number)%value(timestep)=	&
                             amplitude*gaussian_function*sin(2d0*pi*frequency*(time-delay)-phase)

! half timestep evaluation	
	if (abs(time+dt/2d0-delay).lt.width*5) then
	  gaussian_function=( exp(-((time+dt/2d0-delay)/width)**2) )
	else
	  gaussian_function=0D0
	end if

	excitation_functions(excitation_number)%value_face(timestep)=	&
                             amplitude*gaussian_function*sin(2d0*pi*frequency*(time+dt/2d0)-phase)

	
      else if (excitation_functions(excitation_number)%type.EQ.excitation_function_type_gaussian_step_sinusoid) then
      
        amplitude=excitation_functions(excitation_number)%parameters(1)
        width=excitation_functions(excitation_number)%parameters(2)
        delay=excitation_functions(excitation_number)%parameters(3)
        frequency=excitation_functions(excitation_number)%parameters(4)
        phase=excitation_functions(excitation_number)%parameters(5)*pi/180d0
	
! evaluate gaussian function, setting result to zero if we are too far into the tails of the function 	
	if (time.LT.delay) then
	
	  if (abs(time-delay).lt.width*5) then
	    gaussian_function=( exp(-((time-delay)/width)**2) )
	  else
	    gaussian_function=0D0
	  end if
	  
	else
	  gaussian_function=1D0
	end if

	excitation_functions(excitation_number)%value(timestep)=	&
                             amplitude*gaussian_function*sin(2d0*pi*frequency*time-phase)

! half timestep evaluation	
	
	if (time-dt/2d0.LT.delay) then
	
	  if (abs(time-dt/2d0-delay).lt.width*5) then
	    gaussian_function=( exp(-((time-dt/2d0-delay)/width)**2) )
	  else
	    gaussian_function=0D0
	  end if
	  
	else
	  gaussian_function=1D0
	end if

	excitation_functions(excitation_number)%value_face(timestep)=	&
                             amplitude*gaussian_function*sin(2d0*pi*frequency*(time+dt/2d0)-phase)

      else if (excitation_functions(excitation_number)%type.EQ.excitation_function_type_double_exponential) then

        if (time.ge.0d0) then
    
          amplitude=1d0/excitation_functions(excitation_number)%parameters(1)
          alpha=1d0/excitation_functions(excitation_number)%parameters(2)
          beta=1d0/excitation_functions(excitation_number)%parameters(3)
          tpeak=(log(beta)-log(alpha))/(beta-alpha)
          fpeak=(exp(-alpha*tpeak)-exp(-beta*tpeak))
      
          excitation_functions(excitation_number)%value(timestep)=	&
              (amplitude/fpeak)*(exp(-alpha*time)-exp(-beta*time))
      
        else
    
          excitation_functions(excitation_number)%value(timestep)=0d0
      
        end if  

! half timestep evaluation	

        if (time-dt/2d0.ge.0d0) then
    
          amplitude=1d0/excitation_functions(excitation_number)%parameters(1)
          alpha=1d0/excitation_functions(excitation_number)%parameters(2)
          beta=1d0/excitation_functions(excitation_number)%parameters(3)
          tpeak=(log(beta)-log(alpha))/(beta-alpha)
          fpeak=(exp(-alpha*tpeak)-exp(-beta*tpeak))
                                 
          excitation_functions(excitation_number)%value_face(timestep)=	&
              (amplitude/fpeak)*(exp(-alpha*time-dt/2d0)-exp(-beta*time-dt/2d0))
      
        else
    
          excitation_functions(excitation_number)%value_face(timestep)=0d0
      
        end if  
	
      else if (excitation_functions(excitation_number)%type.EQ.excitation_function_type_differential_gaussian) then
      
        amplitude=excitation_functions(excitation_number)%parameters(1)
        width    =excitation_functions(excitation_number)%parameters(2)
        delay    =excitation_functions(excitation_number)%parameters(3)
	
        tpeak=-width/sqrt(2d0)
        fpeak=-2d0*(tpeak/(width**2))*( exp(-(tpeak/width)**2) )

! evaluate differential gaussian function, setting result to zero if we are too far into the tails of the function 	
	if (abs(time-delay).lt.width*5d0) then
	  excitation_functions(excitation_number)%value(timestep)=	&
		              -2d0*((time-delay)/width**2)*amplitude*( exp(-((time-delay)/width)**2) )/fpeak
	else
	  excitation_functions(excitation_number)%value(timestep)=0D0
	end if
	
	if (abs(time+dt/2d0-delay).lt.width*5) then
	  excitation_functions(excitation_number)%value_face(timestep)=	&
		              -2d0*((time+dt/2d0-delay)/width**2)*amplitude*( exp(-((time+dt/2d0-delay)/width)**2) )/fpeak
	else
	  excitation_functions(excitation_number)%value_face(timestep)=0D0
	end if
	
      else if (excitation_functions(excitation_number)%type.EQ.excitation_function_type_noise) then
      
        amplitude=excitation_functions(excitation_number)%parameters(1)
	
        CALL random_number(r_random)
	excitation_functions(excitation_number)%value(timestep)=amplitude*(r_random-0.5d0)
	  
        CALL random_number(r_random)
	excitation_functions(excitation_number)%value_face(timestep)=	amplitude*(r_random-0.5d0)
	
      else if (excitation_functions(excitation_number)%type.EQ.excitation_function_type_sinusoidal_pulse) then
      
        amplitude=excitation_functions(excitation_number)%parameters(1)
        frequency=excitation_functions(excitation_number)%parameters(2)
        phase=excitation_functions(excitation_number)%parameters(3)*pi/180d0
        t1=excitation_functions(excitation_number)%parameters(4)
        t2=excitation_functions(excitation_number)%parameters(5)
        t3=excitation_functions(excitation_number)%parameters(6)
        t4=excitation_functions(excitation_number)%parameters(7)

! work out the trapezoidal pulse amplitude (between 0 and 1)	
	pulse_amplitude=0d0
	if ( (time.ge.t1).AND.(time.lt.t2) ) then
	  pulse_amplitude=(time-t1)/(t2-t1)
	else if ( (time.ge.t2).AND.(time.le.t3) ) then
	  pulse_amplitude=1d0
	else if ( (time.ge.t3).AND.(time.le.t4) ) then
	  pulse_amplitude=(t4-time)/(t4-t3)
	end if
	  
	excitation_functions(excitation_number)%value(timestep)=	&
	     amplitude*pulse_amplitude*sin(2d0*pi*frequency*time-phase)
	  
! half timestep excitation	  
! work out the trapezoidal pulse amplitude (between 0 and 1)	
	pulse_amplitude=0d0
	if ( (time+dt/2d0.ge.t1).AND.(time+dt/2d0.lt.t2) ) then
	  pulse_amplitude=(time+dt/2d0-t1)/(t2-t1)
	else if ( (time+dt/2d0.ge.t2).AND.(time+dt/2d0.le.t3) ) then
	  pulse_amplitude=1d0
	else if ( (time+dt/2d0.ge.t3).AND.(time+dt/2d0.le.t4) ) then
	  pulse_amplitude=(t4-(time+dt/2d0))/(t4-t3)
	end if
	  
	excitation_functions(excitation_number)%value_face(timestep)=	&
	     amplitude*pulse_amplitude*sin(2d0*pi*frequency*(time+dt/2d0)-phase)

      else if (excitation_functions(excitation_number)%type.EQ.excitation_function_type_file) then
      
        amplitude=excitation_functions(excitation_number)%parameters(1)
        delay=excitation_functions(excitation_number)%parameters(2)

! evaluate excitation function	

	excitation_functions(excitation_number)%value(timestep)=					&
	amplitude*interpolate_function(time-delay,							&
	                               excitation_functions(excitation_number)%n_values_from_file,	&
	                               excitation_functions(excitation_number)%time_values_from_file,	&
	                               excitation_functions(excitation_number)%function_values_from_file)
	
	excitation_functions(excitation_number)%value_face(timestep)=					&
	amplitude*interpolate_function(time+dt/2d0-delay,						&
	                               excitation_functions(excitation_number)%n_values_from_file,	&
	                               excitation_functions(excitation_number)%time_values_from_file,	&
	                               excitation_functions(excitation_number)%function_values_from_file)
 	
      end if  ! excitation_function_type
      
      write(excitation_output_unit,time_domain_output_format)time,excitation_number,	&
                                   excitation_functions(excitation_number)%value(timestep)

    end do ! next timestep

  end do ! next excitation function

  CLOSE(unit=excitation_output_unit)
  
  CALL write_line('FINISHED: initialise_excitation_functions',0,output_to_screen_flag)

  RETURN

END SUBROUTINE initialise_excitation_functions
!
! SUBROUTINE excitation
!
! NAME
!     excitation
!
! DESCRIPTION
!     
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!     Huygens surface 4/12/2012 CJS
!
!
SUBROUTINE excitation

USE TLM_general
USE mesh
USE TLM_excitation
USE constants

IMPLICIT NONE

! local variables

  integer 	:: excitation_point
  integer 	:: excitation_surface
  integer 	:: huygens_face
  integer 	:: face
  integer 	:: cz
  integer 	:: number_of_faces
  integer 	:: function_number
  integer 	:: excitation_array_point
  integer 	:: field_component
  real*8	:: offset,offset_min
  real*8 	:: value  
  
  real*8	:: Js(3),Ms(3)
  real*8	:: normx,normy,normz

! START
  
  CALL write_line('CALLED: excitation',0,timestepping_output_to_screen_flag)

! reset all excitations   
  if (allocated( cell_excitation_field )) cell_excitation_field(:,:)=0d0
  if (allocated( face_excitation_field )) face_excitation_field(:,:,:)=0d0
 
  CALL Point_excitation()
  
  CALL Surface_excitation()
  
  CALL Huygens_surface_excitation()
  
  CALL Mode_excitation()
  
  CALL write_line('FINISHED: excitation',0,timestepping_output_to_screen_flag)
  
  RETURN

END SUBROUTINE excitation
!
! SUBROUTINE get_interpolated_excitation_value
!
! NAME
!     get_interpolated_excitation_value
!
! DESCRIPTION
!     Interpolate the excitation value from the 
!     
! COMMENTS
!     
!  
!
! HISTORY
!
!     started  4/12/2012 CJS
!
!
SUBROUTINE get_interpolated_excitation_value(offset,offset_min,function_number,value)


USE TLM_general
USE TLM_excitation
USE constants

IMPLICIT NONE

  real*8	:: offset,offset_min
  integer	:: function_number
  real*8	:: value

! local variables

  integer 	:: timestep1,timestep2
  real*8	:: delta_t,t_offset,local_time
  real*8 	:: value1,value2
    
! START

  t_offset=(offset-offset_min)/c0
  local_time=(dble(timestep-1))*dt-t_offset  ! note this is the time at the start of the timestep, 
                                             ! not at the half timestep connect
  
  timestep1=NINT(local_time/dt)
  timestep2=timestep1+1
  delta_t=local_time-dble(timestep1)*dt
  
  if ((timestep1.ge.1).AND.(timestep1.le.n_timesteps) ) then
    value1=excitation_functions(function_number)%value_face(timestep1)
  else
    value1=0d0
  end if
  
  if ((timestep2.ge.1).AND.(timestep2.le.n_timesteps) ) then
    value2=excitation_functions(function_number)%value_face(timestep2)
  else
    value2=0d0
  end if
  
  value=value1+delta_t*(value2-value1)/dt
  
  RETURN

END SUBROUTINE get_interpolated_excitation_value
!
! FUNCTION interpolate_function
!
! NAME
!     interpolate_function
!
! DESCRIPTION
!     interpolate/ extyrapolate excitation values defined in a file
!
!
! COMMENTS
!     
!
! HISTORY
!
!     started 23/9/2013
!
!

FUNCTION interpolate_function(time,n_values,time_values,function_values)

real*8 	:: interpolate_function

real*8 	:: time
integer	:: n_values
real*8 	:: time_values(1:n_values)
real*8 	:: function_values(1:n_values)

! local variables

integer 	:: i
real*8		:: tmin,tmax
real*8        	:: delta_t,time1,time2
real*8        	:: value1,value2

! START

  tmin=time_values(1)
  tmax=time_values(n_values)
  
  if (time.LE.tmin) then
    interpolate_function=function_values(1)
    RETURN
  else if (time.GE.tmax) then
    interpolate_function=function_values(n_values)
    RETURN
  else
  
    do i=1,n_values-1
    
      time1=time_values(i)
      time2=time_values(i+1)
      
      if ( (time.GE.time1).AND.(time.LE.time2) ) then
!  calculate interpolated value and return
  
        value1=function_values(i)
        value2=function_values(i+1)
  
        interpolate_function=value1+(time-time1)*(value2-value1)/(time2-time1)
	
	RETURN
	
      end if

    end do  ! next value in file data

  end if
  
  write(*,*)'Error in interpolate_function'
  write(*,*)'No value found'
  write(*,*)'time=',time
  write(*,*)'tmin=',tmin
  write(*,*)'tmax=',tmax
  STOP


END FUNCTION interpolate_function
