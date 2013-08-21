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
! SUBROUTINE read_problem_name
!
! NAME
!     read problem name
!
! DESCRIPTION
!     read problem name
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/08/2012 CJS
!
!
SUBROUTINE read_problem_name

USE TLM_general

IMPLICIT NONE
 
! local variables

  integer :: proc,nchar,talk_to

! START

  if (rank.eq.0) then

    write(*,*)'Enter the name of the TLM input file (without.inp extension)'
    read(*,'(A256)')problem_name
    
  end if ! rank.eq.0
    
#if defined(MPI)

  if (rank.eq.0) then
    
    if (np.gt.1) then
    
! send problem name to other processors   
      do proc=1,np-1
      
        talk_to=proc
        nchar=256
        CALL MPI_SEND(problem_name,nchar,MPI_CHARACTER,talk_to,0,MPI_COMM_WORLD,status,ierror)

      end do
    
    end if
    
  else
  
! read problem name from processor 0  
    talk_to=0 
    nchar=256
    CALL MPI_RECV(problem_name,nchar,MPI_CHARACTER,talk_to,0,MPI_COMM_WORLD,status,ierror)

  end if

#endif
  
  write(*,*)'rank=',rank,' Problem name:',trim(problem_name)
  
  CALL write_line('Problem name:'//trim(problem_name),0,output_to_screen_flag)

  RETURN

END
