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
! SUBROUTINE diode_calc
!
! NAME
!     diode_calc
!
! DESCRIPTION
!     
! diode update
!     
!     
! COMMENTS
!     Uses a Newton Raphson method to calculate the voltage and current on
!     a diode with series voltage source and resistance. 
!     Called from Connect for surface impedance diode model
!     Iterative solution for the diode voltage - two iterative formulae are used to try and ensure
!     convergence in all circumstances
!
! HISTORY
!
!     started 3/09/2014 CJS
!     9/9/2014  Include junction capacitance model, also improve the iterative solution so that the 
!               convergence is more robust
!
!
SUBROUTINE diode_calc(Vs_in,R_in,Is,nVt,Rs,Zc,Vcs,Vd,Id,Ic,direction)


IMPLICIT NONE

real*8 	:: Vs_in,R_in,Is,nVt,Rs,Zc,Vcs,Vd,Id,Ic
integer :: direction

! local parameters

real*8 	:: Vs,R
real*8 	:: cvg_limit=1d-8
integer,parameter	:: itloop_max=100

! local variables

real*8 	:: Vk
real*8 	:: last_Vd
real*8 	:: cvg
real*8 	:: Zd

real*8	:: Vc,VRs

integer :: itloop

! START
  
! build Thevenin equivalent of the linear TLM components  
  
  Vs=((Vs_in*direction)/(R_in+Rs)+Vcs/Zc)/	&
      (1d0/(R_in+Rs)+1d0/(Zc))  ! reverse diode if requested
  
  R=(R_in+Rs)*Zc/(R_in+Rs+Zc)
  
! Assume I=0 if |V|<some small value- prevents problems with the numerical convergence in this case

  if (abs(Vs).LT.cvg_limit) then
  
    Vd=0d0
    
  else
   
! solve for the diode voltage turning point, vk, assumed to be the voltage at which the 
! current is 10**8 times Is

    Vk=nVt*log(1d8)

! intial guess at solution

    if (Vs.LE.Vk) then ! assume reverse bias diode, Vd=Vs  
    
      Vd=Vs
  
    else ! assume a conducting diode, Vd =Vk
  
      Vd=Vk
  
    end if

! iterative solution

    do itloop=1,itloop_max

      last_Vd=Vd
	    
      if ( ((Vs-Vd)/(R*Is))+1d0.LE.0d0 ) then
! This  problems for the fast iterative solution, try a different form of the iteration

        Vd=Vd-(Vs-R*Is*(exp(Vd/nVt)-1d0)-Vd)/	&
              (-R*Is*exp(Vd/nVt)/nVt-1d0)

      else

! FAST CONVERGENT SOLUTION
        Vd=Vd-(nVt*log((Vs-Vd)/(R*Is)+1d0)-Vd)/	    &
    	      (-nVt/(R*Is*((Vs-Vd)/(R*Is)+1d0))-1d0)
	    
      end if

      cvg=abs(Vd-last_Vd)
      if (abs(Vd).GT.1d0) then
        cvg=cvg/abs(Vd)
      end if

      if (cvg.LT.cvg_limit) GOTO 1000

    end do
  
    write(*,*)'Iterative diode solution not converged...'
    write(*,*)'Vd	    =',Vd
    write(*,*)'last_Vd      =',last_Vd
    write(*,*)'cvg          =',cvg
    write(*,*)'cvg_limit    =',cvg_limit
    write(*,*)'Vs	    =',Vs
    write(*,*)'R	    =',R
    write(*,*)'Vk	    =',Vk
    write(*,*)'Is	    =',Is
    write(*,*)'nVt	    =',nVt
    
    write(*,*)'Re-run process for de-bugging purposes...'    

! intial guess at solution
    if (Vs.LE.Vk) then ! assume reverse bias diode, Vd=Vs  
      Vd=Vs
    else ! assume a conducting diode, Vd =Vk  
      Vd=Vk  
    end if
    
    write(*,*)0,Vd,Vs
    
    write(*,*)Vd-(nVt*log((Vs-Vd)/(R*Is)+1d0)-Vd)/	    &
    	    (-nVt/(R*Is*((Vs-Vd)/(R*Is)+1d0))-1d0)
    write(*,*)(Vs-Vd)
    write(*,*)(Vs-Vd)/(R*Is)+1d0
    write(*,*)nVt*log((Vs-Vd)/(R*Is)+1d0)
    write(*,*)(nVt*log((Vs-Vd)/(R*Is)+1d0)-Vd)
    write(*,*) (-nVt/(R*Is*((Vs-Vd)/(R*Is)+1d0))-1d0)

! iterative solution
    do itloop=1,itloop_max

      last_Vd=Vd
	    
      if ( ((Vs-Vd)/(R*Is))+1d0.LE.0d0 ) then
! This  problems for the fast iterative solution, try a different form of the iteration

        Vd=Vd-(Vs-R*Is*(exp(Vd/nVt)-1d0)-Vd)/	&
              (-R*Is*exp(Vd/nVt)/nVt-1d0)

      else

! FAST CONVERGENT SOLUTION
        Vd=Vd-(nVt*log((Vs-Vd)/(R*Is)+1d0)-Vd)/	    &
    	      (-nVt/(R*Is*((Vs-Vd)/(R*Is)+1d0))-1d0)
	    
      end if

      cvg=abs(Vd-last_Vd)
      if (abs(Vd).GT.1d0) then
        cvg=cvg/abs(Vd)
      end if
      
      write(*,*)itloop,Vd,cvg

    end do
    
    STOP
  
    1000 CONTINUE ! jump here when converged
  
  end if

! Diode current  
  Id=(Vs-Vd)/R

! check that the impedance of the diode is positive i.e. a lossy device.  
  if (id.NE.0d0) then
  
   Zd=Vd/Id
  
    if (Zd.lt.0D0) then
      write(*,*)'Diode impedance lt 0',Zd
      write(*,*)'Vs	  =',Vs
      write(*,*)'R	  =',R
      write(*,*)'Vk	  =',Vk
      write(*,*)'last_Vd  =',last_Vd
      write(*,*)'Vd	  =',Vd
      write(*,*)'Id       =',Id
      STOP
    end if
    
  end if
  
! From the ideal diode voltage and current, return the voltage and current for the complete diode model

  Vc=Vd
  Ic=(Vc-Vcs)/Zc
  VRs=(Ic+Id)*Rs
  
  Vd=Vd+VRs  ! include voltage across series resistance in the diode voltage 
  Id=Id+Ic   ! include the junction capacitor current in the diode current
    
  Id=Id*direction  ! reverse diode if requested
  Vd=Vd*direction  ! reverse diode if requested

  RETURN

END SUBROUTINE diode_calc
