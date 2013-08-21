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
       PROGRAM RCS_sphere

USE file_information
USE filter_types
USE filter_functions
USE constants

       implicit none
       
! calculate the RCS of a sphere from the Mie series solution

! file names
       character*256 ipfile,opfile

! inner PEC sphere radius
       real*8 a
       
! outer coating radius
       real*8 a2
       
! material and thin layer file names
       character*256 mat_file,thin_layer_file
       
! material and thin layer filter data
  type(Sfilter)	:: epsr_filter
  real*8	:: sigmae
  
  type(Sfilter)	:: mur_filter
  real*8	:: sigmam

  type(Sfilter)	:: z11_filter
  type(Sfilter)	:: z12_filter
  type(Sfilter)	:: z21_filter
  type(Sfilter)	:: z22_filter
       
! coating relative permittivity and permeability
       complex*16 epsr2,mur2
       
! coating permittivity and permeability
       complex*16 eps2,mu2
       
! impedance parameters
       integer thin_layer 
       complex*16 z11,z12,z21,z22       

! frequency range
       real*8 fmin,fmax,df,f
       integer nf,ifreq

! angle range, theta 
       real*8 thetamin,thetamax,dtheta,theta
       integer ntheta,itheta

! angle range phi       
       real*8 phimin,phimax,dphi,phi
       integer nphi,iphi

! constants       
       real*8 w,beta
       
       complex*16 Etheta,Ephi,Htheta,Hphi
       
       real*8 sigma3d
       real*8 sigma3d_db
       
       integer n
       real*8 re,im
              
! START	      
       
       open(unit=99,file='progress')
       write(99,'(A)')'STARTED: GGI_RCS_sphere'
       close(unit=99)
	      
       print*,'Problem Filename'	      
       read(*,'(A256)')ipfile
       open(unit=8,file=ipfile,status='OLD')	      
	      
! set radius of PEC sphere              
       read(8,*)
       read(8,*)a
       print*,'Inner PEC sphere radius=',a
       
! set radius of coating
       read(8,*)
       read(8,*)a2
       print*,'Coating radius=',a2

! read material filename       
       read(8,*)
       read(8,'(A256)')mat_file
       
        
! read material filter data 
       open(UNIT=volume_material_file_unit,						&
           FILE=trim(mat_file),	&
           STATUS='old')
      
       read(volume_material_file_unit,*)    ! read material label      
       read(volume_material_file_unit,*)    ! read frequency range of validity
       read(volume_material_file_unit,*)    ! read comment line
       
! read permittivity filter
       call read_Sfilter(epsr_filter,volume_material_file_unit) 

! read electric conductivity
       read(volume_material_file_unit,*)sigmae

       read(volume_material_file_unit,*)    ! read comment line
      
! read permeability filter
       call read_Sfilter(mur_filter,volume_material_file_unit) 

! read magnetic conductivity
       read(volume_material_file_unit,*)sigmam
      
       close(UNIT=volume_material_file_unit)

! read thin_layer constant. set to 0 for no sheet, 1 to include IBC sheet, 2 for SIBC on PEC                   
       read(8,*)
       read(8,*)thin_layer
! read thin layer filename       
       read(8,*)
       read(8,'(A256)')thin_layer_file 
       
! read thin layer filters  
     
       if (thin_layer.eq.1) then
       
         open(UNIT=surface_material_file_unit,						&
              FILE=trim(thin_layer_file),	&
              STATUS='old')
      
         read(surface_material_file_unit,*)    ! read material label
         read(surface_material_file_unit,*)    ! read frequency range of validity

         read(surface_material_file_unit,*)    ! read comment line	
         call read_Sfilter(z11_filter,surface_material_file_unit) ! read Z11 filter

         read(surface_material_file_unit,*)    ! read comment line
         call read_Sfilter(z12_filter,surface_material_file_unit) ! read Z12 filter

         read(surface_material_file_unit,*)    ! read comment line
         call read_Sfilter(z21_filter,surface_material_file_unit) ! read Z21 filter

         read(surface_material_file_unit,*)    ! read comment line
         call read_Sfilter(z22_filter,surface_material_file_unit) ! read Z22 filter
      
         close(UNIT=surface_material_file_unit)
         
       else if (thin_layer.eq.2) then
       
         open(UNIT=surface_material_file_unit,						&
              FILE=trim(thin_layer_file),	&
              STATUS='old')
      
         read(surface_material_file_unit,*)    ! read material label
         read(surface_material_file_unit,*)    ! read frequency range of validity

         read(surface_material_file_unit,*)    ! read comment line	
         call read_Sfilter(z11_filter,surface_material_file_unit) ! read Z11 filter
      
         close(UNIT=surface_material_file_unit)
         
       end if
                     
! set frequency range     
       read(8,*)
       read(8,*)fmin,fmax,nf 
       print*,'fmin=',fmin      
       print*,'fmax=',fmax      
       print*,'nf  =',nf      
      
! set theta range      
       read(8,*)
       read(8,*)thetamin,thetamax,ntheta
       print*,'thetamin=',thetamin
       print*,'thetamax=',thetamax
       print*,'ntheta  =',ntheta
       
! set phi range
       read(8,*)
       read(8,*)phimin,phimax,nphi
       print*,'phimin=',phimin
       print*,'phimax=',phimax
       print*,'nphi  =',nphi
            
! set number of terms in series solution
       read(8,*)
       read(8,*)n
       print*,'number of terms in expansion=',n
       
! set output filename
       read(8,*)
       read(8,'(A256)')opfile
       print*,'output filename=',opfile
                    
! work out frequency step       
       if (fmin.ne.fmax) then
         df=(fmax-fmin)/nf
       else
         df=1d0
       end if     
       
! convert to radians and work out angle steps
       thetamin=thetamin*pi/180d0       
       thetamax=thetamax*pi/180d0       
       if (thetamin.ne.thetamax) then
         dtheta=(thetamax-thetamin)/(ntheta-1)
       else
         dtheta=1d0
       end if     
       
       phimin=phimin*pi/180d0       
       phimax=phimax*pi/180d0 
       if (phimin.ne.phimax) then
         dphi=(phimax-phimin)/(nphi-1)
       else
         dphi=1d0
       end if
                     
!       call test_legendre()
       
!       call test_bessel()
       
! open output files
      
       open(unit=10,file=opfile)
       
       if (thin_layer.eq.1) then
         open(unit=11,file='z11.fout')
         open(unit=12,file='z12.fout')
         open(unit=21,file='z21.fout')
         open(unit=22,file='z22.fout')
       else if (thin_layer.eq.2) then
         open(unit=11,file='z11.fout')
       end if
       
       open(unit=30,file='epsr.fout')
       open(unit=31,file='mur.fout')
       
       if (a2.eq.a) then
! calculate RCS for PEC sphere only
         if (thin_layer.eq.0) then
 	   print*,'calling rcs'
	 else if (thin_layer.eq.1) then
 	   print*,'Cannot have thin layer on PEC'
	   stop
	 else if (thin_layer.eq.2) then
 	   print*,'calling rcs_ibc'
	 end if
       else
         if (a.eq.0d0) then
	   print*,'calling rcs_dielectric_ibc'
 	 else if(thin_layer.eq.0) then
! calculate RCS for coated sphere only
 	   print*,'calling rcs_coated'
 	 else if(thin_layer.eq.2) then
! calculate RCS for coated sphere only
 	   print*,'Cannot have SIBC on dielectric'
	   STOP
 	 else
! calculate RCS for coated sphere with thin layer
 	   print*,'calling rcs_coated_ibc'
 	 end if
       end if

! loop over frequency       

       do ifreq=1,nf
       
         f=fmin+(ifreq-1)*df
	 
! set parameters       
         w=2d0*pi*f       
         beta=w/c0

! calculate material parameters at this frequency

	 epsr2=evaluate_Sfilter_frequency_response(epsr_filter,f)-j*sigmae/(w*eps0)
! write  material parameters to file 
         write(30,1000)f,real(epsr2),aimag(epsr2)
	 
	 mur2=evaluate_Sfilter_frequency_response(mur_filter,f)-j*sigmam/(w*mu0)
! write  material parameters to file 
         write(31,1000)f,real(mur2),aimag(mur2)
	 
! constants for coating
         eps2=eps0*epsr2
         mu2=mu0*mur2 
             
! calculate impedance boundary parameters at this frequency
         if (thin_layer.eq.1) then
       
           z11=evaluate_Sfilter_frequency_response(z11_filter,f)
           z12=evaluate_Sfilter_frequency_response(z12_filter,f)
           z21=evaluate_Sfilter_frequency_response(z21_filter,f)
           z22=evaluate_Sfilter_frequency_response(z22_filter,f)
              
! change the signs of z12 and z22, this is required due to the definition of
! z parameters used i.e. normal convention for 2 port networks
           z12=-z12
           z22=-z22   
	   
! write impedance boundary parameters to file 
           write(11,1000)f,real(z11),aimag(z11)
           write(12,1000)f,real(z12),aimag(z12)
           write(21,1000)f,real(z21),aimag(z21)
           write(22,1000)f,real(z22),aimag(z22)
	 
1000       format(3E14.6)

         else if (thin_layer.eq.2) then
       
           z11=evaluate_Sfilter_frequency_response(z11_filter,f)

! write impedance boundary parameters to file 
           write(11,1000)f,real(z11),aimag(z11)

         end if  
	                                  
! loop over angles
	
	 do itheta=1,ntheta
	 
	   theta=thetamin+(itheta-1)*dtheta
	   
	   do iphi=1,nphi
	   
	   phi=phimin+(iphi-1)*dphi
	   
!	     write(*,1001)f,fmax,theta*180d0/pi,thetamax*180d0/pi,phi*180d0/pi,phimax*180d0/pi
1001         format (2E14.7,4F10.4)  

             Etheta=(0d0,0d0)
             Ephi  =(0d0,0d0)
             Htheta=(0d0,0d0)
             Hphi  =(0d0,0d0)

             if (a2.eq.a) then
! calculate RCS for PEC sphere only
               if (thin_layer.eq.0) then
                 call rcs(theta,phi,beta,a,n,sigma3d)
	       else
                 call rcs_ibc(theta,phi,w,a,n,z11,sigma3d)
	       end if
	     else if (a.eq.0d0) then
	       if(thin_layer.eq.0) then
       print*,'No dielectric sphere code without resistive sheet as yet'
		 stop 
	       else
                 call rcs_dielectric_ibc(theta,phi,w,a2,eps2,mu2,	&
                 z11,z12,z21,z22,n,Etheta,Ephi,Htheta,Hphi,sigma3d)
	       end if	     
	     else
	       if(thin_layer.eq.0) then
! calculate RCS for coated sphere only
                 call rcs_coated(theta,phi,w,a,a2,eps2,mu2,n		&
                                ,Etheta,Ephi,Htheta,Hphi,sigma3d)
	       else
! calculate RCS for coated sphere with thin layer
                 call rcs_coated_ibc(theta,phi,w,a,a2,eps2,mu2,		&
                  z11,z12,z21,z22,n,Etheta,Ephi,Htheta,Hphi,sigma3d)
	       end if
	     end if

! convert to dB	 
             if (sigma3d.gt.0d0) then
	       sigma3d_db=10d0*log10(sigma3d)
	     else
	       sigma3d_db=-200.
	     end if	 

! write to file	 
             if (a.ne.0d0) then
	       write(10,1020)f,180d0-theta*180./pi,180d0-phi*180./pi,	&
     			  sigma3d_db,				       &
     		      real(Etheta),imag(Etheta),real(Ephi),imag(Ephi)
     	     else
               write(10,1020)f,180d0-theta*180./pi,180d0-phi*180./pi,   &
     		      sigma3d_db,				       &
     		      real(Etheta),imag(Etheta),real(Ephi),imag(Ephi)
     	     end if
1010         format(4E14.6)
1020         format(8E14.5)

! next phi
           end do
! next theta
         end do
! next frequency
       end do
       
! close output files       
       close(unit=10)
       
       close(unit=11)
       if (thin_layer.eq.1) then
         close(unit=12)
         close(unit=21)
         close(unit=22)
       end if
       
       close(unit=30)
       close(unit=31)
             
       open(unit=99,file='progress')
       write(99,'(A)')'FINISHED: GGI_RCS_sphere'
       close(unit=99)

       END
