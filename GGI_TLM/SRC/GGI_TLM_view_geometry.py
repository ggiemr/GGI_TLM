from paraview.simple import *
from math import *
import os
import sys
import time
import glob

# Add the mesh boundary 
bgeom_file=glob.glob("*.b_mesh.vtk.0")

name="Outer boundary"

print "Processing file" , bgeom_file , " name=" , name

bgeom=LegacyVTKReader(guiName=name,FileNames=bgeom_file) 

sDataRepresentation= GetDisplayProperties( bgeom ) 
sDataRepresentation.Opacity = 1.0
sDataRepresentation.Representation = 'Wireframe'

sDataRepresentation.DiffuseColor = [1.0, 1.0, 1.0]	

Show( bgeom )   

# Add the volume geom for each volume material

volume_material_geom_files= glob.glob("./*.v_geom.vtk.*")
nvol=len(volume_material_geom_files)

print "Number of volumes=" , nvol

# loop over the volumes and add gto the image 

for vol in range (0,nvol):

	number=vol+1
	number_string=str(number)

	vgeom_file=glob.glob("*.v_geom.vtk." +number_string)

        name="volume_"+number_string
	
	print "Processing file" , vgeom_file , " name=" , name

	vgeom=LegacyVTKReader(guiName=name,FileNames=vgeom_file) 
	
	sDataRepresentation= GetDisplayProperties( vgeom ) 
	sDataRepresentation.Opacity = 0.7
	sDataRepresentation.Representation = 'Surface'

	colour=float(vol+1)/float(nvol)

	sDataRepresentation.DiffuseColor = [0.0, 0.0 , colour ]	

	Show( vgeom )   

# Add the surface geom for each surface material

surface_material_geom_files= glob.glob("./*.s_geom.vtk.*")
nsurf=len(surface_material_geom_files)

print "Number of surfaces=" , nsurf

# loop over the surfaces and add to the image 

for surf in range (0,nsurf):

	number=surf+1
	number_string=str(number)

	sgeom_file=glob.glob("*.s_geom.vtk." +number_string)

        name="surface_"+number_string
	
	print "Processing file" , sgeom_file , " name=" , name

	sgeom=LegacyVTKReader(guiName=name,FileNames=sgeom_file) 
	
	sDataRepresentation= GetDisplayProperties( sgeom ) 
	sDataRepresentation.Opacity = 0.7
	sDataRepresentation.Representation = 'Surface'

	colour=float(surf+1)/float(nsurf)

	sDataRepresentation.DiffuseColor = [1.0, colour, 0.0]	

	Show( sgeom )   
	

# Add the line geometry

line_geom_files= glob.glob("./*.l_geom.vtk.*")
nlines=len(line_geom_files)

print "Number of lines=" , nlines

# loop over the lines and add to the image

for line in range (0,nlines):

	number=line+1
	number_string=str(number)

	lgeom_file=glob.glob("*.l_geom.vtk." +number_string)

        name="line_"+number_string
	
	print "Processing file" , lgeom_file , " name=" , name

	lgeom=LegacyVTKReader(guiName=name,FileNames=lgeom_file) 
	
	sDataRepresentation= GetDisplayProperties( lgeom ) 
	sDataRepresentation.Opacity = 1.0
	sDataRepresentation.Representation = 'Surface'

	colour=float(line+1)/float(nlines)

	sDataRepresentation.DiffuseColor = [ colour, 0.0, 0.0]	

	Show( lgeom )   
	

# Add the point geometry

point_geom_files= glob.glob("./*.p_geom.vtk.*")
npoints=len(point_geom_files)

print "Number of points=" , npoints

# loop over the points and add to the image

for point in range (0,npoints):

	number=point+1
	number_string=str(number)

	pgeom_file=glob.glob("*.p_geom.vtk." +number_string)

        name="point_"+number_string
	
	print "Processing file" , pgeom_file , " name=" , name

	pgeom=LegacyVTKReader(guiName=name,FileNames=pgeom_file) 
	
	sDataRepresentation= GetDisplayProperties( pgeom ) 
	sDataRepresentation.Opacity = 1.0
	sDataRepresentation.Representation = 'Surface'

	colour=float(point+1)/float(npoints)

	sDataRepresentation.DiffuseColor = [ colour, 0.0, 0.0]	

	Show( pgeom )   
	

# set background colour to blue    		
RenderView1 = GetRenderView()   
RenderView1.Background = [0.0, 0.3333333333333333, 1.0]

# get camera attributes	    
camera = GetActiveCamera() 

# set view angle 1, note the constant 57.29578 converts from degrees to radians
r=camera.GetDistance()
theta=30.0/57.29578
phi=60.0/57.29578
dx=r*sin(theta)*cos(phi)
dy=r*sin(theta)*sin(phi)
dz=r*cos(theta)
camera.SetPosition(dx,dy,dz)

Render()

