from paraview.simple import *
from math import *
import os
import sys
import time
import glob

# Add the mesh boundary 
bmesh_file=glob.glob("*.b_mesh.vtk.0")

name="Outer boundary"

print "Processing file" , bmesh_file , " name=" , name

bmesh=LegacyVTKReader(guiName=name,FileNames=bmesh_file) 

sDataRepresentation= GetDisplayProperties( bmesh ) 
sDataRepresentation.Opacity = 1.0
sDataRepresentation.Representation = 'Wireframe'

sDataRepresentation.DiffuseColor = [1.0, 1.0, 1.0]	

Show( bmesh )   

# Add the volume mesh for each volume material

volume_material_mesh_files= glob.glob("./*.v_mesh.vtk.*")
nvol=len(volume_material_mesh_files)

print "Number of volumes=" , nvol

# loop over the volumes and add gto the image 

for vol in range (0,nvol):

	number=vol+1
	number_string=str(number)

	vmesh_file=glob.glob("*.v_mesh.vtk." +number_string)

        name="volume_"+number_string
	
	print "Processing file" , vmesh_file , " name=" , name

	vmesh=LegacyVTKReader(guiName=name,FileNames=vmesh_file) 
	
	sDataRepresentation= GetDisplayProperties( vmesh ) 
	sDataRepresentation.Opacity = 0.7
	sDataRepresentation.Representation = 'Surface'

	colour=float(vol+1)/float(nvol)

	sDataRepresentation.DiffuseColor = [0.0, 0.0 , colour ]	

	Show( vmesh )   

# Add the surface mesh for each surface material

surface_material_mesh_files= glob.glob("./*.s_mesh.vtk.*")
nsurf=len(surface_material_mesh_files)

print "Number of surfaces=" , nsurf

# loop over the surfaces and add to the image 

for surf in range (0,nsurf):

	number=surf+1
	number_string=str(number)

	smesh_file=glob.glob("*.s_mesh.vtk." +number_string)

        name="surface_"+number_string
	
	print "Processing file" , smesh_file , " name=" , name

	smesh=LegacyVTKReader(guiName=name,FileNames=smesh_file) 
	
	sDataRepresentation= GetDisplayProperties( smesh ) 
	sDataRepresentation.Opacity = 0.7
	sDataRepresentation.Representation = 'Surface'

	colour=float(surf+1)/float(nsurf)

	sDataRepresentation.DiffuseColor = [1.0, colour, 0.0]	

	Show( smesh )   
	

# Add the line meshetry

line_mesh_files= glob.glob("./*.l_mesh.vtk.*")
nlines=len(line_mesh_files)

print "Number of lines=" , nlines

# loop over the lines and add to the image

for line in range (0,nlines):

	number=line+1
	number_string=str(number)

	lmesh_file=glob.glob("*.l_mesh.vtk." +number_string)

        name="line_"+number_string
	
	print "Processing file" , lmesh_file , " name=" , name

	lmesh=LegacyVTKReader(guiName=name,FileNames=lmesh_file) 
	
	sDataRepresentation= GetDisplayProperties( lmesh ) 
	sDataRepresentation.Opacity = 1.0
	sDataRepresentation.Representation = 'Surface'

	colour=float(line+1)/float(nlines)

	sDataRepresentation.DiffuseColor = [ colour, 0.0, 0.0]	

	Show( lmesh )   
	

# Add the point meshetry

point_mesh_files= glob.glob("./*.p_mesh.vtk.*")
npoints=len(point_mesh_files)

print "Number of points=" , npoints

# loop over the points and add to the image

for point in range (0,npoints):

	number=point+1
	number_string=str(number)

	pmesh_file=glob.glob("*.p_mesh.vtk." +number_string)

        name="point_"+number_string
	
	print "Processing file" , pmesh_file , " name=" , name

	pmesh=LegacyVTKReader(guiName=name,FileNames=pmesh_file) 
	
	sDataRepresentation= GetDisplayProperties( pmesh ) 
	sDataRepresentation.Opacity = 1.0
	sDataRepresentation.Representation = 'Surface'

	colour=float(point+1)/float(npoints)

	sDataRepresentation.DiffuseColor = [ colour, 0.0, 0.0]	

	Show( pmesh )   
	

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

