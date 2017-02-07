from paraview.simple import *
from math import *
import os
import sys
import time
import glob

# Add the surface mesh for each surface material
bmesh_file=glob.glob("*.b_mesh.vtk.0")

name="Outer boundary"

print "Processing file" , bmesh_file , " name=" , name

bmesh=LegacyVTKReader(guiName=name,FileNames=bmesh_file) 

sDataRepresentation= GetDisplayProperties( bmesh ) 
sDataRepresentation.Opacity = 1.0
sDataRepresentation.Representation = 'Wireframe'

sDataRepresentation.DiffuseColor = [1.0, 1.0, 1.0]	

Show( bmesh )   

# Add the surface mesh for each surface material

surface_material_mesh_files= glob.glob("./*.smat_faces.vtk.*")
nsurf=len(surface_material_mesh_files)

print "Number of surface materials=" , nsurf

# loop over the surfaces and add the mesh image 

for surf in range (0,nsurf):

	number=surf+1
	number_string=str(number)

	smesh_file=glob.glob("*.smat_faces.vtk." +number_string)

        name="surface_material_"+number_string
	
	print "Processing file" , smesh_file , " name=" , name

	smesh=LegacyVTKReader(guiName=name,FileNames=smesh_file) 
	
	sDataRepresentation= GetDisplayProperties( smesh ) 
	sDataRepresentation.Opacity = 0.7
	sDataRepresentation.Representation = 'Surface'

	colour=float(surf+1)/float(nsurf)

	sDataRepresentation.DiffuseColor = [1.0, colour, 0.0]	

	Show( smesh )   

# Add the volume mesh for each volume material

volume_material_mesh_files= glob.glob("./*.vmat_cells.vtk.*")
nvol=len(volume_material_mesh_files)

print "Number of volume materials=" , nvol

# loop over the volumes and add the mesh image 

for vol in range (0,nvol):

	number=vol+1
	number_string=str(number)

	vmesh_file=glob.glob("*.vmat_cells.vtk." +number_string)

        name="volume_material_"+number_string
	
	print "Processing file" , vmesh_file , " name=" , name

	vmesh=LegacyVTKReader(guiName=name,FileNames=vmesh_file) 
	
	sDataRepresentation= GetDisplayProperties( vmesh ) 
	sDataRepresentation.Opacity = 0.7
	sDataRepresentation.Representation = 'Surface'

	colour=float(vol+1)/float(nvol)

	sDataRepresentation.DiffuseColor = [0.0, 0.0 , colour ]	

	Show( vmesh )   
	

# Add the volume mesh for each volume material

cable_mesh_files= glob.glob("./*.c_mesh.vtk.*")
ncables=len(cable_mesh_files)

print "Number of cables=" , ncables

# loop over the cables and add the mesh image 

for cable in range (0,ncables):

	number=cable+1
	number_string=str(number)

	cmesh_file=glob.glob("*.c_mesh.vtk." +number_string)

        name="cable_"+number_string
	
	print "Processing file" , cmesh_file , " name=" , name

	cmesh=LegacyVTKReader(guiName=name,FileNames=cmesh_file) 
	
	sDataRepresentation= GetDisplayProperties( cmesh ) 
	sDataRepresentation.Opacity = 1.0
	sDataRepresentation.Representation = 'Surface'

	colour=float(cable+1)/float(ncables)

	sDataRepresentation.DiffuseColor = [ colour, 0.0, 0.0]	

	Show( cmesh )   
	

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

