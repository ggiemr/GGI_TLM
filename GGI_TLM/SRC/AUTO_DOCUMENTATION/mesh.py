#  WRITE THE MESHING SECTION

from paraview.simple import *
from math import *
import os
import sys
import time
import glob
from PyRTF import *
from general_subroutines import *

def Build_mesh_section( name , info_file_contents , doc ) :

	ss      = doc.StyleSheet

	mesh_section = Section(break_type=3)

	p = Paragraph( ss.ParagraphStyles.Heading2 )
	p.append( 'Mesh')
	mesh_section.append( p )

#  include the meshing information and mesh image(s)

	start_line="#START OF MESH DESCRIPTION"
	end_line="#END OF MESH DESCRIPTION"

	output_line_list=return_intermediate_lines( start_line , end_line , info_file_contents )

	for line in output_line_list:
		p = Paragraph( ss.ParagraphStyles.Normal )
		p.append( line )
		mesh_section.append( p )

#  plot the mesh files using calls to paraview

	boundary_mesh_files= glob.glob("./*b_mesh.vtk.*")
	volume_mesh_files= glob.glob("./*v_mesh.vtk.*")
	surface_mesh_files= glob.glob("./*s_mesh.vtk.*")
	line_mesh_files= glob.glob("./*l_mesh.vtk.*")


	for b_mesh_file in boundary_mesh_files:
		print b_mesh_file
		boundary = OpenDataFile(b_mesh_file)
		
	RenderView1 = GetRenderView()
	dp1 = Show()
	dp1.Representation = 'Wireframe'
	
	RenderView1.Background = [0.0, 0.3333333333333333, 1.0]

	Show(boundary)

	vmeshlist=[]
	vDataRepresentationlist=[]

	nvol=len(volume_mesh_files)
	vol=0

	for v_mesh_file in volume_mesh_files:
		print v_mesh_file

		vmeshlist.append( OpenDataFile(v_mesh_file) )

		vDataRepresentationlist.append( GetDisplayProperties( vmeshlist[vol] ) )
		vDataRepresentationlist[vol].Opacity = 0.5
		vDataRepresentationlist[vol].Representation = 'Surface With Edges'
	
		colour=float(vol+1)/float(nvol)
	
		vDataRepresentationlist[vol].DiffuseColor = [colour, 0.0, 0.0]	
	
		Show( vmeshlist[vol] )	
		vol=vol+1
	


	smeshlist=[]
	sDataRepresentationlist=[]

	nsurf=len(surface_mesh_files)
	surf=0

	for s_mesh_file in surface_mesh_files:
		print s_mesh_file

		smeshlist.append( OpenDataFile(s_mesh_file) )

		sDataRepresentationlist.append( GetDisplayProperties( smeshlist[surf] ) )
		sDataRepresentationlist[surf].Opacity = 0.5
		sDataRepresentationlist[surf].Representation = 'Surface With Edges'
	
		colour=float(surf+1)/float(nsurf)
	
		sDataRepresentationlist[surf].DiffuseColor = [0.0, colour, 0.0]	
	
		Show( smeshlist[surf] )	
		surf=surf+1
	
	lmeshlist=[]
	lDataRepresentationlist=[]

	nline=len(line_mesh_files)
	line=0

	for l_mesh_file in line_mesh_files:
		print l_mesh_file

		lmeshlist.append( OpenDataFile(l_mesh_file) )

		lDataRepresentationlist.append( GetDisplayProperties( lmeshlist[line] ) )
		lDataRepresentationlist[line].Opacity = 0.5
		lDataRepresentationlist[line].Representation = 'Surface With Edges'
	
		colour=float(line+1)/float(nline)
	
		lDataRepresentationlist[line].DiffuseColor = [0.0, 0.0, colour]	
	
		Show( lmeshlist[line] )	
		line=line+1

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

	mesh_image_filename=name+"_1_mesh.png"
	WriteImage(mesh_image_filename)
	image = Image( mesh_image_filename , scale_x=60 , scale_y=60 )
	mesh_section.append( Paragraph( image ) )

# set view angle 2, note the constant 57.29578 converts from degrees to radians
	r=camera.GetDistance()
	theta=70.0/57.29578
	phi=165.0/57.29578
	dx=r*sin(theta)*cos(phi)
	dy=r*sin(theta)*sin(phi)
	dz=r*cos(theta)
	camera.SetPosition(dx,dy,dz)

	Render()

	mesh_image_filename=name+"_2_mesh.png"
	WriteImage(mesh_image_filename)
	image = Image( mesh_image_filename , scale_x=60 , scale_y=60 )
	mesh_section.append( Paragraph( image ) )

# set view angle 3, note the constant 57.29578 converts from degrees to radians
	r=camera.GetDistance()
	theta=120.0/57.29578
	phi=285.0/57.29578
	dx=r*sin(theta)*cos(phi)
	dy=r*sin(theta)*sin(phi)
	dz=r*cos(theta)
	camera.SetPosition(dx,dy,dz)

	Render()

	mesh_image_filename=name+"_3_mesh.png"
	WriteImage(mesh_image_filename)
	image = Image( mesh_image_filename , scale_x=60 , scale_y=60 )
	mesh_section.append( Paragraph( image ) )

# delete the mesh objects from the view

	for vol in range (0,nvol):
		Delete( vmeshlist[vol] )

	for surf in range (0,nsurf):
		Delete( smeshlist[surf] )

	for line in range (0,nline):
		Delete( lmeshlist[line] )

	doc.Sections.append( mesh_section )
	
	return doc
