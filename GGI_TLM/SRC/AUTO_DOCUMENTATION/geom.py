
#  WRITE THE GEOMETRY SECTION
from paraview.simple import *
from math import *
import os
import sys
import time
import glob
from PyRTF import *
from general_subroutines import *

def Build_geometry_section( name , info_file_contents , doc ) :

	ss      = doc.StyleSheet

	geom_section = Section(break_type=3)

	p = Paragraph( ss.ParagraphStyles.Heading2 )
	p.append( 'Geometry')
	geom_section.append( p )

#  include the geometry files
#  include the geometry information and geometry image(s)

	start_line="#START OF GEOMETRY DESCRIPTION"
	end_line="#END OF GEOMETRY DESCRIPTION"

	output_line_list=return_intermediate_lines( start_line , end_line , info_file_contents )

	for line in output_line_list:
		p = Paragraph( ss.ParagraphStyles.Normal )
		p.append( line )
		geom_section.append( p )

#  plot the geom files using calls to paraview

	volume_geom_files= glob.glob("./*.v_geom.vtk.*")
	surface_geom_files= glob.glob("./*.s_geom.vtk.*")
	line_geom_files= glob.glob("./*.l_geom.vtk.*")

	boundary_mesh_files= glob.glob("./*.b_mesh.vtk.*")

	for b_mesh_file in boundary_mesh_files:
		print b_mesh_file
		boundary = OpenDataFile(b_mesh_file)
		
	dp1 = Show()
	dp1.Representation = 'Wireframe'

	Show(boundary)

	vgeomlist=[]
	vDataRepresentationlist=[]

	nvol=len(volume_geom_files)
	vol=0

	for v_geom_file in volume_geom_files:
		print v_geom_file

		vgeomlist.append( OpenDataFile(v_geom_file) )

		vDataRepresentationlist.append( GetDisplayProperties( vgeomlist[vol] ) )
		vDataRepresentationlist[vol].Opacity = 0.5
		vDataRepresentationlist[vol].Representation = 'Surface With Edges'
	
		colour=float(vol+1)/float(nvol)
	
		vDataRepresentationlist[vol].DiffuseColor = [colour, 0.0, 0.0]	
	
		Show( vgeomlist[vol] )	
		vol=vol+1

	sgeomlist=[]
	sDataRepresentationlist=[]

	nsurf=len(surface_geom_files)
	surf=0

	for s_geom_file in surface_geom_files:
		print s_geom_file

		sgeomlist.append( OpenDataFile(s_geom_file) )

		sDataRepresentationlist.append( GetDisplayProperties( sgeomlist[surf] ) )
		sDataRepresentationlist[surf].Opacity = 0.5
		sDataRepresentationlist[surf].Representation = 'Surface With Edges'
	
		colour=float(surf+1)/float(nsurf)
	
		sDataRepresentationlist[surf].DiffuseColor = [0.0, colour, 0.0]	
	
		Show( sgeomlist[surf] )	
		surf=surf+1
	
	lgeomlist=[]
	lDataRepresentationlist=[]

	nline=len(line_geom_files)
	line=0

	for l_geom_file in line_geom_files:
		print l_geom_file

		lgeomlist.append( OpenDataFile(l_geom_file) )

		lDataRepresentationlist.append( GetDisplayProperties( lgeomlist[line] ) )
		lDataRepresentationlist[line].Opacity = 0.5
		lDataRepresentationlist[line].Representation = 'Surface With Edges'
	
		colour=float(line+1)/float(nline)
	
		lDataRepresentationlist[line].DiffuseColor = [0.0, 0.0, colour]	
	
		Show( lgeomlist[line] )	
		line=line+1

# set background colour to blue    		
	RenderView1 = GetRenderView()	
	RenderView1.Background = [0.0, 0.3333333333333333, 1.0]

# get camera attributes	    
	camera = GetActiveCamera() 

# set view angle 1, note the constant 57.29578 converts from degrees to radians
	r=camera.GetDistance()
	theta=30.0/57.29578
	phi=45.0/57.29578
	dx=r*sin(theta)*cos(phi)
	dy=r*sin(theta)*sin(phi)
	dz=r*cos(theta)
	camera.SetPosition(dx,dy,dz)
	Render()

	geometry_image_filename=name+"_1_geometry.png"
	WriteImage(geometry_image_filename)

	image = Image( geometry_image_filename , scale_x=60 , scale_y=60 )
	geom_section.append( Paragraph( image ) )

# set view angle 2, note the constant 57.29578 converts from degrees to radians
	r=camera.GetDistance()
	theta=70.0/57.29578
	phi=165.0/57.29578
	dx=r*sin(theta)*cos(phi)
	dy=r*sin(theta)*sin(phi)
	dz=r*cos(theta)
	camera.SetPosition(dx,dy,dz)
	Render()

	geometry_image_filename=name+"_2_geometry.png"
	WriteImage(geometry_image_filename)

	image = Image( geometry_image_filename , scale_x=60 , scale_y=60 )
	geom_section.append( Paragraph( image ) )

# set view angle 3, note the constant 57.29578 converts from degrees to radians
	r=camera.GetDistance()
	theta=120.0/57.29578
	phi=285.0/57.29578
	dx=r*sin(theta)*cos(phi)
	dy=r*sin(theta)*sin(phi)
	dz=r*cos(theta)
	camera.SetPosition(dx,dy,dz)
	Render()

	geometry_image_filename=name+"_3_geometry.png"
	WriteImage(geometry_image_filename)

	image = Image( geometry_image_filename , scale_x=60 , scale_y=60 )
	geom_section.append( Paragraph( image ) )

# delete the geometry objects from the view

	for vol in range (0,nvol):
		Delete( vgeomlist[vol] )

	for surf in range (0,nsurf):
		Delete( sgeomlist[surf] )

	for line in range (0,nline):
		Delete( lgeomlist[line] )

	doc.Sections.append( geom_section )
	
	return doc
