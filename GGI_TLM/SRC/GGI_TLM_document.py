#
#    GGI_TLM Time domain electromagnetic field solver based on the TLM method
#    Copyright (C) 2014  Chris Smartt
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.   
#
#
# NAME
#     GGI_TLM_document
#
# DESCRIPTION
#     Python code to automatically create GGI_TLM simulation documents
#     run with the command:
#     pvpython GGI_TLM_simulation_document.py
#    
#     The code uses the PyRTF code by Simon Cusack which is under the GPL license
#     
#     
# COMMENTS
#     Makes calls to paraview to get the png files for geometry and mesh
#     Call to gnuplot to generate results files
#     Information is derived from the .inp and .info files
#     Still missing the following information:
#           Cable information
#           Post processing information
#
#     We also need to sort out the format for the runtime memory usage etc
#
# HISTORY
#
#     started 10/10/2014 CJS
#     15/10/2014 CJS  Some improvement to the quality of the figures - add some colour, change the view to off axis etc
#
#


def return_intermediate_lines( start_line , end_line , contents ) :

	output_line_list=[]

	include_line = False

	for line in contents:
		
		if ( start_line in line ):
#			print( "Start_line found:" + line )
			include_line = True
		if ( end_line in line ):
#			print( "End_line found:" + line )
			include_line = False
	
		if ((include_line) & ( not start_line in line)):
			output_line_list.append(line)

	return output_line_list

# put the required imports here...
from paraview.simple import *
from math import *
import os
import sys
import time
import glob
from PyRTF import *


#  get the problem name from the input file

input_filename_list= glob.glob("*.inp")

# check that there is only one input file
if len(input_filename_list)!=1:
	print "Error: There should be one and only one input file in the directory..."
	progress_file = open('./progress', 'w')
	progress_file.write("Error: There should be one and only one input file in the directory...")
	progress_file.close()
	exit
	
input_filename=input_filename_list[0]
namelen=len(input_filename)

name=input_filename[0:namelen-4]

print "Problem name is : " + name

#  Open the input file to get the required information for the simulation document

inp_filename = name + ".inp"
inp_file = open(inp_filename, 'r')

inp_file_contents = inp_file.readlines()

#  Open the info file to get the required information for the simulation document

info_filename = name + ".info"
info_file = open(info_filename, 'r')

info_file_contents = info_file.readlines()

#  Start to build the simulation document
DR = Renderer()

# Start to create the document using the PyRTF stuff

doc     = Document()
ss      = doc.StyleSheet
section = Section()
doc.Sections.append( section )

title = 'Simulation name: ' + name

p = Paragraph( ss.ParagraphStyles.Heading1 )
p.append( title )
section.append( p )

p = Paragraph( ss.ParagraphStyles.Heading2 )
p.append( 'Model description')
section.append( p )

#  WRITE THE PROBLEM DESCRIPTION SECTION

start_line="#START OF DESCRIPTION"
end_line="#END OF DESCRIPTION"

output_line_list=return_intermediate_lines( start_line , end_line , inp_file_contents )

for line in output_line_list:
	p = Paragraph( ss.ParagraphStyles.Normal )
	p.append( line )
	section.append( p )

#  WRITE THE GEOMETRY SECTION

p = Paragraph( ss.ParagraphStyles.Heading2 )
p.append( 'Geometry')
section.append( p )

#  include the geometry files
#  include the geometry information and geometry image(s)

start_line="#START OF GEOMETRY DESCRIPTION"
end_line="#END OF GEOMETRY DESCRIPTION"

output_line_list=return_intermediate_lines( start_line , end_line , info_file_contents )

for line in output_line_list:
	p = Paragraph( ss.ParagraphStyles.Normal )
	p.append( line )
	section.append( p )

#  plot the geom files using calls to paraview

volume_geom_files= glob.glob("./*v_geom.vtk.*")
surface_geom_files= glob.glob("./*s_geom.vtk.*")
line_geom_files= glob.glob("./*l_geom.vtk.*")

boundary_mesh_files= glob.glob("./*b_mesh.vtk.*")

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

	lDataRepresentationlist.append( GetDisplayProperties( sgeomlist[line] ) )
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
section.append( Paragraph( image ) )

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
section.append( Paragraph( image ) )

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
section.append( Paragraph( image ) )

# delete the geometry objects from the view

for vol in range (0,nvol):
	Delete( vgeomlist[vol] )

for surf in range (0,nsurf):
	Delete( sgeomlist[surf] )

for line in range (0,nline):
	Delete( lgeomlist[line] )


#  WRITE THE MESHING SECTION

p = Paragraph( ss.ParagraphStyles.Heading2 )
p.append( 'Mesh')
section.append( p )

#  include the meshing information and mesh image(s)

start_line="#START OF MESH DESCRIPTION"
end_line="#END OF MESH DESCRIPTION"

output_line_list=return_intermediate_lines( start_line , end_line , info_file_contents )

for line in output_line_list:
	p = Paragraph( ss.ParagraphStyles.Normal )
	p.append( line )
	section.append( p )

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

	lDataRepresentationlist.append( GetDisplayProperties( smeshlist[line] ) )
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

mesh_image_filename=name+"_2_mesh.png"
WriteImage(mesh_image_filename)
image = Image( mesh_image_filename , scale_x=60 , scale_y=60 )
section.append( Paragraph( image ) )

# set view angle 2, note the constant 57.29578 converts from degrees to radians
r=camera.GetDistance()
theta=70.0/57.29578
phi=165.0/57.29578
dx=r*sin(theta)*cos(phi)
dy=r*sin(theta)*sin(phi)
dz=r*cos(theta)
camera.SetPosition(dx,dy,dz)

Render()

mesh_image_filename=name+"_1_mesh.png"
WriteImage(mesh_image_filename)
image = Image( mesh_image_filename , scale_x=60 , scale_y=60 )
section.append( Paragraph( image ) )

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
section.append( Paragraph( image ) )

# INCLUDE THE VOLUME MATERIAL INFORMATION

p = Paragraph( ss.ParagraphStyles.Heading2 )
p.append( 'Volume materials')
section.append( p )

#  include the meshing information and mesh image(s)

start_line="#START OF VOLUME MATERIAL DESCRIPTION"
end_line="#END OF VOLUME MATERIAL DESCRIPTION"

output_line_list=return_intermediate_lines( start_line , end_line , info_file_contents )

for line in output_line_list:
	p = Paragraph( ss.ParagraphStyles.Normal )
	p.append( line )
	section.append( p )


# run the GGI_TLM_material_model_checks process to get all the volume material frequency responses in jpeg format

os.system("rm *.jpg")

GGI_TLM_material_model_checks_file = open('./GGI_TLM_material_model_checks_in.txt', 'w')
GGI_TLM_material_model_checks_file.write(name + "\n")
GGI_TLM_material_model_checks_file.write("1   : option to output volume material frequency responses \n")
GGI_TLM_material_model_checks_file.write("0   : option to output responses for all volumes \n")
GGI_TLM_material_model_checks_file.write("0   : exit \n")
GGI_TLM_material_model_checks_file.close()

os.system("GGI_TLM_material_model_checks < GGI_TLM_material_model_checks_in.txt")

result_file_list= glob.glob("*.jpg")

for result_file in result_file_list:

	image = Image( result_file , scale_x=30 , scale_y=30 )
	section.append( Paragraph( image ) )

# INCLUDE THE SURFACE MATERIAL INFORMATION

p = Paragraph( ss.ParagraphStyles.Heading2 )
p.append( 'Surface materials')
section.append( p )

#  include the meshing information and mesh image(s)

start_line="#START OF SURFACE MATERIAL DESCRIPTION"
end_line="#END OF SURFACE MATERIAL DESCRIPTION"

output_line_list=return_intermediate_lines( start_line , end_line , info_file_contents )

for line in output_line_list:
	p = Paragraph( ss.ParagraphStyles.Normal )
	p.append( line )
	section.append( p )

# run the GGI_TLM_material_model_checks process to get all the volume material frequency responses in jpeg format

os.system("rm *.jpg")

GGI_TLM_material_model_checks_file = open('./GGI_TLM_material_model_checks_in.txt', 'w')
GGI_TLM_material_model_checks_file.write(name + "\n")
GGI_TLM_material_model_checks_file.write("3   : option to output surface material frequency responses \n")
GGI_TLM_material_model_checks_file.write("0   : option to output responses for all surfaces \n")
GGI_TLM_material_model_checks_file.write("0   : exit \n")
GGI_TLM_material_model_checks_file.close()

os.system("GGI_TLM_material_model_checks < GGI_TLM_material_model_checks_in.txt")

result_file_list= glob.glob("*.jpg")

for result_file in result_file_list:

	image = Image( result_file , scale_x=30 , scale_y=30 )
	section.append( Paragraph( image ) )

#  WRITE THE RUNTIME INFORMATION

p = Paragraph( ss.ParagraphStyles.Heading2 )
p.append( 'Run information')
section.append( p )

start_line="#START OF RUN INFORMATION"
end_line="#END OF RUN INFORMATION"

output_line_list=return_intermediate_lines( start_line , end_line , info_file_contents )

for line in output_line_list:
	p = Paragraph( ss.ParagraphStyles.Normal )
	p.append( line )
	section.append( p )

##  WRITE THE POST PROCESSING INFORMATION
#
#p = Paragraph( ss.ParagraphStyles.Heading2 )
#p.append( 'Post Processing')
#section.append( p )
#
#p = Paragraph( ss.ParagraphStyles.Normal )
#p.append( 'text....' )
#section.append( p )

#  write the results section

p = Paragraph( ss.ParagraphStyles.Heading1 )
p.append( 'Results')
section.append( p )

#  INCLUDE THE RESUTLS .JPG FILES

# run the plot_jpg process to produce the result plots in jpeg format

os.system("rm *.jpg")

sed_command_file = open('./sed_command', 'w')
sed_command_file.write("{s/#JPG/ /g \n")
sed_command_file.write("s/#OUTPUT_TO_FILE/ /g } \n")
sed_command_file.close()

os.system("sed -f sed_command PROBLEM_SPECIFICATION_FILES/plot_result.plt > plot_result.plt")

os.system("gnuplot plot_result.plt")

result_file_list= glob.glob("*.jpg")

for result_file in result_file_list:

	image = Image( result_file , scale_x=50 , scale_y=50 )
	section.append( Paragraph( image ) )

# ADD THE INPUT FILE AS AN APPENDIX

p = Paragraph( ss.ParagraphStyles.Heading1 )
p.append( 'APPENDIX 1: GGI_TLM input file')
section.append( p )

for line in inp_file_contents:
	p = Paragraph( ss.ParagraphStyles.Normal )
	p.append( line )
	section.append( p )

# WRITE DOCUMENT TO FILE

document_name=name + ".rtf"

print "Simulation document name is : " + document_name

document_filename = open(document_name, 'w')

DR.Write( doc, document_filename )

document_filename.close()

progress_file = open('./progress', 'w')
progress_file.write("GGI_TLM_document: Finished Correctly")
progress_file.close()
