# INCLUDE THE VOLUME MATERIAL INFORMATION

from paraview.simple import *
from math import *
import os
import sys
import time
import glob
from PyRTF import *
from general_subroutines import *

def Build_vmat_section( name , info_file_contents , doc ) :

	ss      = doc.StyleSheet

	vmat_section = Section(break_type=3)

	p = Paragraph( ss.ParagraphStyles.Heading2 )
	p.append( 'Volume materials')
	vmat_section.append( p )

#  include the meshing information and mesh image(s)

	start_line="#START OF VOLUME MATERIAL DESCRIPTION"
	end_line="#END OF VOLUME MATERIAL DESCRIPTION"

	output_line_list=return_intermediate_lines( start_line , end_line , info_file_contents )

	for line in output_line_list:
		p = Paragraph( ss.ParagraphStyles.Normal )
		p.append( line )
		vmat_section.append( p )


# run the GGI_TLM_material_model_checks process to get all the volume material frequency responses in jpeg format

	os.system("rm *.jpg*")

	GGI_TLM_material_model_checks_file = open('./GGI_TLM_material_model_checks_in.txt', 'w')
	GGI_TLM_material_model_checks_file.write(name + "\n")
	GGI_TLM_material_model_checks_file.write("5   : option to output all volume material data to file \n")
	GGI_TLM_material_model_checks_file.write("0   : exit \n")
	GGI_TLM_material_model_checks_file.close()

	os.system("GGI_TLM_material_model_checks < GGI_TLM_material_model_checks_in.txt")

# Put the volume mesh for each volume material in the document along with the associated 
# frequency response plots

	volume_material_mesh_files= glob.glob("./*.vmat_cells.vtk.*")
	nvol=len(volume_material_mesh_files)

	volume_epsr_response_file_list= glob.glob("*eps.jpg.*")
	volume_mur_response_file_list= glob.glob("*mu.jpg.*")

# check that we have the correct number of files
	if len(volume_epsr_response_file_list) != nvol:
		print "ERROR: There is a problem with the number of files relating to volume materials"
		print "nvol=" +str(nvol) + " n_epsr files="+ str(len(volume_epsr_response_file_list) )
		quit
	
	if len(volume_mur_response_file_list) != nvol:
		print "ERROR: There is a problem with the number of files relating to volume materials"
		print "nvol=" +str(nvol) + " n_mur files="+ str(len(volume_mur_response_file_list) )
		quit

# loop over the volumes and add the mesh image, premittivity and permeability responses for each
# volume material 

	for vol in range (0,nvol):

		number=vol+1
		number_string=str(number)
	
		vmesh_file=glob.glob("*.vmat_cells.vtk." +number_string)
	
# First plot the mesh to file	
	        vmesh = OpenDataFile(vmesh_file)
		vDataRepresentation = GetDisplayProperties( vmesh )
		vDataRepresentation.Opacity = 0.5
		vDataRepresentation.Representation = 'Surface With Edges'
	
		colour=float(vol+1)/float(nvol)
	
		vDataRepresentation.DiffuseColor = [ colour, 0.0, 0.0]	
	
		Show( vmesh )	
		vol=vol+1
	
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

		mesh_image_filename=name+"_volume_material_mesh.png"
		WriteImage(mesh_image_filename)
		image = Image( mesh_image_filename , scale_x=60 , scale_y=60 )
		vmat_section.append( Paragraph( image ) )
	
		Delete( vmesh )
	
# put the volume material frequency responses into the document

		base_filename_list=[ "*eps.jpg." , "*mu.jpg." ]
	
		for filename in base_filename_list:

			imagefile =glob.glob(filename+number_string)
# move to temporary .jpg file (PyRTF expects jpeg files to have the extension .jpg not .jpg.n) then add to document	
# note we are assuming that there is only a single file matching the pattern here...
			if len(imagefile) == 1:
				os.system("mv " + imagefile[0] + " temp_eps.jpg")	
				image = Image( "temp_eps.jpg" , scale_x=30 , scale_y=30 )
				vmat_section.append( Paragraph( image ) )

	doc.Sections.append( vmat_section )
	
	return doc
