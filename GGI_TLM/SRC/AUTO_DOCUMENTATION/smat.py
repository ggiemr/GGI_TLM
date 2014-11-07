# INCLUDE THE SURFACE MATERIAL INFORMATION

from paraview.simple import *
from math import *
import os
import sys
import time
import glob
from PyRTF import *
from general_subroutines import *

def Build_smat_section( name , info_file_contents , doc ) :

	ss      = doc.StyleSheet

	smat_section = Section(break_type=3)

	p = Paragraph( ss.ParagraphStyles.Heading2 )
	p.append( 'Surface materials')
	smat_section.append( p )

#  include the meshing information and mesh image(s)

	start_line="#START OF SURFACE MATERIAL DESCRIPTION"
	end_line="#END OF SURFACE MATERIAL DESCRIPTION"

	output_line_list=return_intermediate_lines( start_line , end_line , info_file_contents )

	for line in output_line_list:
		p = Paragraph( ss.ParagraphStyles.Normal )
		p.append( line )
		smat_section.append( p )

# run the GGI_TLM_material_model_checks process to get all the volume material frequency responses in jpeg format

	os.system("rm *.jpg*")

	GGI_TLM_material_model_checks_file = open('./GGI_TLM_material_model_checks_in.txt', 'w')
	GGI_TLM_material_model_checks_file.write(name + "\n")
	GGI_TLM_material_model_checks_file.write("6   : option to output all surface material data to file \n")
	GGI_TLM_material_model_checks_file.write("0   : exit \n")
	GGI_TLM_material_model_checks_file.close()

	os.system("GGI_TLM_material_model_checks < GGI_TLM_material_model_checks_in.txt")

# Put the surface mesh for each surface material in the document

	surface_material_mesh_files= glob.glob("./*.smat_faces.vtk.*")
	nsurf=len(surface_material_mesh_files)

# loop over the surfaces and add the mesh image and IBC matrix frequency responses for each

	for surf in range (0,nsurf):

		number=surf+1
		number_string=str(number)
	
		smesh_file=glob.glob("*.smat_faces.vtk." +number_string)

		smesh=OpenDataFile(smesh_file) 

		sDataRepresentation= GetDisplayProperties( smesh ) 
		sDataRepresentation.Opacity = 0.5
		sDataRepresentation.Representation = 'Surface With Edges'
	
		colour=float(surf+1)/float(nsurf)
	
		sDataRepresentation.DiffuseColor = [0.0, colour, 0.0]	
	
		Show( smesh )	

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

		mesh_image_filename=name+"_surface_material_mesh.png"
		WriteImage(mesh_image_filename)
		image = Image( mesh_image_filename , scale_x=60 , scale_y=60 )
		smat_section.append( Paragraph( image ) )

# delete the mesh objects from the view

		Delete( smesh )

# put the surface material frequency responses into the document

		anisotropic_test=glob.glob("*:X_polarisation_z11.jpg."+number_string)
	
		if len(anisotropic_test) == 0:
# Assume we have a normal isotropic IBC model
			base_filename_list=[ "*_z11.jpg." , 
			                     "*_z12.jpg." ,
			                     "*_z21.jpg." ,
			                     "*_z22.jpg." ]

	        else:
# Assume we have an aisotropic IBC model	
			base_filename_list=[ "*:X_polarisation_z11.jpg." , 
			                     "*:X_polarisation_z12.jpg." ,
			                     "*:X_polarisation_z21.jpg." ,
			                     "*:X_polarisation_z22.jpg." ,
		                             "*:Y_polarisation_z11.jpg." , 
			                     "*:Y_polarisation_z12.jpg." ,
			                     "*:Y_polarisation_z21.jpg." ,
			                     "*:Y_polarisation_z22.jpg." ,
		                             "*:Z_polarisation_z11.jpg." , 
			                     "*:Z_polarisation_z12.jpg." ,
			                     "*:Z_polarisation_z21.jpg." ,
			                     "*:Z_polarisation_z22.jpg." ]
	
		for filename in base_filename_list:
			imagefile =glob.glob(filename+number_string)
# move to temporary .jpg file (PyRTF expects jpeg files to have the extension .jpg not .jpg.n) then add to document	
# note we are assuming that there is only a single file matching the pattern here...
			if len(imagefile) == 1:
				os.system("mv " + imagefile[0] + " temp_IBC.jpg")	
				image = Image( "temp_IBC.jpg" , scale_x=30 , scale_y=30 )
				smat_section.append( Paragraph( image ) )

	doc.Sections.append( smat_section )
	
	return doc

