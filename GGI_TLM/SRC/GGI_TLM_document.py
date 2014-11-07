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
#     22/10/2014 CJS  Split script into modules - modules are in AUTO_DOCUMENTATION directory
#
#

# put the required imports here...
from paraview.simple import *
from math import *
import os
import sys
import time
import glob
from PyRTF import *

from AUTO_DOCUMENTATION import *

# Open all the files required to construct the document text
# Note that image files are created and read as required for each section. 

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

#  Open the cable.info file to get the required information for the simulation document if it exists

cable_info_file_contents = []
cable_info_filename = name + ".cable_info"
cable_info_filename_list= glob.glob(cable_info_filename)
if len(cable_info_filename_list)==1:
	cable_info_filename = cable_info_filename_list[0]
	cable_info_file = open(cable_info_filename, 'r')
	cable_info_file_contents = cable_info_file.readlines()

#  Open the post procesing info file to get the required information for the simulation document

post_filename_list= glob.glob("GGI_TLM_post_process.info")

# check that there is only one input file	
post_file_contents = []
if len(post_filename_list)==1:
	post_filename = "GGI_TLM_post_process.info"
	post_file = open(post_filename, 'r')
	post_file_contents = post_file.readlines()

#  Start to build the simulation document
DR = Renderer()

# Start to create the document using the PyRTF stuff

doc     = Document()

ss      = doc.StyleSheet

# Build the individual document sections

Build_introduction_section( name , inp_file_contents , doc )
Build_geometry_section( name , info_file_contents , doc  )
Build_mesh_section( name , info_file_contents , doc  )
Build_vmat_section( name , info_file_contents , doc  )
Build_smat_section( name , info_file_contents , doc  )
Build_cable_section( name , cable_info_file_contents , doc  )
Build_run_section( name , info_file_contents , doc  )
Build_output_section( name , info_file_contents , doc  )
Build_post_section( name , post_file_contents , doc  )
Build_results_section( name , info_file_contents , doc  )
Build_appendix1_section( name , inp_file_contents , doc  )
	
# WRITE DOCUMENT TO FILE

document_name=name + ".rtf"

print "Simulation document name is : " + document_name

document_filename = open(document_name, 'w')

DR.Write( doc, document_filename )

document_filename.close()

progress_file = open('./progress', 'w')
progress_file.write("GGI_TLM_document: Finished Correctly")
progress_file.close()
