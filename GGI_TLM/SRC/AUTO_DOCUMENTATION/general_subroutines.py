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
#     general_subroutines.py
#     GGI_TLM_document
#
# DESCRIPTION
#     General subroutines for automatic document generation
#     The code uses the PyRTF code by Simon Cusack which is under the GPL license
#     
#	Includes the subroutines:
#		return_intermediate_lines
#     
# COMMENTS
#     
#
# HISTORY
#
#     started 10/10/2014 CJS
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

