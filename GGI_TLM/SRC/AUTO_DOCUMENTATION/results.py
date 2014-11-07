
#  write the results section

import os
import sys
import time
import glob
from PyRTF import *
from general_subroutines import *

def Build_results_section( name , info_file_contents , doc ) :

	ss      = doc.StyleSheet

	results_section = Section(break_type=3)

	p = Paragraph( ss.ParagraphStyles.Heading1 )
	p.append( 'Results')
	results_section.append( p )

#  INCLUDE THE RESUTLS .JPG FILES

# run the plot_jpg process to produce the result plots in jpeg format

	os.system("rm *.jpg*")

	sed_command_file = open('./sed_command', 'w')
	sed_command_file.write("{s/#JPG/ /g \n")
	sed_command_file.write("s/#OUTPUT_TO_FILE/ /g } \n")
	sed_command_file.close()

	os.system("sed -f sed_command PROBLEM_SPECIFICATION_FILES/plot_result.plt > plot_result.plt")

	os.system("gnuplot plot_result.plt")

	result_file_list= glob.glob("*.jpg")

	for result_file in result_file_list:

		image = Image( result_file , scale_x=50 , scale_y=50 )
		results_section.append( Paragraph( image ) )

	doc.Sections.append( results_section )
	
	return doc
