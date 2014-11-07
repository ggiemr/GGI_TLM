#  WRITE THE CABLE INFORMATION

from PyRTF import *
from general_subroutines import *

def Build_cable_section( name , cable_info_file_contents , doc ) :

	ss      = doc.StyleSheet

	cable_section = Section(break_type=3)

	p = Paragraph( ss.ParagraphStyles.Heading2 )
	p.append( 'Cable information')
	cable_section.append( p )

	start_line="#START OF CABLE SUMMARY INFORMATION"
	end_line="#END OF CABLE SUMMARY INFORMATION"

	output_line_list=return_intermediate_lines( start_line , end_line , cable_info_file_contents )

	for line in output_line_list:
		p = Paragraph( ss.ParagraphStyles.Normal )
		p.append( line )
		cable_section.append( p )

	doc.Sections.append( cable_section )
	
	return doc
