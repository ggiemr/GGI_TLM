#  WRITE THE OUTPUT INFORMATION

from PyRTF import *
from general_subroutines import *

def Build_output_section( name , info_file_contents , doc ) :

	ss      = doc.StyleSheet

	output_section = Section(break_type=3)

	p = Paragraph( ss.ParagraphStyles.Heading2 )
	p.append( 'Output information')
	output_section.append( p )

	start_line="#START OF OUTPUT INFORMATION"
	end_line="#END OF OUTPUT INFORMATION"

	output_line_list=return_intermediate_lines( start_line , end_line , info_file_contents )

	for line in output_line_list:
		p = Paragraph( ss.ParagraphStyles.Normal )
		p.append( line )
		output_section.append( p )


	doc.Sections.append( output_section )
	
	return doc
