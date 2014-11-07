
#  WRITE THE RUNTIME INFORMATION

from PyRTF import *
from general_subroutines import *

def Build_run_section( name , info_file_contents , doc ) :

	ss      = doc.StyleSheet

	run_section = Section(break_type=3)

	p = Paragraph( ss.ParagraphStyles.Heading2 )
	p.append( 'Run information')
	run_section.append( p )

	start_line="#START OF RUN INFORMATION"
	end_line="#END OF RUN INFORMATION"

	output_line_list=return_intermediate_lines( start_line , end_line , info_file_contents )

	for line in output_line_list:
		p = Paragraph( ss.ParagraphStyles.Normal )
		p.append( line )
		run_section.append( p )


	doc.Sections.append( run_section )
	
	return doc
