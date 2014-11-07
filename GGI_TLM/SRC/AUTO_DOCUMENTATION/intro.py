from PyRTF import *
from general_subroutines import *

def Build_introduction_section( name , inp_file_contents , doc ) :

	ss      = doc.StyleSheet

	intro_section = Section(break_type=None)

	title = 'Simulation name: ' + name

	p = Paragraph( ss.ParagraphStyles.Heading1 )
	p.append( title )
	intro_section.append( p )

	p = Paragraph( ss.ParagraphStyles.Heading2 )
	p.append( 'Model description')
	intro_section.append( p )

#  WRITE THE PROBLEM DESCRIPTION SECTION

	start_line="#START OF DESCRIPTION"
	end_line="#END OF DESCRIPTION"

	output_line_list=return_intermediate_lines( start_line , end_line , inp_file_contents )

	for line in output_line_list:
		p = Paragraph( ss.ParagraphStyles.Normal )
		p.append( line )
		intro_section.append( p )
		
	doc.Sections.append( intro_section )
	
	return doc
