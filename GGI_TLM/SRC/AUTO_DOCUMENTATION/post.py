#  WRITE THE POST PROCESSING INFORMATION

from PyRTF import *
from general_subroutines import *

def Build_post_section( name , post_file_contents , doc ) :

	ss      = doc.StyleSheet

	post_section = Section(break_type=3)

	p = Paragraph( ss.ParagraphStyles.Heading2 )
	p.append( 'Post Processing')
	post_section.append( p )

	start_line="#START OF POST_PROCESSING DESCRIPTION"
	end_line="#END OF POST_PROCESSING DESCRIPTION"

	output_line_list=return_intermediate_lines( start_line , end_line , post_file_contents )

	for line in output_line_list:
		p = Paragraph( ss.ParagraphStyles.Normal )
		p.append( line )
		post_section.append( p )


	doc.Sections.append( post_section )
	
	return doc
