
# ADD THE INPUT FILE AS AN APPENDIX

from PyRTF import *
from general_subroutines import *

def Build_appendix1_section( name , inp_file_contents , doc ) :

	ss      = doc.StyleSheet

	appendix1_section = Section(break_type=3)

	p = Paragraph( ss.ParagraphStyles.Heading1 )
	p.append( 'APPENDIX 1: GGI_TLM input file')
	appendix1_section.append( p )

	for line in inp_file_contents:
		p = Paragraph( ss.ParagraphStyles.Normal )
		p.append( line )
		appendix1_section.append( p )

	doc.Sections.append( appendix1_section )
	
	return doc
