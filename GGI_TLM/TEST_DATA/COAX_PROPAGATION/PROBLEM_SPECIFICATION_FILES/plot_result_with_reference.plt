#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_43_COAX_PROPAGATION_ref.jpg"

set title "Cable current"

set xlabel "Time (s)"
set ylabel "A"
plot "PROBLEM_SPECIFICATION_FILES/coax.cable_current.tout_ref" u 1:3 title "GGI_TLM reference" w p,\
     "coax.cable_current.tout" u 1:3 title "GGI_TLM" w p
#PAUSE pause -1

