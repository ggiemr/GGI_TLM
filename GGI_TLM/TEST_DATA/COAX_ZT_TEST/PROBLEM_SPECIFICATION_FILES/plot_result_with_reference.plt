#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_29_coax_ZT_test_current_ref.jpg"
set xlabel "Time (s)"
set ylabel "A"
plot "PROBLEM_SPECIFICATION_FILES/coax_original_TLM.cable_current.tout" u 1:3 title "Coax current: Fieldsolve Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/coax.cable_current.tout_ref" u 1:3 title "Coax current: GGI_TLM Reference" w p ,\
     "coax.cable_current.tout" u 1:3 title "Coax current: GGI_TLM" w l
     
#PAUSE pause -1

set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_29_coax_ZT_test_field_ref.jpg"
set xlabel "Time (s)"
set ylabel "A"
plot "PROBLEM_SPECIFICATION_FILES/coax_original_TLM.field.tout" u 1:3 title "Coax field: Fieldsolve Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/coax.field.tout_ref" u 1:3 title "Coax field: GGI_TLM Reference" w p ,\
     "coax.field.tout" u 1:3 title "Coax field: GGI_TLM" w l
     
#PAUSE pause -1
