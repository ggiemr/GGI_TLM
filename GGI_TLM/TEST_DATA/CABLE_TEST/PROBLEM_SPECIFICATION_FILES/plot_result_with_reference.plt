#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_27_cable_test_current_ref.jpg"

set title "Cable Current"
set xlabel "Time (s)"
set ylabel "A"
plot "PROBLEM_SPECIFICATION_FILES/Reference_current.reference" u 1:3 title "Cable current: GGI_TLM Reference" w p ,\
     "cable_test.cable_current.tout" u 1:3 title "Cable current: GGI_TLM" w l 
     
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_27_cable_test_field_ref.jpg"
set title "Electric Field"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot "PROBLEM_SPECIFICATION_FILES/cable_test_fieldsolve_TLM.field.tout" u 1:3 title "E field: Fieldsolve Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/Reference_field.reference" u 1:3 title "E field: GGI_TLM Reference" w p ,\
     "cable_test.field.tout" u 1:3 title "E field: GGI_TLM" w l
#PAUSE pause -1
