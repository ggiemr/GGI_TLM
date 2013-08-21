#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_31_dipole_current_ref.jpg"
set xlabel "Time (s)"
set ylabel "A"
plot "PROBLEM_SPECIFICATION_FILES/dipole_original_TLM.cable_current.tout" u 1:3 title "Dipole current: Fieldsolve Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/dipole.cable_current.tout_ref" u 1:3 title "Dipole current: GGI_TLM Reference" w p ,\
     "dipole.cable_current.tout" u 1:3 title "Dipole current: GGI_TLM" w l 
     
#PAUSE pause -1

set autoscale x
set autoscale y


#OUTPUT_TO_FILE set output "Test_case_31_dipole_field_ref.jpg"
set xlabel "Time (s)"
set ylabel "V/m"
plot "PROBLEM_SPECIFICATION_FILES/dipole_original_TLM.field.tout" u 1:3 title "Dipole field: Fieldsolve Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/dipole.field.tout_ref" u 1:3 title "Dipole field: GGI_TLM Reference" w p ,\
     "dipole.field.tout" u 1:3 title "Dipole field: GGI_TLM" w l 
     
#PAUSE pause -1
