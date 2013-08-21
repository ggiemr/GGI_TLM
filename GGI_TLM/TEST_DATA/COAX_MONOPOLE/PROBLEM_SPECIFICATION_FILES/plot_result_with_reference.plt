#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_30_monopole_current_ref.jpg"
set xlabel "Time (s)"
set ylabel "A"
plot "PROBLEM_SPECIFICATION_FILES/coax_monopole_original_TLM.cable_current.tout" u 1:3 title "Coax monopole current: Fieldsolve Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/coax_monopole.cable_current.tout_ref" u 1:3 title "Coax monopole current: GGI_TLM Reference" w p ,\
     "coax_monopole.cable_current.tout" u 1:3 title "Coax monopole current: GGI_TLM" w l    
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_30_monopole_field_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot "PROBLEM_SPECIFICATION_FILES/coax_monopole_original_TLM.field.tout" u 1:3 title "Coax monopole field: Fieldsolve Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/coax_monopole.field.tout_ref" u 1:3 title "Coax monopole field: GGI_TLM Reference" w p ,\
     "coax_monopole.field.tout" u 1:3 title "Coax monopole field: GGI_TLM" w l 
#PAUSE pause -1
