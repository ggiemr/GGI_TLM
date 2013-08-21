#JPG set term jpeg

set title "Monopole current"

set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_9_MONOPOLE_current_ref.jpg"

set xlabel "Time (s)"
set ylabel "A"

plot "PROBLEM_SPECIFICATION_FILES/monopole_original_TLM.cable_current.tout" u 1:3 title "Fieldsolve reference" w p,\
     "PROBLEM_SPECIFICATION_FILES/Reference_monopole_current.reference" u 1:3 title "GGI_TLM reference" w p,\
     "monopole.cable_current.tout" u 1:3 title "GGI_TLM" w l
#PAUSE pause -1


set title "Monopole field"

set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_9_MONOPOLE_field_ref.jpg"

set xlabel "Time (s)"
set ylabel "V/m"
plot "PROBLEM_SPECIFICATION_FILES/monopole_original_TLM.field.tout" u 1:3 title "Fieldsolve reference" w p,\
     "PROBLEM_SPECIFICATION_FILES/Reference_monopole_field.reference" u 1:3 title "GGI_TLM reference" w p,\
     "monopole.field.tout" u 1:3 title "GGI_TLM" w l
pause -1
