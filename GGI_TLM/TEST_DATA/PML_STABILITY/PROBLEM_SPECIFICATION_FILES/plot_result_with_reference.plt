#JPG set term jpeg

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V/m"

#OUTPUT_TO_FILE set output "Test_case_8_dipole_field_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/cavity.field.tout_ref" u 1:3 title "E field: GGI_TLM Reference" w p ,\
     "cavity.field.tout" u 1:3 title "E field: GGI_TLM" w l 
     
#PAUSE pause -1

