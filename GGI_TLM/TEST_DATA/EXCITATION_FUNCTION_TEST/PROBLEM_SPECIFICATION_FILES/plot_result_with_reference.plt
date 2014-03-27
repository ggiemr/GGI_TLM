#JPG set term jpeg

#OUTPUT_TO_FILE set output "EXCITATION_TEST_Ex_field_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot \
"PROBLEM_SPECIFICATION_FILES/excitation.tout_ref" u 1:3 title "GGI_TLM reference" w p,\
"excitation_test.field.tout" u 1:3 title "GGI_TLM" w l
#PAUSE pause -1

