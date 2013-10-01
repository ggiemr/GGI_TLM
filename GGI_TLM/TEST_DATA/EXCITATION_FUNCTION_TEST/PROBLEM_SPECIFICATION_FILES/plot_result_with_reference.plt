#JPG set term jpeg

#OUTPUT_TO_FILE set output "EXCITATION_TEST_Excitation_function_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot \
"excitation_test.excitation.tout" u 1:3 w l,\
"PROBLEM_SPECIFICATION_FILES/gaussian_from_file.dat" u 1:2 w p
#PAUSE pause -1


#OUTPUT_TO_FILE set output "EXCITATION_TEST_Ex_field_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot \
"PROBLEM_SPECIFICATION_FILES/excitation.tout_ref" u 1:3 w p,\
"excitation_test.field.tout" u 1:3 w l
#PAUSE pause -1

