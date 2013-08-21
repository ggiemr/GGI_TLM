#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_34_MATERIAL_TEST_field_ref.jpg"
set xlabel "Time (s)"
set ylabel "V/m"
plot "PROBLEM_SPECIFICATION_FILES/mat_test.field.tout_ref" u 1:3 title "E field: GGI_TLM Reference" w p ,\
     "PROBLEM_SPECIFICATION_FILES/mat_test_original_TLM.field.tout" u 1:3 title "E field: Fieldsolve Reference" w p ,\
     "mat_test.field.tout" u 1:3 title "E field: GGI_TLM" w l 
     
#PAUSE pause -1
