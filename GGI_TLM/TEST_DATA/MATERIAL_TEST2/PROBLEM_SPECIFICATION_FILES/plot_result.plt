#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_34_MATERIAL_TEST_field.jpg"
set xlabel "Time (s)"
set ylabel "V/m"
plot "mat_test.field.tout" u 1:3 title "E field: GGI_TLM" w l
#PAUSE pause -1
