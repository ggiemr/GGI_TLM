
#JPG set term jpeg

set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_6_MATERIAL_STABILITY_TEST_field_ref.jpg"

set xlabel "Time (s)"
set ylabel "V/m"
plot "PROBLEM_SPECIFICATION_FILES/cavity.field.tout_ref" u 1:3 title "Ey: GGI_TLM reference" w p,\
     "cavity.field.tout" u 1:3 title "Ey: GGI_TLM" w l
#PAUSE pause -1

