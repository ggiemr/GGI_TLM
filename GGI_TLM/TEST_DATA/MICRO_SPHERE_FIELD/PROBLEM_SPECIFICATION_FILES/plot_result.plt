#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_47_MICRO_SPHERE_FIELD_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "E field"
plot "sphere.field.tout" u 1:3 title "GGI_TLM" w l,\
     "PROBLEM_SPECIFICATION_FILES/micro_sphere.field_ref" u 1:2 title "Analytic reference" w p 
#PAUSE pause -1
