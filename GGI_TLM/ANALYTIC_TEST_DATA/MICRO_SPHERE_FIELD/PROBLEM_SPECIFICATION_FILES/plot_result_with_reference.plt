#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_3_MICRO_SPHERE_FIELD_field_ref.jpg"
set xlabel "Time (s)"
set ylabel "Electric Field (V/m)"
plot "micro_sphere.field" u 1:2 title "E field: GGI_RCS_sphere" w l
     "PROBLEM_SPECIFICATION_FILES/micro_sphere.field_ref" u 1:2 title "E field: GGI_RCS_sphere reference" w p
     
#PAUSE pause -1
