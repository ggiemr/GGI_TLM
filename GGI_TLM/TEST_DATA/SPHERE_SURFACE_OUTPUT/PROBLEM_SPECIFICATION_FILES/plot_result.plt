#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_38_SPHERE_SURFACE_OUTPUT_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "E field"
plot "sphere.field.tout" u 1:3 title "GGI_TLM" w l
#PAUSE pause -1
