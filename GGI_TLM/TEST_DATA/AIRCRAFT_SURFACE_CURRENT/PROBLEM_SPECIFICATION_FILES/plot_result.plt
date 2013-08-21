#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_45_AIRCRAFT_SURFACE_CURRENT_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "E field"
plot "aircraft.field.tout" u 1:3 title "GGI_TLM" w l
#PAUSE pause -1
