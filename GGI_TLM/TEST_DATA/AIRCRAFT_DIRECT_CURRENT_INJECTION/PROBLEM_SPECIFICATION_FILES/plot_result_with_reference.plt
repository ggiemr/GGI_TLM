#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_45_AIRCRAFT_SURFACE_CURRENT_field_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "E field"
plot "aircraft.field.tout" u 1:3 title "GGI_TLM" w l,\
     "PROBLEM_SPECIFICATION_FILES/aircraft.field.tout_ref" u 1:3 title "GGI_TLM reference" w p 
#PAUSE pause -1

