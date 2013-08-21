#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_55_HUYGENS_SURFACE_Ex_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot "PROBLEM_SPECIFICATION_FILES/Ex_centre.tout_ref" u 1:3 title "GGI_TLM: reference" w p, \
"Ex_centre.tout" u 1:3 w l 
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_55_HUYGENS_SURFACE_Ey_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot "PROBLEM_SPECIFICATION_FILES/Ey_centre.tout_ref" u 1:3 title "GGI_TLM: reference" w p, \
"Ey_centre.tout" u 1:3 w l 
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_55_HUYGENS_SURFACE_Ez_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot "PROBLEM_SPECIFICATION_FILES/Ez_centre.tout_ref" u 1:3 title "GGI_TLM: reference" w p, \
"Ez_centre.tout" u 1:3 w l 
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_55_HUYGENS_SURFACE_Hx_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "A/m"
plot "PROBLEM_SPECIFICATION_FILES/Hx_centre.tout_ref" u 1:3 title "GGI_TLM: reference" w p, \
"Hx_centre.tout" u 1:3 w l 
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_55_HUYGENS_SURFACE_Hy_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "A/m"
plot "PROBLEM_SPECIFICATION_FILES/Hy_centre.tout_ref" u 1:3 title "GGI_TLM: reference" w p, \
"Hy_centre.tout" u 1:3 w l 
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_55_HUYGENS_SURFACE_Hz_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "A/m"
plot "PROBLEM_SPECIFICATION_FILES/Hz_centre.tout_ref" u 1:3 title "GGI_TLM: reference" w p, \
"Hz_centre.tout" u 1:3 w l 
#PAUSE pause -1

