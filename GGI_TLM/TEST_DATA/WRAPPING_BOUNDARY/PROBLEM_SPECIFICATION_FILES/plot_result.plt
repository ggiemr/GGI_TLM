#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_56_WRAPPING_SURFACE_Ex_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot \
"Ex_centre.tout" u 1:3 w l 
#PAUSE pause -1

