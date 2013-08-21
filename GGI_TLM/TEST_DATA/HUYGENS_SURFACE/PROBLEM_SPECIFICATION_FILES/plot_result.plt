#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_55_HUYGENS_SURFACE_Ex_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot \
"Ex_centre.tout" u 1:3 w l ,\
"Ex_ymin.tout" u 1:3 w l ,\
"Ex_ymax.tout" u 1:3 w l ,\
"Ex_zmin.tout" u 1:3 w l ,\
"Ex_zmax.tout" u 1:3 w l 
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_55_HUYGENS_SURFACE_Ey_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot \
"Ey_centre.tout" u 1:3 w l ,\
"Ey_xmin.tout" u 1:3 w l ,\
"Ey_xmax.tout" u 1:3 w l ,\
"Ey_zmin.tout" u 1:3 w l ,\
"Ey_zmax.tout" u 1:3 w l 
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_55_HUYGENS_SURFACE_Ez_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot \
"Ez_centre.tout" u 1:3 w l ,\
"Ez_xmin.tout" u 1:3 w l ,\
"Ez_xmax.tout" u 1:3 w l ,\
"Ez_ymin.tout" u 1:3 w l ,\
"Ez_ymax.tout" u 1:3 w l 
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_55_HUYGENS_SURFACE_Hx_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "A/m"
plot \
"Hx_centre.tout" u 1:3 w l ,\
"Hx_ymin.tout" u 1:3 w l ,\
"Hx_ymax.tout" u 1:3 w l ,\
"Hx_zmin.tout" u 1:3 w l ,\
"Hx_zmax.tout" u 1:3 w l 
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_55_HUYGENS_SURFACE_Hy_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "A/m"
plot \
"Hy_centre.tout" u 1:3 w l ,\
"Hy_xmin.tout" u 1:3 w l ,\
"Hy_xmax.tout" u 1:3 w l ,\
"Hy_zmin.tout" u 1:3 w l ,\
"Hy_zmax.tout" u 1:3 w l 
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_55_HUYGENS_SURFACE_Hz_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "A/m"
plot \
"Hz_centre.tout" u 1:3 w l ,\
"Hz_xmin.tout" u 1:3 w l ,\
"Hz_xmax.tout" u 1:3 w l ,\
"Hz_ymin.tout" u 1:3 w l ,\
"Hz_ymax.tout" u 1:3 w l 
#PAUSE pause -1

