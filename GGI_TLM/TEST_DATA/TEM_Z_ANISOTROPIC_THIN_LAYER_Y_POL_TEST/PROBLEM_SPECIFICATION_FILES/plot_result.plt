#JPG set term jpeg

set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_4_TEM_THIN_LAYER_field.jpg"

set xlabel "Time (s)"
set ylabel "V/m"
plot "Ey.tout" u 1:3 title "Ey: GGI_TLM" w l 
#PAUSE pause -1


set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_4_TEM_THIN_LAYER_beta.jpg"

set xlabel "Frequency (Hz)"
set ylabel "Propagetion constant (m^-1)"
plot "TEM_material.fout" u 1:2 title "Re{beta}, GGI_TLM" w l ,\
     "TEM_material.fout" u 1:3 title "Im{beta}, GGI_TLM" w l 
#PAUSE pause -1

set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_4_TEM_THIN_LAYER_Z0.jpg"

set xlabel "Frequency (Hz)"
set ylabel "Wave impedance (ohms)"
plot "TEM_material.fout" u 1:4 title "Re{Z}, GGI_TLM" w l ,\
     "TEM_material.fout" u 1:5 title "Im{Z}, GGI_TLM" w l 
#PAUSE pause -1

set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "Reflection coefficient "

#OUTPUT_TO_FILE set output "Test_case_4_TEM_THIN_LAYER_R.jpg"

plot "TEM_material.fout" u 1:6 title "Re{R}, GGI_TLM" w l ,\
     "TEM_material.fout" u 1:7 title "Im{R}, GGI_TLM" w l 
#PAUSE pause -1
