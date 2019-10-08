#JPG set term jpeg

set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_5_TEM_MATERIAL_field.jpg"

set xlabel "Time (s)"
set ylabel "V/m"
plot "Ey.tout" u 1:3 title "Ey: GGI_TLM" w l 
#PAUSE pause -1


set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_5_TEM_MATERIAL_S11.jpg"

set xlabel "Frequency (Hz)"
set ylabel "S11(dB)"
plot "TEM_material.fout" u 1:9 title "S11(dB), GGI_TLM" w l 
#PAUSE pause -1
