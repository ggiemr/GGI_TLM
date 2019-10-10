#JPG set term jpeg

set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_5_TEM_MATERIAL_field_ref.jpg"

set xlabel "Time (s)"
set ylabel "V/m"
plot "PROBLEM_SPECIFICATION_FILES/Ey.tout_ref" u 1:3 title "Ey: GGI_TLM reference" w p,\
     "Ey.tout" u 1:3 title "Ey: GGI_TLM" w l
#PAUSE pause -1


set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_5_TEM_MATERIAL_R_ref.jpg"

set xlabel "Frequency (Hz)"
set ylabel "Reflection coefficient (dB)"
plot "PROBLEM_SPECIFICATION_FILES/TEM_material.fout_ref" u 1:9 title "S11, GGI_TLM reference" w p ,\
     "TEM_material.fout" u 1:9 title "S11, GGI_TLM" w l 
#PAUSE pause -1
