#JPG set term jpeg

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V/m"

#OUTPUT_TO_FILE set output "Test_case_4_TEM_THIN_LAYER_field_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/Ey.tout_ref" u 1:3 title "Ey: GGI_TLM reference" w p,\
     "Ey.tout" u 1:3 title "Ey: GGI_TLM" w l
#PAUSE pause -1


set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_4_TEM_THIN_LAYER_beta_ref.jpg"

set xlabel "Frequency (Hz)"
set ylabel "Propagetion constant (m^-1)"
plot "PROBLEM_SPECIFICATION_FILES/TEM_material.fout_ref" u 1:2 title "Re{beta}, GGI_TLM reference" w p ,\
     "PROBLEM_SPECIFICATION_FILES/TEM_material.fout_ref" u 1:3 title "Im{beta}, GGI_TLM reference" w p ,\
     "TEM_material.fout" u 1:2 title "Re{beta}, GGI_TLM" w l ,\
     "TEM_material.fout" u 1:3 title "Im{beta}, GGI_TLM" w l 
#PAUSE pause -1

set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_4_TEM_THIN_LAYER_Z0_ref.jpg"

set xlabel "Frequency (Hz)"
set ylabel "Wave impedance (ohms)"
plot "PROBLEM_SPECIFICATION_FILES/TEM_material.fout_ref" u 1:4 title "Re{Z}, GGI_TLM reference" w p ,\
     "PROBLEM_SPECIFICATION_FILES/TEM_material.fout_ref" u 1:5 title "Im{Z}, GGI_TLM reference" w p ,\
     "TEM_material.fout" u 1:4 title "Re{Z}, GGI_TLM" w l ,\
     "TEM_material.fout" u 1:5 title "Im{Z}, GGI_TLM" w l 
#PAUSE pause -1

set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_4_TEM_THIN_LAYER_R_ref.jpg"

set xlabel "Frequency (Hz)"
set ylabel "Reflection coefficient "
plot "PROBLEM_SPECIFICATION_FILES/TEM_material.fout_ref" u 1:6 title "Re{R}, GGI_TLM reference" w p ,\
     "PROBLEM_SPECIFICATION_FILES/TEM_material.fout_ref" u 1:7 title "Im{R}, GGI_TLM reference" w p ,\
     "TEM_material.fout" u 1:6 title "Re{R}, GGI_TLM" w l ,\
     "TEM_material.fout" u 1:7 title "Im{R}, GGI_TLM" w l 
#PAUSE pause -1
