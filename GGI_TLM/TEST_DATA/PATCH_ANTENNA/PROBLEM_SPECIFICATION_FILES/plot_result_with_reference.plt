#JPG set term jpeg

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "A"

#OUTPUT_TO_FILE set output "Test_case_8_patch_current_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/patch.cable_current.tout_ref" u 1:3 title "patch current: GGI_TLM Reference" w p ,\
     "patch.cable_current.tout" u 1:3 title "patch current: GGI_TLM" w l 
     
#PAUSE pause -1

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V/m"

#OUTPUT_TO_FILE set output "Test_case_8_patch_field_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/patch.field.tout_ref" u 1:3 title "patch field: GGI_TLM Reference" w p ,\
     "patch.field.tout" u 1:3 title "patch field: GGI_TLM" w l 
     
#PAUSE pause -1

set autoscale x
set autoscale y

set title "patch reflection coefficient"
set xlabel "Frequency(Hz)"
set ylabel "dB"

#OUTPUT_TO_FILE set output "Test_case_8_patch_S11.jpg"

plot "PROBLEM_SPECIFICATION_FILES/patch.S11.fout_ref" u 1:7 title "|S11|: GGI_TLM reference" w p ,\
     "patch.S11.fout" u 1:7 title "|S11|: GGI_TLM " w l

#PAUSE pause -1

