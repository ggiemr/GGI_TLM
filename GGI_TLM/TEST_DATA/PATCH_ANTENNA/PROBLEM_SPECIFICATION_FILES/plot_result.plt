
#JPG set term jpeg

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "A"

#OUTPUT_TO_FILE set output "Test_case_8_patch_current.jpg"

plot "patch.cable_current.tout" u 1:3 title "patch current: GGI_TLM" w l 
     
#PAUSE pause -1

#___________________________________________________________________#

set xlabel "Time (s)"
set ylabel "V/m"

set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_8_patch_field.jpg"

plot "patch.field.tout" u 1:3 title "patch field: GGI_TLM" w l 
     
#PAUSE pause -1


#___________________________________________________________________#


set autoscale x
set autoscale y

set title "patch reflection coefficient"
set xlabel "Frequency(Hz)"
set ylabel "dB"

#OUTPUT_TO_FILE set output "Test_case_8_patch_S11.jpg"

plot "patch.S11.fout" u 1:7 title "|S11|: GGI_TLM " w l

#PAUSE pause -1


