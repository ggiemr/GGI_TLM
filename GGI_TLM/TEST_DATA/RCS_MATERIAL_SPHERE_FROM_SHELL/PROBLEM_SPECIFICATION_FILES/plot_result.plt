#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_11_RCS_MATERIAL_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "E field"
plot "sphere.field.tout" u 1:3 title "GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_11_RCS_MATERIAL_far_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "Normalised far field"
plot "sphere.rcs.tout" u 2:3 title "GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_11_RCS_MATERIAL_RCS.jpg"
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "RCS (dBsm)"
plot "sphere.rcs.fout" u 1:5 title "GGI_TLM" w l
#PAUSE pause -1
