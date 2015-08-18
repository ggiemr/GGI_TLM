#JPG set term jpeg


#OUTPUT_TO_FILE set output "Test_case_10_FAR_FIELD.jpg"
set autoscale x
set autoscale y
set xlabel "theta (degrees)"
set ylabel "Far field"
plot "plate.far_field.fout.1" u 1:3 title "E theta: GGI_TLM" w l,\
     "plate.far_field.fout.1" u 1:4 title "E phi: GGI_TLM" w l
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_10_RCS_PEC_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "E field"
plot "plate.field.tout" u 1:3 title "GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_10_RCS_PEC_far_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "Normalised far field"
plot "plate.rcs.tout" u 2:3 title "GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_10_RCS_PEC_RCS.jpg"
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "RCS (dBsm)"
plot "plate.rcs.fout" u 1:5 title "GGI_TLM" w l
#PAUSE pause -1
