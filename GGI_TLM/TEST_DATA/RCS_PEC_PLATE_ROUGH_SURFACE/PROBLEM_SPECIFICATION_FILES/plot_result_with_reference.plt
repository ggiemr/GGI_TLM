#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_10_FAR_FIELD_ref.jpg"
set autoscale x
set autoscale y
set xlabel "theta (degrees)"
set ylabel "Far field"
plot "PROBLEM_SPECIFICATION_FILES/plate.far_field.fout.1_ref" u 1:3 title "E theta: GGI_TLM" w p,\
     "PROBLEM_SPECIFICATION_FILES/plate.far_field.fout.1_ref" u 1:4 title "E phi: GGI_TLM" w p,\
     "plate.far_field.fout.1" u 1:3 title "E theta: GGI_TLM" w l,\
     "plate.far_field.fout.1" u 1:4 title "E phi: GGI_TLM" w l
     
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_10_RCS_PEC_field_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "E field"
plot      "PROBLEM_SPECIFICATION_FILES/plate.field.tout_ref" u 1:3 title "GGI_TLM reference" w p ,\
"plate.field.tout" u 1:3 title "GGI_TLM" w l

#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_10_RCS_PEC_far_field_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "Normalised far field"
plot      "PROBLEM_SPECIFICATION_FILES/plate.rcs.tout_ref" u 2:3 title "GGI_TLM reference" w p ,\
"plate.rcs.tout" u 2:3 title "GGI_TLM" w l

#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_10_RCS_PEC_RCS_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "RCS (dBsm)"
plot      "PROBLEM_SPECIFICATION_FILES/plate.rcs.fout_ref" u 1:5 title "GGI_TLM reference" w p ,\
"plate.rcs.fout" u 1:5 title "GGI_TLM" w l

#PAUSE pause -1
