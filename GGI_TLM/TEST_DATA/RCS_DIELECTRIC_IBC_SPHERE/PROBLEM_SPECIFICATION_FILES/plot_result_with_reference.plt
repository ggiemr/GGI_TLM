#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_49_RCS_DIELECTRIC_IBC_far_field_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "Normalised far field"
plot "sphere.rcs.tout" u 2:3 title "GGI_TLM" w l,\
     "PROBLEM_SPECIFICATION_FILES/sphere.rcs.tout_ref" u 2:3 title "GGI_TLM reference" w p 
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_49_RCS_DIELECTRIC_IBC_rcs_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "RCS (dBsm)"
plot "sphere.rcs.fout" u 1:5 title "GGI_TLM" w l,\
     "PROBLEM_SPECIFICATION_FILES/sphere.rcs.fout_ref" u 1:5 title "GGI_TLM reference" w p ,\
     "PROBLEM_SPECIFICATION_FILES/dielectric_ibc.rcs_ref" u 1:4 title "Analytic result" w p
#PAUSE pause -1
