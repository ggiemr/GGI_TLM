#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_11_RCS_MATERIAL_field_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "E field"
plot "sphere.field.tout" u 1:3 title "GGI_TLM" w l,\
     "PROBLEM_SPECIFICATION_FILES/sphere.field.tout_ref" u 1:3 title "Fieldsolve reference" w p ,\
     "PROBLEM_SPECIFICATION_FILES/sphere_original_TLM.field.tout" u 1:3 title "Fieldsolve reference" w p
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_11_RCS_MATERIAL_far_field_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "Normalised far field"
plot "sphere.rcs.tout" u 2:3 title "GGI_TLM" w l,\
     "PROBLEM_SPECIFICATION_FILES/sphere.rcs.tout_ref" u 2:3 title "GGI_TLM reference" w p ,\
     "PROBLEM_SPECIFICATION_FILES/sphere_original_TLM.mfar" u 2:3 title "Fieldsolve reference" w p
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_11_RCS_MATERIAL_RCS_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "RCS (dBsm)"
plot "sphere.rcs.fout" u 1:5 title "GGI_TLM" w l,\
     "PROBLEM_SPECIFICATION_FILES/sphere.rcs.fout_ref" u 1:5 title "GGI_TLM reference" w p ,\
     "PROBLEM_SPECIFICATION_FILES/sphere_original_TLM.mrcs" u 1:5 title "Fieldsolve reference" w p
#PAUSE pause -1
