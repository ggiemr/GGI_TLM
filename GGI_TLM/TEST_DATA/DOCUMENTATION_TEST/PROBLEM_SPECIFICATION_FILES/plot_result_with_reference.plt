#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_DOCUMENTATION_TEST_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "RCS (dBsm)"
plot "PROBLEM_SPECIFICATION_FILES/geometry.rcs.fout_ref" u 1:5 title "GGI_TLM reference" w p ,\
     "geometry.rcs.fout" u 1:5 title "GGI_TLM" w l 
#PAUSE pause -1
