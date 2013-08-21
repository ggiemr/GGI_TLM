#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_23_WAVEGUIDE_IRIS_S11_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "S11 (dB)"
plot "waveguide_iris.fout" u 1:5 title "Waveguide iris S11: GGI_TLM" w l,\
"PROBLEM_SPECIFICATION_FILES/waveguide_iris.fout_ref" u 1:5 title "Waveguide iris S11: GGI_TLM reference" w p,\
"PROBLEM_SPECIFICATION_FILES/waveguide_iris1_original_TLM.fout.S11" u 1:5 title "Waveguide iris S11: Fieldsolve" w p
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_23_WAVEGUIDE_IRIS_S21_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "S21 (dB)"
plot "waveguide_iris.fout" u 1:9 title "Waveguide iris S21: GGI_TLM" w l,\
"PROBLEM_SPECIFICATION_FILES/waveguide_iris.fout_ref" u 1:9 title "Waveguide iris S21: GGI_TLM reference" w p,\
"PROBLEM_SPECIFICATION_FILES/waveguide_iris1_original_TLM.fout.S21" u 1:5 title "Waveguide iris S21: Fieldsolve" w p
#PAUSE pause -1

