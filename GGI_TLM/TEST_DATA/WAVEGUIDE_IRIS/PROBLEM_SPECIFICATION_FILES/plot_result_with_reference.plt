#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_23_WAVEGUIDE_IRIS_S11_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "S11 (dB)"
plot "PROBLEM_SPECIFICATION_FILES/waveguide_iris.fout_ref" u 1:5 title "Waveguide iris S11: GGI_TLM reference" w p,\
"waveguide_iris.fout" u 1:5 title "Waveguide iris S11: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_23_WAVEGUIDE_IRIS_S21_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "S21 (dB)"
plot "PROBLEM_SPECIFICATION_FILES/waveguide_iris.fout_ref" u 1:9 title "Waveguide iris S21: GGI_TLM reference" w p,\
"waveguide_iris.fout" u 1:9 title "Waveguide iris S21: GGI_TLM" w l
#PAUSE pause -1

