#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_23_WAVEGUIDE_IRIS_field.jpg"
set xlabel "Time (s)"
set ylabel "Field"
plot "waveguide_iris.field.tout" u 1:3 title "Waveguide field: GGI_TLM" w p
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_23_WAVEGUIDE_IRIS_S11.jpg"
set xlabel "Frequency (Hz)"
set ylabel "S11 (dB)"
plot "waveguide_iris.fout" u 1:5 title "S11: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_23_WAVEGUIDE_IRIS_S21.jpg"
set xlabel "Frequency (Hz)"
set ylabel "S21 (dB)"
plot "waveguide_iris.fout" u 1:9 title "S21: GGI_TLM" w l
#PAUSE pause -1

