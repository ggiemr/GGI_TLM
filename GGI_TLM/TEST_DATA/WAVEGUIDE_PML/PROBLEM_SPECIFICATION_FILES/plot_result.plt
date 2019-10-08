#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_23_waveguide_pml_field.jpg"
set xlabel "Time (s)"
set ylabel "Field"
plot "Ex1.tout" u 1:3 title "Waveguide E field: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_23_waveguide_pml_S11.jpg"
set xlabel "Frequency (Hz)"
set ylabel "S11 (dB)"
plot "waveguide_pml.fout" u 1:9 title "S11: GGI_TLM" w l
#PAUSE pause -1

