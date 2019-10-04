#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_23_waveguide_pml_field.jpg"
set xlabel "Time (s)"
set ylabel "Field"
plot "waveguide_pml.field.tout" u 1:3 title "Waveguide field: GGI_TLM" w p
  pause -1

#OUTPUT_TO_FILE set output "Test_case_23_waveguide_pml_S11.jpg"
set xlabel "Frequency (Hz)"
set ylabel "S11 (dB)"
plot "waveguide_pml.fout" u 1:9 title "S11: GGI_TLM" w l
  pause -1

