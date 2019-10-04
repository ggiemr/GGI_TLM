#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_23_waveguide_pml_S11_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "S11 (dB)"
plot "PROBLEM_SPECIFICATION_FILES/waveguide_pml.fout_ref" u 1:9 title "Waveguide iris S11: GGI_TLM reference" w p,\
"waveguide_pml.fout" u 1:9 title "Waveguide iris S11: GGI_TLM" w l
#PAUSE pause -1


