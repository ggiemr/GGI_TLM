#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_23_waveguide_pml_Ex_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "V/m"
plot "PROBLEM_SPECIFICATION_FILES/Ex1.tout_ref" u 1:3 title "Waveguide PML Ex: GGI_TLM reference" w p,\
"Ex1.tout" u 1:3 title "Waveguide PML Ex: GGI_TLM" w l
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_23_waveguide_pml_S11_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "S11 (dB)"
plot "PROBLEM_SPECIFICATION_FILES/waveguide_pml.fout_ref" u 1:9 title "Waveguide PML S11: GGI_TLM reference" w p,\
"waveguide_pml.fout" u 1:9 title "Waveguide PML S11: GGI_TLM" w l
#PAUSE pause -1


