#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_20_WAVEGUIDE_RESONATOR_field_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "V/m"
plot "PROBLEM_SPECIFICATION_FILES/waveguide.field.fout_ref" u 1:5 title "Cavity E field: GGI_TLM reference" w p ,\
"waveguide.field.fout" u 1:5 title "Cavity E field: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_20_WAVEGUIDE_RESONATOR_mode_field_ref.jpg"
set xlabel "x cell"
set ylabel "y cell"
splot "PROBLEM_SPECIFICATION_FILES/WR90.mode_ref" u 1:2:6 title "Mode E field:  GGI_TLM reference" w p ,\
"WR90.mode" u 1:2:6 title "Mode E field: GGI_TLM" w l
#PAUSE pause -1
