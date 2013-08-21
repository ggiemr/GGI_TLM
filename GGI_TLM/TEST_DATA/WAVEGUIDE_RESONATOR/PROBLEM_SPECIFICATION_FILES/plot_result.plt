#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_20_WAVEGUIDE_RESONATOR_field.jpg"
set xlabel "Frequency (Hz)"
set ylabel "V/m"
plot "waveguide.field.fout" u 1:5 title "Cavity E field: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_20_WAVEGUIDE_RESONATOR_mode_field.jpg"
set xlabel "x cell"
set ylabel "y cell"
splot "WR90.mode" u 1:2:6 title "Mode E field: Fieldsolve" w p
#PAUSE pause -1

