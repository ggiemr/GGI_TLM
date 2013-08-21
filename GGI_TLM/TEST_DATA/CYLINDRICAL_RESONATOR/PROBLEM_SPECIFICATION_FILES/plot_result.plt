#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_1_cavity_free_space.jpg"

set xlabel "Frequency (Hz)"
set ylabel "V/m"
plot "cylinder.fout" u 1:5 title "Cavity E field: GGI_TLM" w l

#PAUSE pause -1
