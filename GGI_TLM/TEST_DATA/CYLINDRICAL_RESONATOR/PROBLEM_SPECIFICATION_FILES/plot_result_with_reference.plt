#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_1_cavity_free_space_ref.jpg"

set xlabel "Frequency (Hz)"
set ylabel "V/m"

plot "cylinder.fout" u 1:5 title "Cavity E field: GGI_TLM" w l,\
"PROBLEM_SPECIFICATION_FILES/cylinder.fout_ref" u 1:5 title "Cavity E field: GGI_TLM reference" w p

#PAUSE pause -1

