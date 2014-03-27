#JPG set term jpeg

#OUTPUT_TO_FILE set output "parallel_test_ref.jpg"

set xlabel "Frequency (Hz)"
set ylabel "V/m"

plot "PROBLEM_SPECIFICATION_FILES/cavity.fout_ref" u 1:5 title "Cavity E field: GGI_TLM reference" w p,\
"cavity.fout" u 1:5 title "Cavity E field: GGI_TLM" w l

#PAUSE pause -1

