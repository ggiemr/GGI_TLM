#JPG set term jpeg

#OUTPUT_TO_FILE set output "parallel_test.jpg"

set xlabel "Frequency (Hz)"
set ylabel "V/m"
plot "cavity.fout" u 1:5 title "Cavity E field: GGI_TLM" w l

#PAUSE pause -1
