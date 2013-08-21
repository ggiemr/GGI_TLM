#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_25_dielectric_filter_S_params.jpg"
set autoscale x
set autoscale y
set xlabel "Frequency(Hz)"
set ylabel "dB"
plot "filter_S11.fout" u 1:7 title "S11" w l,\
     "filter_S21.fout" u 1:7 title "S21" w l
     
#PAUSE pause -1

