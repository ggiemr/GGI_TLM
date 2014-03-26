#JPG set term jpeg

#OUTPUT_TO_FILE set output "S21_S11_dB_ref.jpg"
set autoscale x
set autoscale y
set yrange [-40:5]
set xlabel "Frequency (Hz)"
set ylabel "S11/S21 (dB)"
plot "PROBLEM_SPECIFICATION_FILES/S21.fout_ref" u 1:5 title "  |S11|: GGI_TLM reference" w p,\
     "PROBLEM_SPECIFICATION_FILES/S21.fout_ref" u 1:9 title "  |S21|: GGI_TLM reference" w p,\
     "S21.fout" u 1:5 title "  |S11|: GGI_TLM" w l,\
     "S21.fout" u 1:9 title "  |S21|: GGI_TLM" w l
#PAUSE pause -1
