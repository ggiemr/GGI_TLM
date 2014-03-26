
#JPG set term jpeg
#OUTPUT_TO_FILE set output 'S11_dB.jpg'
set autoscale x
set autoscale y
set title 'S11 '
set xlabel 'Frequency (Hz)'
set ylabel 'S11 (dB)'
plot \
 'RESULTS/S21.fout_1' u 1:5 title '  |S11|: Run number 1' w l \
,'RESULTS/S21.fout_2' u 1:5 title '  |S11|: Run number 2' w l \
,'RESULTS/S21.fout_3' u 1:5 title '  |S11|: Run number 3' w l \
,'RESULTS/S21.fout_4' u 1:5 title '  |S11|: Run number 4' w l \

  #PAUSE pause -1 
  
#OUTPUT_TO_FILE set output 'S21_dB.jpg'
set autoscale x
set autoscale y
set title 'S21 '
set xlabel 'Frequency (Hz)'
set ylabel 'S21 (dB)'
plot \
 'RESULTS/S21.fout_1' u 1:9 title '  |S21|: Run number 1' w l \
,'RESULTS/S21.fout_2' u 1:9 title '  |S21|: Run number 2' w l \
,'RESULTS/S21.fout_3' u 1:9 title '  |S21|: Run number 3' w l \
,'RESULTS/S21.fout_4' u 1:9 title '  |S21|: Run number 4' w l \

  #PAUSE pause -1 
