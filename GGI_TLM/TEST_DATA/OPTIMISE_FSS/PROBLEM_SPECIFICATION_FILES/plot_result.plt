
#JPG set term jpeg
#OUTPUT_TO_FILE set output 'S11_dB.jpg'
set autoscale x
set autoscale y
set title 'S11 '
set xlabel 'Frequency (Hz)'
set ylabel 'S11 (dB)'
plot 'S21.fout' u 1:5 title '  |S11|' w l 

#PAUSE pause -1 
  
#OUTPUT_TO_FILE set output 'S21_dB.jpg'
set autoscale x
set autoscale y
set yrange [-40:5]
set title 'S21 '
set xlabel 'Frequency (Hz)'
set ylabel 'S21 (dB)'
plot 'S21.fout' u 1:9 title '  |S21|' w l 

#PAUSE pause -1 
