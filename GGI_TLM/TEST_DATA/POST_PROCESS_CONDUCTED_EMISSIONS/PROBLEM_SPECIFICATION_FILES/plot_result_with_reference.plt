#JPG set term jpeg

set autoscale x
set autoscale y
set logscale x
set grid

set xlabel "Frequency (Hz)"
set ylabel "dbuV"

#OUTPUT_TO_FILE set output "Conducted_voltage_ref.jpg"

plot "ce336_600w.favg_ref" u 1:4 title "Conducted voltage: reference" w l ,\
     "ce336_600w.favg" u 1:4 title "Conducted voltage with detector bandwidth" w l ,\
     "sub_segment_freq_average.dat" u 1:6 title "Conducted voltage" w l 
pause -1

#PAUSE pause -1
