#JPG set term jpeg

set autoscale x
set autoscale y
set logscale x

set xlabel "Frequency (Hz)"
set ylabel "dbuV"

#OUTPUT_TO_FILE set output "Conducted_voltage.jpg"

plot "ce336_600w.favg" u 1:4 title "Conducted voltage" w l 
     
#PAUSE pause -1


