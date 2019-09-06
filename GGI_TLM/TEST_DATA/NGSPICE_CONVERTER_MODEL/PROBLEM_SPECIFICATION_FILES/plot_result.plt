#JPG set term jpeg

set autoscale x

set xlabel "Time (s)"
#JPG set term jpeg

set ylabel "V"
set autoscale x 
set autoscale y
set yrange[-20:50]

#OUTPUT_TO_FILE set output "Converter_voltages_1.jpg"

plot "v1.tout" u 1:3 title "Source voltage" w l ,\
     "v2.tout" u 1:3 title "Switch voltage" w l ,\
     "v3.tout" u 1:3 title "Diode voltage" w l ,\
     "v4.tout" u 1:3 title "Inductor voltage" w l  ,\
     "v6.tout" u 1:3 title "Load voltage" w l 
     
#PAUSE pause -1

set xlabel "Time (s)"
set ylabel "V"
set autoscale x 
set autoscale y
set yrange[-50:50]

#OUTPUT_TO_FILE set output "Converter_voltages_2.jpg"

plot "v5.tout" u 1:3 title "Capacitor voltage" w l ,\
     "v6.tout" u 1:3 title "Load voltage" w l 
     
#PAUSE pause -1

quit

