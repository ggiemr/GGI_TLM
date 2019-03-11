#JPG set term jpeg

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V"

#OUTPUT_TO_FILE set output "Source_voltage.jpg"

plot "two_single_wires.excitation.tout" u 1:3 title "source voltage: GGI_TLM" w l 
     
#PAUSE pause -1

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "A"

#OUTPUT_TO_FILE set output "Source_currrent.jpg"

plot "two_single_wires.cable_current.tout" u 1:3 title "source current: GGI_TLM" w l 
     
#PAUSE pause -1

#JPG set term jpeg

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V"

#OUTPUT_TO_FILE set output "Source_voltage.jpg"

plot "source_voltage.tout" u 1:3 title "source voltage: GGI_TLM" w l 
     
#PAUSE pause -1
