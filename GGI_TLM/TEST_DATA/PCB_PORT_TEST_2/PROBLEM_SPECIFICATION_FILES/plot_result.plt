#JPG set term jpeg

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V"

#OUTPUT_TO_FILE set output "Source_voltage.jpg"

plot "source_voltage.tout" u 1:3 title "source voltage: GGI_TLM" w l ,\
     "load_voltage.tout" u 1:3 title "load voltage: GGI_TLM" w l ,\
     "v3.tout" u 1:3 title "v3: GGI_TLM" w l 
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Ey.jpg"

plot "two_wire_diode_test.field.tout" u 1:3 title "Ey: GGI_TLM" w l   
#PAUSE pause -1
