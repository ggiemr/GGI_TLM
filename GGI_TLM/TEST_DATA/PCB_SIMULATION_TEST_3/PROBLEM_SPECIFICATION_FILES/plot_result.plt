#JPG set term jpeg

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V"

#OUTPUT_TO_FILE set output "Source_voltage.jpg"

plot "source_voltage1.tout" u 1:3 title "source voltage1: GGI_TLM" w l ,\
     "load_voltage1.tout" u 1:3 title "load voltage1: GGI_TLM" w l ,\
    "source_voltage2.tout" u 1:3 title "source voltage2: GGI_TLM" w l ,\
     "load_voltage2.tout" u 1:3 title "load voltage2: GGI_TLM" w l 
     
#PAUSE pause -1
