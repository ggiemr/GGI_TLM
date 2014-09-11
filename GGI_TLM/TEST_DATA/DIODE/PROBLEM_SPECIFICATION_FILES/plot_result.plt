#JPG set term jpeg

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "A"

#OUTPUT_TO_FILE set output "Test_case_DIODE.jpg"

plot "microstrip.cable_current.tout" u 1:3 title "microstrip current: GGI_TLM" w l 
     
#PAUSE pause -1
