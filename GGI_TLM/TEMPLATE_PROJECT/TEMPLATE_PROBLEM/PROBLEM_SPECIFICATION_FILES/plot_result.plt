#JPG set term jpeg

#OUTPUT_TO_FILE set output "problem_name_wire_current.jpg"
set xlabel "Time (s)"
set ylabel "A"
plot "problem_name.cable_current.tout" u 1:3 title "Port 2 current: GGI_TLM" w l 
#PAUSE pause -1

#OUTPUT_TO_FILE set output "problem_name_S21.jpg"
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "S21 (dB)"
plot  "problem_name_S21.fout" u 1:7 title "S21: GGI_TLM" w l 
#PAUSE pause -1


