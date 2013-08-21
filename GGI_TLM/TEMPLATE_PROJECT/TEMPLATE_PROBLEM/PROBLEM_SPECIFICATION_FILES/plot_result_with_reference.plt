#JPG set term jpeg

#OUTPUT_TO_FILE set output "problem_name_current_ref.jpg"
set xlabel "Time (s)"
set ylabel "A"
plot "PROBLEM_SPECIFICATION_FILES/problem_name.cable_current.tout_ref" u 1:3 title "Port 2 current: GGI_TLM Reference" w p ,\
     "problem_name.cable_current.tout" u 1:3 title "Port 2 current: GGI_TLM" w l 
#PAUSE pause -1

#OUTPUT_TO_FILE set output "problem_name_S21_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "S21 (dB)"
plot "PROBLEM_SPECIFICATION_FILES/problem_name_S21.fout_ref" u 1:7 title "S21: GGI_TLM reference" w p,\
     "problem_name_S21.fout" u 1:7 title "S21: GGI_TLM" w l 
#PAUSE pause -1

