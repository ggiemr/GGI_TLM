
#JPG set term jpeg


set xlabel "Time (s)"
set ylabel "A"

#OUTPUT_TO_FILE set output "monopole_current_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/monopole.cable_current.tout_ref" u 1:3 title "monopole current: GGI_TLM Reference" w p ,\
     "monopole.cable_current.tout" u 1:3 title "monopole current: GGI_TLM" w l 
     
#PAUSE pause -1

set autoscale x
set yrange[-3000:5000]

set title "monopole impedance"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"

#OUTPUT_TO_FILE set output "monopole_impedance_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/monopole.impedance.fout_ref" u 1:3 title "Re{Z}: GGI_TLM Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/monopole.impedance.fout_ref" u 1:4 title "Im{Z}: GGI_TLM Reference " w p,\
     "monopole.impedance.fout" u 1:3 title "Re{Z}: GGI_TLM " w l,\
     "monopole.impedance.fout" u 1:4 title "Im{Z}: GGI_TLM " w l

#PAUSE pause -1

