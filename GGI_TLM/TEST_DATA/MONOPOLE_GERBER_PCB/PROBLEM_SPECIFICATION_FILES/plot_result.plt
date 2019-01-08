#JPG set term jpeg
#___________________________________________________________________#

set xlabel "Time (s)"
set ylabel "A"

#OUTPUT_TO_FILE set output "monopole_current.jpg"

plot "monopole.cable_current.tout" u 1:3 title "monopole current: GGI_TLM" w l 
     
#PAUSE pause -1

#___________________________________________________________________#

set autoscale x
set yrange[-3000:5000]

set title "monopole impedance"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"
set grid

#OUTPUT_TO_FILE set output "monopole_impedance.jpg"

plot "monopole.impedance.fout" u 1:3 title "Re{Z}: GGI_TLM " w l,\
     "monopole.impedance.fout" u 1:4 title "Im{Z}: GGI_TLM " w l

#PAUSE pause -1

#___________________________________________________________________#
