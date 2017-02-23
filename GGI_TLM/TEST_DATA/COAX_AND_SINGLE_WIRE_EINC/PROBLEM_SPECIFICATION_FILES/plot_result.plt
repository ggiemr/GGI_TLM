#JPG set term jpeg

#OUTPUT_TO_FILE set output "source_current.jpg"

set xlabel "time (s)"
set ylabel "Current (A)"
plot "coax_and_wire.cable_current.tout" u 1:3 w l
#PAUSE pause -1

##OUTPUT_TO_FILE set output "excitation.jpg"
#
#set xlabel "Frequency (Hz)"
#set ylabel "V/m"
#plot "excitation.fout" u 1:5 w l
##PAUSE pause -1

#OUTPUT_TO_FILE set output "source_current.jpg"

set xlabel "Frequency (Hz)"
set ylabel "A"
plot "source_current.fout" u 1:5 w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "normalised_source_voltage.jpg"

set logscale x
set xlabel "Frequency (Hz)"
set ylabel "V(dB)"
plot "normalised_source_voltage.fout" u 1:7 w l
#PAUSE pause -1
