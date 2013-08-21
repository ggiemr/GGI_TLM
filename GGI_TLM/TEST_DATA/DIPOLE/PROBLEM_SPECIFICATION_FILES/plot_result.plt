#JPG set term jpeg

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V/m"

#OUTPUT_TO_FILE set output "Test_case_8_dipole_field.jpg"

plot "dipole.field.tout" u 1:3 title "Dipole field: GGI_TLM" w l 
     
#PAUSE pause -1

#___________________________________________________________________#

set autoscale x
set autoscale y

set xlabel "Angle (degrees)"
set ylabel "Normalised Far field"

#OUTPUT_TO_FILE set output "Test_case_8_dipole_far_field.jpg"

plot  "dipole.far_field.fout" u 1:3 title "Etheta: GGI_TLM" w l,\
      "dipole.far_field.fout" u 1:4 title "Ephi: GGI_TLM" w l
     
#PAUSE pause -1

#___________________________________________________________________#

set xlabel "Time (s)"
set ylabel "A"

#OUTPUT_TO_FILE set output "Test_case_8_dipole_current.jpg"

plot "dipole.cable_current.tout" u 1:3 title "Dipole current: GGI_TLM" w l 
     
#PAUSE pause -1

#___________________________________________________________________#

set autoscale x
set yrange[-3000:5000]

set title "Dipole impedance"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"

#OUTPUT_TO_FILE set output "Test_case_8_dipole_impedance.jpg"

plot "dipole.impedance.fout" u 1:3 title "Re{Z}: GGI_TLM " w l,\
     "dipole.impedance.fout" u 1:4 title "Im{Z}: GGI_TLM " w l

#PAUSE pause -1

#___________________________________________________________________#


