#JPG set term jpeg

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V/m"

#OUTPUT_TO_FILE set output "Test_case_8_helix_field.jpg"

plot "helix.field.tout" u 1:3 title "helix field: GGI_TLM" w l 
     
#PAUSE pause -1

#___________________________________________________________________#


#
# FAR FIELD PLOT NOW DONE IN VTK FORMAT WITH PARAVIEW
#

#set autoscale x
#set autoscale y
#
#set xlabel "Angle (degrees)"
#set ylabel "Normalised Far field"
#
##OUTPUT_TO_FILE set output "Test_case_8_helix_far_field.jpg"
#
#plot  "helix.far_field.fout" u 1:3 title "Etheta: GGI_TLM" w l,\
#      "helix.far_field.fout" u 1:4 title "Ephi: GGI_TLM" w l
#    
##PAUSE pause -1

#___________________________________________________________________#

set xlabel "Time (s)"
set ylabel "A"

#OUTPUT_TO_FILE set output "Test_case_8_helix_current.jpg"

plot "helix.cable_current.tout" u 1:3 title "helix current: GGI_TLM" w l 
     
#PAUSE pause -1

#___________________________________________________________________#

set autoscale x
set yrange[-3000:5000]

set title "helix impedance"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"

#OUTPUT_TO_FILE set output "Test_case_8_helix_impedance.jpg"

plot "helix.impedance.fout" u 1:3 title "Re{Z}: GGI_TLM " w l,\
     "helix.impedance.fout" u 1:4 title "Im{Z}: GGI_TLM " w l

#PAUSE pause -1

#___________________________________________________________________#


