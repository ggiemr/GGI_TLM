#JPG set term jpeg

set autoscale x
set autoscale y
set grid
set xlabel "Time (s)"
set ylabel "A/m"

#OUTPUT_TO_FILE set output "Test_case_8_loop_field.jpg"

plot "loop.field.tout" u 1:3 title "loop field: GGI_TLM" w l 
     
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
##OUTPUT_TO_FILE set output "Test_case_8_loop_far_field.jpg"
#
#plot  "loop.far_field.fout" u 1:3 title "Etheta: GGI_TLM" w l,\
#      "loop.far_field.fout" u 1:4 title "Ephi: GGI_TLM" w l
#    
##PAUSE pause -1

#___________________________________________________________________#

set xlabel "Time (s)"
set ylabel "A"

#OUTPUT_TO_FILE set output "Test_case_8_loop_current.jpg"

plot "loop.cable_current.tout" u 1:3 title "loop current: GGI_TLM" w l 
     
#PAUSE pause -1

#___________________________________________________________________#

set autoscale x
set yrange[-3000:5000]

set title "loop impedance"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"

#OUTPUT_TO_FILE set output "Test_case_8_loop_impedance.jpg"

plot "loop.impedance.fout" u 1:3 title "Re{Z}: GGI_TLM " w l,\
     "loop.impedance.fout" u 1:4 title "Im{Z}: GGI_TLM " w l

#PAUSE pause -1

#___________________________________________________________________#


