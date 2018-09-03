#JPG set term jpeg

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V/m"

#OUTPUT_TO_FILE set output "Test_case_8_helix_field_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/helix.field.tout_ref" u 1:3 title "helix field: GGI_TLM Reference" w p ,\
     "helix.field.tout" u 1:3 title "helix field: GGI_TLM" w l 
     
#PAUSE pause -1

#
# FAR FIELD PLOT NOW DONE IN VTK FORMAT WITH PARAVIEW
#

#set autoscale x
#set autoscale y
#
#set xlabel "Angle (degrees)"
#set ylabel "Normalised Far field"
#
##OUTPUT_TO_FILE set output "Test_case_8_helix_far_field_ref.jpg"
#
#plot  "PROBLEM_SPECIFICATION_FILES/helix.far_field.fout_ref" u 1:3 title "helix Etheta: GGI_TLM Reference" w p ,\
#      "PROBLEM_SPECIFICATION_FILES/helix.far_field.fout_ref" u 1:4 title "helix Ephi: GGI_TLM Reference" w p ,\
#      "helix.far_field.fout" u 1:3 title "Etheta: GGI_TLM" w l,\
#      "helix.far_field.fout" u 1:4 title "Ephi: GGI_TLM" w l
#     
##PAUSE pause -1

set xlabel "Time (s)"
set ylabel "A"

#OUTPUT_TO_FILE set output "Test_case_8_helix_current_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/helix.cable_current.tout_ref" u 1:3 title "helix current: GGI_TLM Reference" w p ,\
     "helix.cable_current.tout" u 1:3 title "helix current: GGI_TLM" w l 
     
#PAUSE pause -1

set autoscale x
set yrange[-3000:5000]

set title "helix impedance"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"

#OUTPUT_TO_FILE set output "Test_case_8_helix_impedance_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/helix.impedance.fout_ref" u 1:3 title "Re{Z}: GGI_TLM Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/helix.impedance.fout_ref" u 1:4 title "Im{Z}: GGI_TLM Reference " w p,\
     "helix.impedance.fout" u 1:3 title "Re{Z}: GGI_TLM " w l,\
     "helix.impedance.fout" u 1:4 title "Im{Z}: GGI_TLM " w l

#PAUSE pause -1

