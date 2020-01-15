#JPG set term jpeg

set autoscale x
set autoscale y
set grid
set xlabel "Time (s)"
set ylabel "V/m"

#OUTPUT_TO_FILE set output "Test_case_8_loop_field_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/loop.field.tout_ref" u 1:3 title "loop field: GGI_TLM Reference" w p ,\
     "loop.field.tout" u 1:3 title "loop field: GGI_TLM" w l 
     
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
##OUTPUT_TO_FILE set output "Test_case_8_loop_far_field_ref.jpg"
#
#plot  "PROBLEM_SPECIFICATION_FILES/loop.far_field.fout_ref" u 1:3 title "loop Etheta: GGI_TLM Reference" w p ,\
#      "PROBLEM_SPECIFICATION_FILES/loop.far_field.fout_ref" u 1:4 title "loop Ephi: GGI_TLM Reference" w p ,\
#      "loop.far_field.fout" u 1:3 title "Etheta: GGI_TLM" w l,\
#      "loop.far_field.fout" u 1:4 title "Ephi: GGI_TLM" w l
#     
##PAUSE pause -1

set xlabel "Time (s)"
set ylabel "A"

#OUTPUT_TO_FILE set output "Test_case_8_loop_current_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/loop.cable_current.tout_ref" u 1:3 title "loop current: GGI_TLM Reference" w p ,\
     "loop.cable_current.tout" u 1:3 title "loop current: GGI_TLM" w l 
     
#PAUSE pause -1

set autoscale x
set yrange[-3000:5000]

set title "loop impedance"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"

#OUTPUT_TO_FILE set output "Test_case_8_loop_impedance_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/loop.impedance.fout_ref" u 1:3 title "Re{Z}: GGI_TLM Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/loop.impedance.fout_ref" u 1:4 title "Im{Z}: GGI_TLM Reference " w p,\
     "loop.impedance.fout" u 1:3 title "Re{Z}: GGI_TLM " w l,\
     "loop.impedance.fout" u 1:4 title "Im{Z}: GGI_TLM " w l

#PAUSE pause -1

