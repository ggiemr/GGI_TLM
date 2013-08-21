#JPG set term jpeg

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V/m"

#OUTPUT_TO_FILE set output "Test_case_8_dipole_field_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/dipole_original_TLM.field.tout" u 1:3 title "Dipole field: Fieldsolve Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/dipole.field.tout_ref" u 1:3 title "Dipole field: GGI_TLM Reference" w p ,\
     "dipole.field.tout" u 1:3 title "Dipole field: GGI_TLM" w l 
     
#PAUSE pause -1

#___________________________________________________________________#

set autoscale x
set autoscale y

set xlabel "Angle (degrees)"
set ylabel "Normalised Far field"

#OUTPUT_TO_FILE set output "Test_case_8_dipole_far_field_ref.jpg"

plot  "PROBLEM_SPECIFICATION_FILES/dipole_original_TLM.ffout" u 1:3 title "Dipole Etheta: Fieldsolve Reference " w p,\
      "PROBLEM_SPECIFICATION_FILES/dipole_original_TLM.ffout" u 1:4 title "Dipole Ephi: Fieldsolve Reference " w p,\
      "PROBLEM_SPECIFICATION_FILES/dipole.far_field.fout_ref" u 1:3 title "Dipole Etheta: GGI_TLM Reference" w p ,\
      "PROBLEM_SPECIFICATION_FILES/dipole.far_field.fout_ref" u 1:4 title "Dipole Ephi: GGI_TLM Reference" w p ,\
      "dipole.far_field.fout" u 1:3 title "Etheta: GGI_TLM" w l,\
      "dipole.far_field.fout" u 1:4 title "Ephi: GGI_TLM" w l
     
#PAUSE pause -1

#___________________________________________________________________#

set xlabel "Time (s)"
set ylabel "A"

#OUTPUT_TO_FILE set output "Test_case_8_dipole_current_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/dipole_original_TLM.cable_current.tout" u 1:3 title "Dipole current: Fieldsolve Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/dipole.cable_current.tout_ref" u 1:3 title "Dipole current: GGI_TLM Reference" w p ,\
     "dipole.cable_current.tout" u 1:3 title "Dipole current: GGI_TLM" w l 
     
#PAUSE pause -1

#___________________________________________________________________#

set autoscale x
set yrange[-3000:5000]

set title "Dipole impedance"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"

#OUTPUT_TO_FILE set output "Test_case_8_dipole_impedance_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/dipole_original_TLM.impedance.fout" u 1:2 title "Re{Z}: Fieldsolve Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/dipole_original_TLM.impedance.fout" u 1:3 title "Im{Z}: Fieldsolve Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/dipole.impedance.fout_ref" u 1:3 title "Re{Z}: GGI_TLM Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/dipole.impedance.fout_ref" u 1:4 title "Im{Z}: GGI_TLM Reference " w p,\
     "dipole.impedance.fout" u 1:3 title "Re{Z}: GGI_TLM " w l,\
     "dipole.impedance.fout" u 1:4 title "Im{Z}: GGI_TLM " w l

#PAUSE pause -1

#___________________________________________________________________#


