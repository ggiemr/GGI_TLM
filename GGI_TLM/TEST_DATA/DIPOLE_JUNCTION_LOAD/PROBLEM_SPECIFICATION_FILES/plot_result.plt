#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_31_dipole_current.jpg"
set xlabel "Time (s)"
set ylabel "A"
plot "dipole.cable_current.tout" u 1:3 title "Dipole current: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_31_dipole_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot "dipole.field.tout" u 1:3 title "Dipole field: GGI_TLM" w l
#PAUSE pause -1
