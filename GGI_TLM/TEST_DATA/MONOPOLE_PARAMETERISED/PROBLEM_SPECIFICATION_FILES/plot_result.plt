#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_9_MONOPOLE_current_ref.jpg"

set title "Monopole current"
set xlabel "Time (s)"
set ylabel "A"
plot "monopole.cable_current.tout" u 1:3 w l
#PAUSE pause -1


set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_9_MONOPOLE_field_ref.jpg"

set title "Monopole field"
set xlabel "Time (s)"
set ylabel "V/m"
plot "monopole.field.tout" u 1:3 w l
#PAUSE pause -1
