#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_27_cable_test_current.jpg"

set title "Cable Current"
set xlabel "Time (s)"
set ylabel "A"
plot "cable_test.cable_current.tout" u 1:3 title "Cable current: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_27_cable_test_field.jpg"
set title "Electric Field"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot "cable_test.field.tout" u 1:3 title "E field: GGI_TLM" w l
#PAUSE pause -1
