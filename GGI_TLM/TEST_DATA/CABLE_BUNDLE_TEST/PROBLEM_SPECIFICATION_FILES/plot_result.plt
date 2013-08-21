#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_61_cable_bundle_test_current.jpg"

set title "Cable Current"
set xlabel "Time (s)"
set ylabel "A"
plot "cable_test.cable_current.tout" u 1:3 title "Cable current: GGI_TLM" w p
#PAUSE pause -1
