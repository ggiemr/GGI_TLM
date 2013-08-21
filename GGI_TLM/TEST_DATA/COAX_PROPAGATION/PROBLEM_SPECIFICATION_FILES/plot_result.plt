#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_44_COAX_PROPAGATION.jpg"

set title "Cable current"

set xlabel "Time (s)"
set ylabel "A"
plot "coax.cable_current.tout" u 1:3 w p
#PAUSE pause -1
