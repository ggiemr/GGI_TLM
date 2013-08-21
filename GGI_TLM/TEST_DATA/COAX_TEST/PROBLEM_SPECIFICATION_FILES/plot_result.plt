#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_26_coax_test_current.jpg"
set xlabel "Time (s)"
set ylabel "A"
plot "coax.cable_current.tout" u 1:3 w p
#PAUSE pause -1

##OUTPUT_TO_FILE set output "Test_case_26_coax_test_field.jpg"
#set autoscale x
#set autoscale y
#set xlabel "Time (s)"
#set ylabel "V/m"
#plot "coax.field.tout" u 1:3 w l
##PAUSE pause -1
