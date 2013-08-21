#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_30_monopole_current.jpg"
set xlabel "Time (s)"
set ylabel "A"
plot "coax_monopole.cable_current.tout" u 1:3 title "Coax_monopole current: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_30_monopole_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot "coax_monopole.field.tout" u 1:3 title "Coax_monopole field: GGI_TLM" w l
#PAUSE pause -1
