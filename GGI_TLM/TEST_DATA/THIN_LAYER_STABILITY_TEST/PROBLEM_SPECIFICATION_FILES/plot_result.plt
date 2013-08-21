#JPG set term jpeg

set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_7_THIN_LAYER_STABILITY_TEST_field.jpg"

set xlabel "Time (s)"
set ylabel "V/m"
plot "cavity.field.tout" u 1:3 title "Ey: GGI_TLM" w l 
#PAUSE pause -1
