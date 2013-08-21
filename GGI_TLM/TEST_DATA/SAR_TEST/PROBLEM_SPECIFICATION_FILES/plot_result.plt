#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_17_SAR_TEST_field.jpg"

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V/m"
plot "SAR_test.field.tout" u 1:3 title "SAR test field: GGI_TLM" w l 
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_17_SAR_TEST_SAR.jpg"

set autoscale x
set autoscale y

set xlabel "SAR output number"
set ylabel "SAR (W/kg)"
plot "SAR_test.SAR.fout" u 2:3 title "GGI_TLM" w l 
#PAUSE pause -1
