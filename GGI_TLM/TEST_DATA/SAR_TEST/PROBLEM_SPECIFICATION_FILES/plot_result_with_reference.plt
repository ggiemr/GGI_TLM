#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_17_SAR_TEST_field_ref.jpg"

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V/m"
plot "PROBLEM_SPECIFICATION_FILES/SAR_test_original_TLM.field.tout" u 1:3 title "SAR test field: Fieldsolve Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/SAR_test.field.tout_ref" u 1:3 title "SAR test field: GGI_TLM Reference" w p ,\
     "SAR_test.field.tout" u 1:3 title "SAR test field: GGI_TLM" w l 
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_17_SAR_TEST_SAR_ref.jpg"

set autoscale x
set yrange[1.0E-029:1.5E-029]

set xlabel "SAR output number"
set ylabel "SAR (W/kg)"
plot "PROBLEM_SPECIFICATION_FILES/SAR_test_original_TLM.SAR" u 2:3 title "SAR : Fieldsolve Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/SAR_test.SAR.fout_ref" u 2:3 title "SAR : GGI_TLM Reference" w p ,\
     "SAR_test.SAR.fout" u 2:3 title "SAR : GGI_TLM" w l 
#PAUSE pause -1
