#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_18_SAR_TEM_CELL_field_ref.jpg"

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V/m"
plot "PROBLEM_SPECIFICATION_FILES/TEM_cell.cable_current.tout_ref" u 1:3 title "SAR test cable_current: GGI_TLM Reference" w p ,\
     "TEM_cell.cable_current.tout" u 1:3 title "SAR test cable_current: GGI_TLM" w l 
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_18_SAR_TEM_CELL_SAR_ref.jpg"
set autoscale x
set autoscale y

set xlabel "SAR output number"
set ylabel "SAR (W/kg)"
plot "PROBLEM_SPECIFICATION_FILES/TEM_cell.SAR.fout_ref" u 2:3 title "SAR : GGI_TLM Reference" w p ,\
     "TEM_cell.SAR.fout" u 2:3 title "SAR : GGI_TLM" w l 
#PAUSE pause -1
