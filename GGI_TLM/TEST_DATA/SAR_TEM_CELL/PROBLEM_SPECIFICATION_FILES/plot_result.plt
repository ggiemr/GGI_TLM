#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_18_SAR_TEM_CELL_field.jpg"

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V/m"
plot "TEM_cell.cable_current.tout" u 1:3 title "SAR test cable_current: GGI_TLM" w l 
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_18_SAR_TEM_CELL_SAR.jpg"
set autoscale x
set autoscale y

set xlabel "SAR output number"
set ylabel "SAR (W/kg)"
plot "TEM_cell.SAR.fout" u 2:3 title "GGI_TLM" w p
#PAUSE pause -1
