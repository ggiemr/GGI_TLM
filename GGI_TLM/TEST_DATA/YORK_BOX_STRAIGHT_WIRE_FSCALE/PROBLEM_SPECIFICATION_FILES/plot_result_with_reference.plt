#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_41_YORK_BOX_STRAIGHT_WIRE_FSCALE_current_ref.jpg"
set xlabel "Time (s)"
set ylabel "A"
plot "PROBLEM_SPECIFICATION_FILES/york_box_original_TLM.cable_current.tout" u 1:3 title "Port 2 current: Fieldsolve Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/york_box.cable_current.tout_ref" u 1:3 title "Port 2 current: GGI_TLM Reference" w p ,\
     "york_box.cable_current.tout" u 1:3 title "Port 2 current: GGI_TLM" w l 
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_41_YORK_BOX_STRAIGHT_WIRE_FSCALE_S21_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "S21 (dB)"
plot "PROBLEM_SPECIFICATION_FILES/york_box_original_S21_TLM.fout" u 1:2 title "S21: Fieldsolve Reference " w p,\
     "PROBLEM_SPECIFICATION_FILES/york_box_S21.fout_ref" u 1:7 title "S21: GGI_TLM Fscale reference" w p,\
     "york_box_S21.fout_1" u 1:7 title "S21: GGI_TLM, no scale" w l ,\
     "york_box_S21.fout" u 1:7 title "S21: GGI_TLM, Fscale" w l 
    
#PAUSE pause -1

