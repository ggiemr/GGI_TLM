#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_42_YORK_BOX_STRAIGHT_WIRE_FWARP_current.jpg"
set xlabel "Time (s)"
set ylabel "A"
plot "york_box.cable_current.tout" u 1:3 title "York box current 3cm mesh: GGI_TLM" w p
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_42_YORK_BOX_STRAIGHT_WIRE_FWARP_S21.jpg"
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "S21 (dB)"
plot "PROBLEM_SPECIFICATION_FILES/york_box_original_S21_TLM.fout" u 1:2 title "S21: Fieldsolve Reference " w p,\
     "york_box_S21.fout_1" u 1:7 title "S21 :3cm mesh GGI_TLM" w l ,\
     "york_box_S21.fout" u 1:7 title "S21 :3cm mesh GGI_TLM, Fwarp" w l
     
  
#PAUSE pause -1




