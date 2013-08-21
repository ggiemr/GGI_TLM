#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_57_SE_BOX_SLOT_field.jpg"
set title "Internal E field"
set xlabel "Frequency (Hz)"
set ylabel "V/m"
plot "box_slot.field.tout" u 1:3 title "GGI_TLM" w l    
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_57_SE_BOX_SLOT_SE.jpg"
set title "Shielding Effectiveness"
set xlabel "Frequency (Hz)"
set ylabel "dB"
plot "Shielding_effectiveness.fout" u 1:7 title "GGI_TLM" w l   
#PAUSE pause -1
