#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_57_SE_BOX_SLOT_SE_ref.jpg"
set title "Shielding Effectiveness"
set xlabel "Frequency (Hz)"
set ylabel "dB"
plot "PROBLEM_SPECIFICATION_FILES/Shielding_effectiveness.fout_ref" u 1:7 title "GGI_TLM reference" w l ,\
     "Shielding_effectiveness.fout" u 1:7 title "GGI_TLM" w p
     
#PAUSE pause -1
