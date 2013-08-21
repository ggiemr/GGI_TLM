#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_51_WIRE_OVER_GROUND_voltage_1_2_ref.jpg"
set logscale x
set logscale y
set xlabel "Frequency (Hz)"
set ylabel "mV"
plot "PROBLEM_SPECIFICATION_FILES/voltage_1_norm.fout_ref" u 1:5 title "GGI_TLM: voltage 1 reference" w p,\
     "PROBLEM_SPECIFICATION_FILES/voltage_2_norm.fout_ref" u 1:5 title "GGI_TLM: voltage 2 reference" w p,\
     "voltage_1_norm.fout" u 1:5 title "GGI_TLM: voltage 1" w l,\
     "voltage_2_norm.fout" u 1:5 title "GGI_TLM: voltage 2" w l
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_51_WIRE_OVER_GROUND_voltage_3_4_ref.jpg"
set logscale x
set logscale y
set xlabel "Frequency (Hz)"
set ylabel "mV"
plot "PROBLEM_SPECIFICATION_FILES/voltage_3_norm.fout_ref" u 1:5 title "GGI_TLM: voltage 3 reference" w p,\
     "PROBLEM_SPECIFICATION_FILES/voltage_4_norm.fout_ref" u 1:5 title "GGI_TLM: voltage 4 reference" w p,\
     "voltage_3_norm.fout" u 1:5 title "GGI_TLM: voltage 3" w l,\
     "voltage_4_norm.fout" u 1:5 title "GGI_TLM: voltage 4" w l
#PAUSE pause -1
