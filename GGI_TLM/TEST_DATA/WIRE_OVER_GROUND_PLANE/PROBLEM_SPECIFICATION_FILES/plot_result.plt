#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_51_WIRE_OVER_GROUND_voltage_1_2.jpg"
set logscale x
set logscale y
set xlabel "Frequency (Hz)"
set ylabel "mV"
plot "voltage_1_norm.fout" u 1:5 title "GGI_TLM: voltage 1" w l,\
     "voltage_2_norm.fout" u 1:5 title "GGI_TLM: voltage 2" w l
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_51_WIRE_OVER_GROUND_voltage_3_4.jpg"
set logscale x
set logscale y
set xlabel "Frequency (Hz)"
set ylabel "mV"
plot "voltage_3_norm.fout" u 1:5 title "GGI_TLM: voltage 3" w l,\
     "voltage_4_norm.fout" u 1:5 title "GGI_TLM: voltage 4" w l
#PAUSE pause -1
