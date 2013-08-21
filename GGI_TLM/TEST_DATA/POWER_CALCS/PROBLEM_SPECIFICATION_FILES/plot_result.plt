#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_37_POWER_CALCS_field.jpg"
set xlabel "Time (s)"
set ylabel "V/m"
plot "power_calc.field.tout" u 1:3 title "E field: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_37_POWER_CALCS_freq_power.jpg"
set xlabel "Frequency (Hz)"
set ylabel "W"
plot "power_calc.frequency_domain_power_surface.fout" u 1:5 title "Power: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_37_POWER_CALCS_TCS.jpg"
set xlabel "Frequency (Hz)"
set ylabel "Transmission Cross Section dB square metres"
plot "TCS.fout" u 1:7 title "Cavity E field: GGI_TLM" w l
#PAUSE pause -1
