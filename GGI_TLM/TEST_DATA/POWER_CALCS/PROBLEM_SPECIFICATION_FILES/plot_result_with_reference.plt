#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_37_POWER_CALCS_field_ref.jpg"
set xlabel "Time (s)"
set ylabel "V/m"
plot "power_calc.field.tout" u 1:3 title "E field: GGI_TLM" w l,\
"PROBLEM_SPECIFICATION_FILES/power_calc.field.tout_ref" u 1:3 title "E field: GGI_TLM reference" w p
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_37_POWER_CALCS_freq_power_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "W"
plot "power_calc.frequency_domain_power_surface.fout" u 1:5 title "Power: GGI_TLM" w l,\
"PROBLEM_SPECIFICATION_FILES/power_calc.frequency_domain_power_surface.fout_ref" u 1:5 title "Power: GGI_TLM reference" w p
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_37_POWER_CALCS_TCS_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "Transmission Cross Section dB square metres"
plot "TCS.fout" u 1:7 title "Cavity E field: GGI_TLM" w l,\
"PROBLEM_SPECIFICATION_FILES/TCS.fout_ref" u 1:7 title "Cavity E field: GGI_TLM reference" w p
#PAUSE pause -1
