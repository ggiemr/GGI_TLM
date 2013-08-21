#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_39_TCS_MODE_STIR_field_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "V/m"
plot "E_cavity.fout" u 1:5 title "Cavity E field: GGI_TLM" w l,\
"PROBLEM_SPECIFICATION_FILES/E_cavity.fout_ref" u 1:5 title "Cavity E field: GGI_TLM reference" w p
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_39_TCS_MODE_STIR_power_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "W"
plot "TCS_mode_stir.frequency_domain_power_surface.fout" u 1:7 title "Transmitted Power: GGI_TLM" w l,\
"TCS_mode_stir.frequency_domain_power_surface.fout_ref" u 1:7 title "Transmitted Power: GGI_TLM reference" w p
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_39_TCS_MODE_STIR_TCS_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "Transmission Cross Section dB square metres"
plot "TCS.fout" u 1:7 title "Transmission cross section: GGI_TLM" w l,\
"PROBLEM_SPECIFICATION_FILES/TCS.fout_ref" u 1:7 title "Transmission cross section: GGI_TLM reference" w p

#PAUSE pause -1
