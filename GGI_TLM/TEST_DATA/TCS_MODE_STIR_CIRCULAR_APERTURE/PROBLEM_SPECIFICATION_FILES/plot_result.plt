#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_40_TCS_MODE_STIR_field.jpg"
set xlabel "Time (s)"
set ylabel "V/m"
plot "TCS_mode_stir.field.tout" u 1:3 title "Cavity E field: GGI_TLM" w l,\
"TCS_mode_stir.volume_average_field.tout" u 1:3 title "Cavity E field: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_40_TCS_MODE_STIR_freq_field.jpg"
set xlabel "Frequency (Hz)"
set ylabel "dB V/m"
plot "E_cavity.fout" u 1:7 title "Cavity E field: GGI_TLM" w l,\
"E_cavity_volume_average.fout" u 1:7 title "Cavity E field: GGI_TLM" w l,\
"E_cavity_freq_stir.fout" u 1:7 title "Cavity E field: GGI_TLM freq_stir" w l,\
"E_cavity_volume_average_freq_stir.fout" u 1:7 title "Cavity E field: GGI_TLM freq_stir" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_40_TCS_MODE_STIR_power.jpg"
set xlabel "Frequency (Hz)"
set ylabel "W"
plot "TCS_mode_stir.frequency_domain_power_surface.fout" u 1:3 title "Transmitted Power: GGI_TLM" w l,\
"TCS_mode_stir.frequency_domain_power_surface_freq_stir.fout" u 1:3 title "Transmitted Power: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_40_TCS_MODE_STIR_power_dB.jpg"
set xlabel "Frequency (Hz)"
set ylabel "dB"
plot "TCS_mode_stir.frequency_domain_power_surface.fout" u 1:7 title "Transmitted Power: GGI_TLM" w l,\
"TCS_mode_stir.frequency_domain_power_surface_freq_stir.fout" u 1:7 title "Transmitted Power: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_40_TCS_MODE_STIR_TCS.jpg"
set xlabel "Frequency (Hz)"
set ylabel "Transmission Cross Section dB square metres"
set xrange [0:10e9]
plot "TCS.fout" u 1:7 title "Cavity E field: GGI_TLM" w l,\
"TCS_volume_average.fout" u 1:7 title "Cavity E field: GGI_TLM volume average" w l,\
"TCS_freq_stir.fout" u 1:7 title "Cavity E field: GGI_TLM freq_stir" w l,\
"TCS_volume_average_freq_stir.fout" u 1:7 title "Cavity E field: GGI_TLM volume average freq_stir" w l,\
"PROBLEM_SPECIFICATION_FILES/TCS_analytic.fout" u 1:2 title "Analytic formula" w p
#PAUSE pause -1

set xlabel "Frequency (Hz)"
set ylabel "Transmission Cross Section dB square metres"
set xrange [0:10e9]
plot "TCS_freq_stir.fout" u 1:7 title "GGI_TLM freq_stir" w l,\
"PROBLEM_SPECIFICATION_FILES/TCS_analytic.fout" u 1:2 title "Analytic formula" w p
#PAUSE pause -1
