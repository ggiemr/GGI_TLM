#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_32_CAVITY_MODE_STIR_field.jpg"
set xlabel "Frequency (Hz)"
set ylabel "V/m"
plot "cavity_1_2_3.fout" u 1:5 title "Cavity E field: GGI_TLM" w l,\
"REFERENCE_RESULT/analytic_cavity_modes" u 4:5 title "Analytic cavity mode frequencies" w l
#PAUSE pause -1
