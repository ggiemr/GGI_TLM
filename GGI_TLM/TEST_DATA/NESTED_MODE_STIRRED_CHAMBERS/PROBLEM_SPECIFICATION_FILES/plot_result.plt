#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_36_NESTED_MODE_STIRRED_CHAMBER_time_field.jpg"
set xlabel "Time (s)"
set ylabel "V/m"
plot "nested_RC.field.tout" u 1:3 title "Cavity E field: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_36_NESTED_MODE_STIRRED_CHAMBER_freq_field.jpg"
set xlabel "Frequency (Hz)"
set ylabel "V/m"
plot "nested_RC.field.fout" u 1:5 title "Cavity E field: GGI_TLM" w l
#PAUSE pause -1
