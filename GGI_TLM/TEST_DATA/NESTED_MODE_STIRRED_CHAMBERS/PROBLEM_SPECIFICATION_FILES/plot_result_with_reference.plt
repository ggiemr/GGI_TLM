#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_36_NESTED_MODE_STIRRED_CHAMBER_field_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "V/m"
plot "nested_RC.field.fout" u 1:5 title "Cavity E field: GGI_TLM" w l,\
"PROBLEM_SPECIFICATION_FILES/nested_RC.field.fout_ref" u 1:5 title "Cavity E field: GGI_TLM reference" w p

#PAUSE pause -1

