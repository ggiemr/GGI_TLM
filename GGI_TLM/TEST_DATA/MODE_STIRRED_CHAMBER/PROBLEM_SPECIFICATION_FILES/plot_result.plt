#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_35_MODE_STIRRED_CHAMBER_field.jpg"
set autoscale x
set autoscale y

set xlabel "Frequency (Hz)"
set ylabel "V/m"
plot "mode_stir1.field.fout" u 1:5 title "Cavity E field1: GGI_TLM" w l,\
"mode_stir2.field.fout" u 1:5 title "Cavity E field2: GGI_TLM" w l,\
"mode_stir3.field.fout" u 1:5 title "Cavity E field3: GGI_TLM" w l,\
"mode_stir4.field.fout" u 1:5 title "Cavity E field4: GGI_TLM" w l,\
"mode_stir5.field.fout" u 1:5 title "Cavity E field5: GGI_TLM" w l,\
"mode_stir6.field.fout" u 1:5 title "Cavity E field6: GGI_TLM" w l,\
"mode_stir7.field.fout" u 1:5 title "Cavity E field7: GGI_TLM" w l,\
"mode_stir8.field.fout" u 1:5 title "Cavity E field8: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_35_MODE_STIRRED_CHAMBER_field_sum.jpg"
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "V/m"
plot "mode_stir_sum.field.fout" u 1:5 title "Cavity E field sum: GGI_TLM" w l
#PAUSE pause -1
