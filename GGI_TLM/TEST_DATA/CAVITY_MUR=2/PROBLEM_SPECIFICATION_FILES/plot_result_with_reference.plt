#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_3_cavity_mur=2_ref.jpg"

set xlabel "Frequency (Hz)"
set ylabel "V/m"

plot "cavity_1_2_3.fout" u 1:5 title "Cavity E field: GGI_TLM" w l,\
"PROBLEM_SPECIFICATION_FILES/cavity_1_2_3.fout_ref" u 1:5 title "Cavity E field: GGI_TLM reference" w p,\
"PROBLEM_SPECIFICATION_FILES/cavity_1_2_3_TLM.field.fout" u 1:4 title "Cavity E field: Fieldsolve reference" w p ,\
"PROBLEM_SPECIFICATION_FILES/analytic_cavity_modes" u 4:5 title "Analytic cavity mode frequencies" w l

#PAUSE pause -1
