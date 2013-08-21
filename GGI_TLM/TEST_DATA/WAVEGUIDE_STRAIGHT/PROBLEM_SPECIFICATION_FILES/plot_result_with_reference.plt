#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_21_WAVEGUIDE_STRAIGHT_mode_ref.jpg"
set xlabel "Time (s)"
set ylabel "V/m"
plot "waveguide_straight.mode.tout" u 1:3 title "Waveguide mode field: GGI_TLM" w p,\
"PROBLEM_SPECIFICATION_FILES/waveguide_straight.mode.tout_ref" u 1:3 title "Waveguide mode field: GGI_TLM reference" w p
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_21_WAVEGUIDE_STRAIGHT_field_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot "waveguide_straight.field.tout" u 1:3 title "Waveguide E field: GGI_TLM" w p,\
"PROBLEM_SPECIFICATION_FILES/waveguide_straight.field.tout_ref" u 1:3 title "Waveguide E field: GGI_TLM reference" w p
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_21_WAVEGUIDE_STRAIGHT_Z_ref.jpg"
set autoscale x
set autoscale y
set title "Wave impedance"
set xlabel "Frequency (Hz)"
set ylabel "Ohms"
plot "PROBLEM_SPECIFICATION_FILES/waveguide_straight_original.fout.Z" u 1:2 title "Re{Z} :Fieldsolve" w p,\
     "PROBLEM_SPECIFICATION_FILES/waveguide_straight_original.fout.Z" u 1:3 title "Im{Z} :Fieldsolve" w p,\
     "PROBLEM_SPECIFICATION_FILES/waveguide_straight.fout_ref" u 1:4 title "Re{Z} :GGI_TLM reference" w p,\
     "PROBLEM_SPECIFICATION_FILES/waveguide_straight.fout_ref" u 1:5 title "Im{Z} :GGI_TLM reference" w p,\
     "waveguide_straight.fout" u 1:4 title "Re{Z} :GGI_TLM" w l,\
     "waveguide_straight.fout" u 1:5 title "Im{Z} :GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_21_WAVEGUIDE_STRAIGHT_beta_ref.jpg"
set autoscale x
set autoscale y
set title "Propagation constant"
set xlabel "Frequency (Hz)"
set ylabel "m^-1"
plot "PROBLEM_SPECIFICATION_FILES/waveguide_straight_original.fout.beta" u 1:2 title "Re{Beta} :Fieldsolve" w p,\
     "PROBLEM_SPECIFICATION_FILES/waveguide_straight_original.fout.beta" u 1:3 title "Im{Beta} :Fieldsolve" w p,\
     "PROBLEM_SPECIFICATION_FILES/waveguide_straight.fout_ref" u 1:2 title "Re{Z} :GGI_TLM reference" w p,\
     "PROBLEM_SPECIFICATION_FILES/waveguide_straight.fout_ref" u 1:3 title "Im{Z} :GGI_TLM reference" w p,\
     "waveguide_straight.fout" u 1:2 title "Re{Beta} :GGI_TLM" w l,\
     "waveguide_straight.fout" u 1:3 title "Im{Beta} :GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_21_WAVEGUIDE_STRAIGHT_R_ref.jpg"
set autoscale x
set autoscale y
set title "Reflection Coefficient"
set xlabel "Frequency (Hz)"
set ylabel "R"
plot "PROBLEM_SPECIFICATION_FILES/waveguide_straight_original.fout.R" u 1:2 title "Re{R} :Fieldsolve" w p,\
     "PROBLEM_SPECIFICATION_FILES/waveguide_straight_original.fout.R" u 1:3 title "Im{R} :Fieldsolve" w p,\
     "PROBLEM_SPECIFICATION_FILES/waveguide_straight.fout_ref" u 1:6 title "Re{Z} :GGI_TLM reference" w p,\
     "PROBLEM_SPECIFICATION_FILES/waveguide_straight.fout_ref" u 1:7 title "Im{Z} :GGI_TLM reference" w p,\
     "waveguide_straight.fout" u 1:6 title "Re{R} :GGI_TLM" w l,\
     "waveguide_straight.fout" u 1:7 title "Im{R} :GGI_TLM" w l
#PAUSE pause -1


