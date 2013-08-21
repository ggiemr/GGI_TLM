#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_22_WAVEGUIDE_STRAIGHT_ANIMATION_mode.jpg"
set xlabel "Time (s)"
set ylabel "V/m"
plot "waveguide_straight.mode.tout" u 1:3 title "Waveguide mode field: GGI_TLM" w p
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_22_WAVEGUIDE_STRAIGHT_ANIMATION_field.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "V/m"
plot "waveguide_straight.field.tout" u 1:3 title "Waveguide E field: GGI_TLM" w p
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_22_WAVEGUIDE_STRAIGHT_ANIMATION_Z.jpg"
set autoscale x
set autoscale y
set title "Wave impedance"
set xlabel "Frequency (Hz)"
set ylabel "Ohms"
plot "waveguide_straight.fout" u 1:4 title "Re{Z} :GGI_TLM" w l,\
     "waveguide_straight.fout" u 1:5 title "Im{Z} :GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_22_WAVEGUIDE_STRAIGHT_ANIMATION_beta.jpg"
set autoscale x
set autoscale y
set title "Propagation constant"
set xlabel "Frequency (Hz)"
set ylabel "m^-1"
plot "waveguide_straight.fout" u 1:2 title "Re{Beta} :GGI_TLM" w l,\
     "waveguide_straight.fout" u 1:3 title "Im{Beta} :GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_22_WAVEGUIDE_STRAIGHT_ANIMATION_R.jpg"
set autoscale x
set autoscale y
set title "Reflection Coefficient"
set xlabel "Frequency (Hz)"
set ylabel "R"
plot "waveguide_straight.fout" u 1:6 title "Re{R} :GGI_TLM" w l,\
     "waveguide_straight.fout" u 1:7 title "Im{R} :GGI_TLM" w l
#PAUSE pause -1



