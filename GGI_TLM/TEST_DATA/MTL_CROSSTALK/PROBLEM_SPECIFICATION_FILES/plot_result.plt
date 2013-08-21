#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_52_MTL_CROSSTALK_current.jpg"
set title "Cable current"
set xlabel "Frequency (Hz) "
set ylabel "Current (A)"
plot "mtl.cable_current.tout" u 1:3 w l
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_52_MTL_CROSSTALK_magnitude.jpg"
set autoscale x
set autoscale y
set title "Near end crosstalk"
set logscale x
set logscale y
set xlabel "Frequency (Hz) "
set ylabel "Magnitude (mV)"
plot "crosstalk.fout" u 1:5 w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_52_MTL_CROSSTALK_phase.jpg"
set autoscale x
set autoscale y
set title "Near end crosstalk"
set logscale x
set nologscale y
set xlabel "Frequency (Hz) "
set ylabel "Phase (radians)"
plot "crosstalk.fout" u 1:6 w l
#PAUSE pause -1
