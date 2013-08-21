#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_52_MTL_CROSSTALK_magnitude_ref.jpg"

set title "Near end crosstalk"
set logscale x
set logscale y
set xlabel "Frequency (Hz) "
set ylabel "Magnitude (mV)"
plot "PROBLEM_SPECIFICATION_FILES/crosstalk.fout_ref" u 1:5 title "Crosstalk: GGI_TLM Reference " w p,\
     "crosstalk.fout" u 1:5 w l
#PAUSE pause -1


#OUTPUT_TO_FILE set output "Test_case_52_MTL_CROSSTALK_phase.jpg"
set autoscale x
set autoscale y
set title "Near end crosstalk"
set logscale x
set nologscale y
set xlabel "Frequency (Hz) "
set ylabel "Phase (radians)"
plot "PROBLEM_SPECIFICATION_FILES/crosstalk.fout_ref" u 1:6 title "Crosstalk: GGI_TLM Reference " w p,\
     "crosstalk.fout" u 1:6 w l
#PAUSE pause -1
