#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_DOCUMENTATION_TEST.jpg"
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "RCS (dBsm)"
plot "geometry.rcs.fout" u 1:5 title "GGI_TLM" w l
#PAUSE pause -1
