#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_43_COAX_COUPLING_crosstalk.jpg"

set xlabel "Frequency (Hz)"
set ylabel "Crosstalk"
plot "crosstalk.fout" u 1:7 w l
#PAUSE pause -1
