#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_43_COAX_COUPLING_crosstalk_ref.jpg"

set xlabel "Frequency (Hz)"
set ylabel "Crosstalk"
plot "PROBLEM_SPECIFICATION_FILES/crosstalk.fout_ref" u 1:7 title "Coax field: Fieldsolve Reference " w p,\
     "crosstalk.fout" u 1:7 w l
#PAUSE pause -1
