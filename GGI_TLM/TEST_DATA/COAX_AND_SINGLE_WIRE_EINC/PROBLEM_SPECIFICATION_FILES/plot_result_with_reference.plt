#JPG set term jpeg

#OUTPUT_TO_FILE set output "COAX_AND_SINGLE_WIRE_EINC_ref.jpg"

set xlabel "Frequency (Hz)"
set ylabel "Crosstalk"
plot "PROBLEM_SPECIFICATION_FILES/normalised_source_voltage.fout_ref" u 1:7 title "Coax field: GGI_TLM Reference " w p,\
     "normalised_source_voltage.fout" u 1:7 w l
#PAUSE pause -1
