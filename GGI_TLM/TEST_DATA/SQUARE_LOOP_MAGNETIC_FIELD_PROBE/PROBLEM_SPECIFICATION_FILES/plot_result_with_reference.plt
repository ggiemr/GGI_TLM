#JPG set term jpeg

#OUTPUT_TO_FILE set output "Magnetic_field_probe_normalised_voltage_ref.jpg"

set logscale x
set xlabel "Frequency (Hz)"
set ylabel "Normalised voltage"
plot "PROBLEM_SPECIFICATION_FILES/normalised_voltage.fout_ref" u 1:7 title "Normalised voltage: GGI_TLM Reference " w p,\
     "normalised_voltage.fout" u 1:7 title "Normalised voltage: GGI_TLM" w l
#PAUSE pause -1
