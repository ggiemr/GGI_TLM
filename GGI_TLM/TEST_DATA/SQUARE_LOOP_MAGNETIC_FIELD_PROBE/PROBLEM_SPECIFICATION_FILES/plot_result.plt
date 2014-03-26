#JPG set term jpeg

#OUTPUT_TO_FILE set output "Magnetic_field_probe_field.jpg"

set xlabel "Time(s)"
set ylabel "Magnetic Field (A/m)"
plot "probe.field.tout" u 1:3 title "Hz at centre of loop" w l,\
     "probe.volume_average_field.tout" u 1:3 title "Hz averaged over loop" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Magnetic_field_probe_current.jpg"

set xlabel "Time(s)"
set ylabel "Probe current (A)"
plot "probe_current.tout" u 1:3 w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Magnetic_field_probe_normalised_voltage.jpg"
set logscale x
set xlabel "Frequency (Hz)"
set ylabel "Normalised Voltage"
plot "normalised_voltage.fout" u 1:5 w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Magnetic_field_probe_normalised_voltage2.jpg"

set xlabel "Frequency (Hz)"
set ylabel "Normalised Voltage (dB)"
plot "normalised_voltage.fout" u 1:7 title "Vprobe/Hincident" w l ,\
     "normalised_voltage2.fout" u 1:7 title "Vprobe/Hzaverage" w l,\
     "PROBLEM_SPECIFICATION_FILES/Vprobe.fout_analytic" u 1:5 title "Analytic reference" w p
#PAUSE pause -1
