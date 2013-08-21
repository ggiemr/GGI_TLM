#JPG set term jpeg

set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_54_RESONANT_ABSORBER_field.jpg"

set xlabel "Time (s)"
set ylabel "V/m"
plot "Ey.tout" u 1:3 title "Ey: GGI_TLM" w l 
#PAUSE pause -1


set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_54_RESONANT_ABSORBER_beta.jpg"

set xlabel "Frequency (Hz)"
set ylabel "Propagetion constant (m^-1)"
plot "Resonant_absorber.fout" u 1:2 title "Re{beta}, GGI_TLM" w l ,\
     "Resonant_absorber.fout" u 1:3 title "Im{beta}, GGI_TLM" w l 
#PAUSE pause -1

set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_54_RESONANT_ABSORBER_Z0.jpg"

set xlabel "Frequency (Hz)"
set ylabel "Wave impedance (ohms)"
plot "Resonant_absorber.fout" u 1:4 title "Re{Z}, GGI_TLM" w l ,\
     "Resonant_absorber.fout" u 1:5 title "Im{Z}, GGI_TLM" w l 
#PAUSE pause -1

set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "Reflection coefficient "

#OUTPUT_TO_FILE set output "Test_case_54_RESONANT_ABSORBER_R.jpg"

plot "Resonant_absorber.fout" u 1:6 title "Re{R}, GGI_TLM" w l ,\
     "Resonant_absorber.fout" u 1:7 title "Im{R}, GGI_TLM" w l  
#PAUSE pause -1

set xrange [5e9:16e9]
set yrange [0:1]

#OUTPUT_TO_FILE set output "Test_case_54_RESONANT_ABSORBER_R2.jpg"
plot  "Resonant_absorber.fout" u 1:8 title "  |R|, GGI_TLM" w l 
#PAUSE pause -1
