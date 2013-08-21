#JPG set term jpeg

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V/m"

#OUTPUT_TO_FILE set output "Test_case_54_RESONANT_ABSORBER_field_ref.jpg"

plot "PROBLEM_SPECIFICATION_FILES/Ey.tout_ref" u 1:3 title "Ey: GGI_TLM reference" w p,\
     "Ey.tout" u 1:3 title "Ey: GGI_TLM" w l
#PAUSE pause -1


set autoscale x
set autoscale y

#OUTPUT_TO_FILE set output "Test_case_54_RESONANT_ABSORBER_R_ref.jpg"
set xrange [5e9:16e9]
set yrange [0:1]
set xlabel "Frequency (Hz)"
set ylabel "Reflection coefficient "
plot "PROBLEM_SPECIFICATION_FILES/Resonant_absorber.fout_ref" u 1:8 title "  |R|, GGI_TLM reference" w p ,\
     "Resonant_absorber.fout" u 1:8 title "  |R|, GGI_TLM" w l  
#PAUSE pause -1
