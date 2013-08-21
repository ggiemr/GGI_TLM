#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_25_dielectric_filter_S11_ref.jpg"

set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "S11 (dB)"
plot "PROBLEM_SPECIFICATION_FILES/filter_S11.fout_ref" u 1:7 title "S11: GGI_TLM reference" w p,\
     "filter_S11.fout" u 1:7 title "S11: GGI_TLM" w l 
    
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_25_dielectric_filter_S21_ref.jpg"

set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "S21 (dB)"
plot "PROBLEM_SPECIFICATION_FILES/filter_S21.fout_ref" u 1:7 title "S21: GGI_TLM reference" w p,\
     "filter_S21.fout" u 1:7 title "S21: GGI_TLM" w l 
    
#PAUSE pause -1

