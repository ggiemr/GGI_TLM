#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_24_WIRE_OVER_LOSSY_GROUND_current_ref.jpg"
set xlabel "Time (s)"
set ylabel "A"
plot "PROBLEM_SPECIFICATION_FILES/Current_1.tout_ref" u 1:3 title "Cable current 1: GGI_TLM Reference" w p ,\
     "PROBLEM_SPECIFICATION_FILES/Current_2.tout_ref" u 1:3 title "Cable current 2: GGI_TLM Reference" w p ,\
     "PROBLEM_SPECIFICATION_FILES/Current_3.tout_ref" u 1:3 title "Cable current 3: GGI_TLM Reference" w p ,\
     "PROBLEM_SPECIFICATION_FILES/Current_4.tout_ref" u 1:3 title "Cable current 4: GGI_TLM Reference" w p ,\
     "Current_1.tout" u 1:3 title "Cable current 1: GGI_TLM" w l,\
     "Current_2.tout" u 1:3 title "Cable current 2: GGI_TLM" w l,\
     "Current_3.tout" u 1:3 title "Cable current 3: GGI_TLM" w l,\
     "Current_4.tout" u 1:3 title "Cable current 4: GGI_TLM" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_24_WIRE_OVER_LOSSY_GROUND_velocity_ref.jpg"
set logscale x
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "Wave velocity"
plot "PROBLEM_SPECIFICATION_FILES/wire_over_ground.fout_ref" u 1:10 title "Re{w/beta}: GGI_TLM Reference" w p ,\
     "wire_over_ground.fout" u 1:10 title " Re{w/beta}" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_24_WIRE_OVER_LOSSY_GROUND_beta_ref.jpg"
set logscale x
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "Imaginary part of Propagation constant"
plot "PROBLEM_SPECIFICATION_FILES/wire_over_ground.fout_ref" u 1:3 title "Re{beta}: GGI_TLM Reference" w p ,\
     "wire_over_ground.fout" u 1:3 title " Im{beta}" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_24_WIRE_OVER_LOSSY_GROUND_Z_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "Impedance"
plot "PROBLEM_SPECIFICATION_FILES/wire_over_ground.fout_ref" u 1:4 title "Re{Z}: GGI_TLM Reference" w p ,\
     "PROBLEM_SPECIFICATION_FILES/wire_over_ground.fout_ref" u 1:5 title "Im{Z}: GGI_TLM Reference" w p ,\
     "wire_over_ground.fout" u 1:4 title " Re{Z}" w l,\
     "wire_over_ground.fout" u 1:5 title " Im{Z}" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_24_WIRE_OVER_LOSSY_GROUND_R_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "R"
plot "PROBLEM_SPECIFICATION_FILES/wire_over_ground.fout_ref" u 1:6 title "Re{R}: GGI_TLM Reference" w p ,\
     "PROBLEM_SPECIFICATION_FILES/wire_over_ground.fout_ref" u 1:7 title "Im{R}: GGI_TLM Reference" w p ,\
     "wire_over_ground.fout" u 1:6 title " Re{R}" w l,\
     "wire_over_ground.fout" u 1:7 title " Im{R}" w l 
#PAUSE pause -1
