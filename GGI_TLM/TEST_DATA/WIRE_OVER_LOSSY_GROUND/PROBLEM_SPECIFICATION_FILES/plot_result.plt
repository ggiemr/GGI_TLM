#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_24_WIRE_OVER_LOSSY_GROUND_current.jpg"
set xlabel "Time (s)"
set ylabel "A"
plot "Current_1.tout" u 1:3 w l,\
     "Current_2.tout" u 1:3 w l,\
     "Current_3.tout" u 1:3 w l,\
     "Current_4.tout" u 1:3 w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_24_WIRE_OVER_LOSSY_GROUND_velocity.jpg"
set logscale x
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "Wave velocity"
plot "wire_over_ground.fout" u 1:10 title " Re{w/beta}" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_24_WIRE_OVER_LOSSY_GROUND_beta.jpg"
set logscale x
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "Imaginary part of Propagation constant"
plot "wire_over_ground.fout" u 1:3 title " Im{beta}" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_24_WIRE_OVER_LOSSY_GROUND_Z.jpg"
set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "Impedance"
plot "wire_over_ground.fout" u 1:4 title " Re{Z}" w l,\
     "wire_over_ground.fout" u 1:5 title " Im{Z}" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_24_WIRE_OVER_LOSSY_GROUND_R.jpg"
set xlabel "Frequency (Hz)"
set ylabel "R"
plot "wire_over_ground.fout" u 1:6 title " Re{R}" w l,\
     "wire_over_ground.fout" u 1:7 title " Im{R}" w l    
#PAUSE pause -1
