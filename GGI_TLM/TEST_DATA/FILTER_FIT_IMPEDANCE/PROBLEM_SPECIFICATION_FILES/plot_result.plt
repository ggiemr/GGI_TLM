#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_14_filter_fit_Impedance.jpg"

set xlabel "Frequency(Hz)"
set ylabel "Impedance (ohms)"
plot "test_impedance.Z.fd_input" u 1:2 title "Re{Z_test_impedance}: Input data" w p,\
     "test_impedance.Z.fd_input" u 1:3 title "Im{Z_test_impedance}: Input data" w p,\
     "test_impedance.2.Z.fd_trial" u 1:2 title "Re{Z_test_impedance}: Model order 2" w l,\
     "test_impedance.2.Z.fd_trial" u 1:3 title "Im{Z_test_impedance}: Model order 2" w l
     
#PAUSE pause -1
#OUTPUT_TO_FILE quit

set title "Filter_convert test S->S_PZ"
set autoscale x
set autoscale y

set xlabel "Frequency(Hz)"
set ylabel "Impedance (ohms)"
plot "test_impedance.2.Z.fd_trial" u 1:2 title "Re{Zr_test_impedance}: S_filter Model order 2" w p,\
     "test_impedance.2.Z.fd_trial" u 1:3 title "Im{Zr_test_impedance}: S_filter Model order 2" w p,\
     "test_impedance.2.Z.Z.fd_trial" u 1:2 title "Re{Zr_test_impedance}: S_filter Model order 2" w p,\
     "test_impedance.2.Z.Z.fd_trial" u 1:3 title "Im{Zr_test_impedance}: S_filter Model order 2" w p,\
     "test_impedance.2.Z.Z.fd_trial" u 1:4 title "Re{Zr_test_impedance}: S_PZ_filter Model order 2" w l,\
     "test_impedance.2.Z.Z.fd_trial" u 1:5 title "Im{Zr_test_impedance}: S_PZ_filter Model order 2" w l
     
pause -1

set title "Filter_convert test S->S_PR"
set autoscale x
set autoscale y

set xlabel "Frequency(Hz)"
set ylabel "Impedance (ohms)"
plot "test_impedance.2.Z.fd_trial" u 1:2 title "Re{Zr_test_impedance}: S_filter Model order 2" w p,\
     "test_impedance.2.Z.fd_trial" u 1:3 title "Im{Zr_test_impedance}: S_filter Model order 2" w p,\
     "test_impedance.2.Z.Z.fd_trial" u 1:2 title "Re{Zr_test_impedance}: S_filter Model order 2" w p,\
     "test_impedance.2.Z.Z.fd_trial" u 1:3 title "Im{Zr_test_impedance}: S_filter Model order 2" w p,\
     "test_impedance.2.Z.Z.fd_trial" u 1:6 title "Re{Zr_test_impedance}: S_PR_filter Model order 2" w l,\
     "test_impedance.2.Z.Z.fd_trial" u 1:7 title "Im{Zr_test_impedance}: S_PR_filter Model order 2" w l
     
pause -1
