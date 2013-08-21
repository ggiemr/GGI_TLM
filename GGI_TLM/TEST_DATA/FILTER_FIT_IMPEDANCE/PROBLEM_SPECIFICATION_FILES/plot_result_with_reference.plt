#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_14_filter_fit_Impedance_ref.jpg"

set xlabel "Frequency(Hz)"
set ylabel "Impedance (ohms)"
plot "test_impedance.Z.fd_input" u 1:2 title "Re{Z_test_impedance}: Input data" w p,\
     "test_impedance.Z.fd_input" u 1:3 title "Im{Z_test_impedance}: Input data" w p,\
     "PROBLEM_SPECIFICATION_FILES/test_impedance.2.Z.fd_trial_ref" u 1:2 title "Re{Z_test_impedance}: Reference model order 2" w p,\
     "PROBLEM_SPECIFICATION_FILES/test_impedance.2.Z.fd_trial_ref" u 1:3 title "Im{Z_test_impedance}: Reference model order 2" w p,\
     "test_impedance.2.Z.fd_trial" u 1:2 title "Re{Z_test_impedance}: Model order 2" w l,\
     "test_impedance.2.Z.fd_trial" u 1:3 title "Im{Z_test_impedance}: Model order 2" w l
     
#PAUSE pause -1
