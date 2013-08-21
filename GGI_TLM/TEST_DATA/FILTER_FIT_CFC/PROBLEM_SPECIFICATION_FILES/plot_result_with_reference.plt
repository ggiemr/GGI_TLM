#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_13_filter_fit_CFC_Z11_ref.jpg"

set title "CFC, Z11"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"
plot "CFC.z11.fd_input" u 1:2 title "Re{Z11_CFC}" w p,\
     "CFC.z11.fd_input" u 1:3 title "Im{Z11_CFC}" w p,\
     "PROBLEM_SPECIFICATION_FILES/CFC.2.z11.fd_trial_ref" u 1:2 title "Re{z11_CFC}, reference model order 2" w p,\
     "PROBLEM_SPECIFICATION_FILES/CFC.2.z11.fd_trial_ref" u 1:3 title "Im{z11_CFC}, reference model order 2" w p,\
     "CFC.2.z11.fd_trial" u 1:2 title "Re{Z11_CFC}, model order 2" w l,\
     "CFC.2.z11.fd_trial" u 1:3 title "Im{Z11_CFC}, model order 2" w l
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_13_filter_fit_CFC_Z12_ref.jpg"

set title "CFC, z12"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"
plot "CFC.z12.fd_input" u 1:2 title "Re{z12_CFC}" w p,\
     "CFC.z12.fd_input" u 1:3 title "Im{z12_CFC}" w p,\
     "PROBLEM_SPECIFICATION_FILES/CFC.2.z12.fd_trial_ref" u 1:2 title "Re{z12_CFC}, reference model order 2" w p,\
     "PROBLEM_SPECIFICATION_FILES/CFC.2.z12.fd_trial_ref" u 1:3 title "Im{z12_CFC}, reference model order 2" w p,\
     "CFC.2.z12.fd_trial" u 1:2 title "Re{z12_CFC}, model order 2" w l,\
     "CFC.2.z12.fd_trial" u 1:3 title "Im{z12_CFC}, model order 2" w l
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_13_filter_fit_CFC_Z22_ref.jpg"

set title "CFC, z22"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"
plot "CFC.z22.fd_input" u 1:2 title "Re{z22_CFC}" w p,\
     "CFC.z22.fd_input" u 1:3 title "Im{z22_CFC}" w p,\
     "PROBLEM_SPECIFICATION_FILES/CFC.2.z22.fd_trial_ref" u 1:2 title "Re{z22_CFC}, reference model order 2" w p,\
     "PROBLEM_SPECIFICATION_FILES/CFC.2.z22.fd_trial_ref" u 1:3 title "Im{z22_CFC}, reference model order 2" w p,\
     "CFC.2.z22.fd_trial" u 1:2 title "Re{z22_CFC}, model order 2" w l,\
     "CFC.2.z22.fd_trial" u 1:3 title "Im{z22_CFC}, model order 2" w l
     
#PAUSE pause -1
