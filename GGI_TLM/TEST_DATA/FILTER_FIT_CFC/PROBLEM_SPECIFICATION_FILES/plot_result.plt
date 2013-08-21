#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_13_filter_fit_CFC_Z11.jpg"

set title "CFC, z11"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"
plot "CFC.z11.fd_input" u 1:2 title "Re{z11_CFC}" w p,\
     "CFC.z11.fd_input" u 1:3 title "Im{z11_CFC}" w p,\
     "CFC.2.z11.fd_trial" u 1:2 title "Re{z11_CFC}, model order 2" w l,\
     "CFC.2.z11.fd_trial" u 1:3 title "Im{z11_CFC}, model order 2" w l
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_13_filter_fit_CFC_Z12.jpg"

set title "CFC, z12"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"
plot "CFC.z12.fd_input" u 1:2 title "Re{z12_CFC}" w p,\
     "CFC.z12.fd_input" u 1:3 title "Im{z12_CFC}" w p,\
     "CFC.2.z12.fd_trial" u 1:2 title "Re{z12_CFC}, model order 2" w l,\
     "CFC.2.z12.fd_trial" u 1:3 title "Im{z12_CFC}, model order 2" w l
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_13_filter_fit_CFC_Z22.jpg"

set title "CFC, z22"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"
plot "CFC.z22.fd_input" u 1:2 title "Re{z22_CFC}" w p,\
     "CFC.z22.fd_input" u 1:3 title "Im{z22_CFC}" w p,\
     "CFC.2.z22.fd_trial" u 1:2 title "Re{z22_CFC}, model order 2" w l,\
     "CFC.2.z22.fd_trial" u 1:3 title "Im{z22_CFC}, model order 2" w l
     
#PAUSE pause -1
#OUTPUT_TO_FILE quit

set title "CFC, z11 Filter_convert test S->S_PZ"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"
plot "CFC.2.z11.fd_trial" u 1:2 title "Re{z11_CFC}, model order 2" w p,\
     "CFC.2.z11.fd_trial" u 1:3 title "Im{z11_CFC}, model order 2" w p,\
     "CFC.2.smat.z11.fd_trial" u 1:2 title "Re{z11_CFC}, S_filter model order 2" w p,\
     "CFC.2.smat.z11.fd_trial" u 1:3 title "Re{z11_CFC}, S_filter model order 2" w p,\
     "CFC.2.smat.z11.fd_trial" u 1:4 title "Re{z11_CFC}, S_PZ_filter model order 2" w l,\
     "CFC.2.smat.z11.fd_trial" u 1:5 title "Re{z11_CFC}, S_PZ_filter model order 2" w l
pause -1

set title "CFC, z12 Filter_convert test S->S_PZ"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"
plot "CFC.2.z12.fd_trial" u 1:2 title "Re{z12_CFC}, model order 2" w p,\
     "CFC.2.z12.fd_trial" u 1:3 title "Im{z12_CFC}, model order 2" w p,\
     "CFC.2.smat.z12.fd_trial" u 1:2 title "Re{z12_CFC}, S_filter model order 2" w p,\
     "CFC.2.smat.z12.fd_trial" u 1:3 title "Re{z12_CFC}, S_filter model order 2" w p,\
     "CFC.2.smat.z12.fd_trial" u 1:4 title "Re{z12_CFC}, S_PZ_filter model order 2" w l,\
     "CFC.2.smat.z12.fd_trial" u 1:5 title "Re{z12_CFC}, S_PZ_filter model order 2" w l
pause -1

set title "CFC, z22 Filter_convert test S->S_PZ"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"
plot "CFC.2.z22.fd_trial" u 1:2 title "Re{z22_CFC}, model order 2" w p,\
     "CFC.2.z22.fd_trial" u 1:3 title "Im{z22_CFC}, model order 2" w p,\
     "CFC.2.smat.z22.fd_trial" u 1:2 title "Re{z22_CFC}, S_filter model order 2" w p,\
     "CFC.2.smat.z22.fd_trial" u 1:3 title "Re{z22_CFC}, S_filter model order 2" w p,\
     "CFC.2.smat.z22.fd_trial" u 1:4 title "Re{z22_CFC}, S_PZ_filter model order 2" w l,\
     "CFC.2.smat.z22.fd_trial" u 1:5 title "Re{z22_CFC}, S_PZ_filter model order 2" w l
pause -1

set title "CFC, z11 Filter_convert test S->S_PR"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"
plot "CFC.2.z11.fd_trial" u 1:2 title "Re{z11_CFC}, model order 2" w p,\
     "CFC.2.z11.fd_trial" u 1:3 title "Im{z11_CFC}, model order 2" w p,\
     "CFC.2.smat.z11.fd_trial" u 1:2 title "Re{z11_CFC}, S_filter model order 2" w p,\
     "CFC.2.smat.z11.fd_trial" u 1:3 title "Re{z11_CFC}, S_filter model order 2" w p,\
     "CFC.2.smat.z11.fd_trial" u 1:6 title "Re{z11_CFC}, S_PR_filter model order 2" w l,\
     "CFC.2.smat.z11.fd_trial" u 1:7 title "Re{z11_CFC}, S_PR_filter model order 2" w l
pause -1

set title "CFC, z12 Filter_convert test S->S_PR"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"
plot "CFC.2.z12.fd_trial" u 1:2 title "Re{z12_CFC}, model order 2" w p,\
     "CFC.2.z12.fd_trial" u 1:3 title "Im{z12_CFC}, model order 2" w p,\
     "CFC.2.smat.z12.fd_trial" u 1:2 title "Re{z12_CFC}, S_filter model order 2" w p,\
     "CFC.2.smat.z12.fd_trial" u 1:3 title "Re{z12_CFC}, S_filter model order 2" w p,\
     "CFC.2.smat.z12.fd_trial" u 1:6 title "Re{z12_CFC}, S_PR_filter model order 2" w l,\
     "CFC.2.smat.z12.fd_trial" u 1:7 title "Re{z12_CFC}, S_PR_filter model order 2" w l
pause -1

set title "CFC, z22 Filter_convert test S->S_PR"
set xlabel "Frequency(Hz)"
set ylabel "Ohms"
plot "CFC.2.z22.fd_trial" u 1:2 title "Re{z22_CFC}, model order 2" w p,\
     "CFC.2.z22.fd_trial" u 1:3 title "Im{z22_CFC}, model order 2" w p,\
     "CFC.2.smat.z22.fd_trial" u 1:2 title "Re{z22_CFC}, S_filter model order 2" w p,\
     "CFC.2.smat.z22.fd_trial" u 1:3 title "Re{z22_CFC}, S_filter model order 2" w p,\
     "CFC.2.smat.z22.fd_trial" u 1:6 title "Re{z22_CFC}, S_PR_filter model order 2" w l,\
     "CFC.2.smat.z22.fd_trial" u 1:7 title "Re{z22_CFC}, S_PR_filter model order 2" w l
pause -1
