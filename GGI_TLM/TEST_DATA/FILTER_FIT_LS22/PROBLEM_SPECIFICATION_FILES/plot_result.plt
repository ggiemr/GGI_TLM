#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_15_filter_fit_ls22.jpg"

set xlabel "Frequency(Hz)"
set ylabel "Relative permittivity"
plot "ls22.eps.fd_input" u 1:2 title "Re{epsr_ls22}: Input data" w p,\
     "ls22.eps.fd_input" u 1:3 title "Im{epsr_ls22}: Input data" w p,\
     "ls22.2.eps.fd_trial" u 1:2 title "Re{epsr_ls22}: Model order 2" w l,\
     "ls22.2.eps.fd_trial" u 1:3 title "Im{epsr_ls22}: Model order 2" w l
     
#PAUSE pause -1
#OUTPUT_TO_FILE quit

set title "Filter_convert test S->S_PZ"
set autoscale x
set autoscale y

set xlabel "Frequency(Hz)"
set ylabel "Relative permittivity"
plot "ls22.2.eps.fd_trial" u 1:2 title "Re{epsr_ls22}: S_filter Model order 2" w p,\
     "ls22.2.eps.fd_trial" u 1:3 title "Im{epsr_ls22}: S_filter Model order 2" w p,\
     "ls22.2.vmat.eps.fd_trial" u 1:2 title "Re{epsr_ls22}: S_filter Model order 2" w p,\
     "ls22.2.vmat.eps.fd_trial" u 1:3 title "Im{epsr_ls22}: S_filter Model order 2" w p,\
     "ls22.2.vmat.eps.fd_trial" u 1:4 title "Re{epsr_ls22}: S_PZ_filter Model order 2" w l,\
     "ls22.2.vmat.eps.fd_trial" u 1:5 title "Im{epsr_ls22}: S_PZ_filter Model order 2" w l
     
pause -1

set title "Filter_convert test S->S_PR"
set autoscale x
set autoscale y

set xlabel "Frequency(Hz)"
set ylabel "Relative permittivity"
plot "ls22.2.eps.fd_trial" u 1:2 title "Re{epsr_ls22}: S_filter Model order 2" w p,\
     "ls22.2.eps.fd_trial" u 1:3 title "Im{epsr_ls22}: S_filter Model order 2" w p,\
     "ls22.2.vmat.eps.fd_trial" u 1:2 title "Re{epsr_ls22}: S_filter Model order 2" w p,\
     "ls22.2.vmat.eps.fd_trial" u 1:3 title "Im{epsr_ls22}: S_filter Model order 2" w p,\
     "ls22.2.vmat.eps.fd_trial" u 1:6 title "Re{epsr_ls22}: S_PR_filter Model order 2" w l,\
     "ls22.2.vmat.eps.fd_trial" u 1:7 title "Im{epsr_ls22}: S_PR_filter Model order 2" w l
     
pause -1
