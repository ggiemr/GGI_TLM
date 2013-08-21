#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_15_filter_fit_ls22_ref.jpg"

set xlabel "Frequency(Hz)"
set ylabel "Relative permittivity"
plot "ls22.eps.fd_input" u 1:2 title "Re{epsr_ls22}: Input data" w p,\
     "ls22.eps.fd_input" u 1:3 title "Im{epsr_ls22}: Input data" w p,\
     "PROBLEM_SPECIFICATION_FILES/ls22.2.eps.fd_trial_ref" u 1:2 title "Re{epsr_ls22}: Reference model order 2" w p,\
     "PROBLEM_SPECIFICATION_FILES/ls22.2.eps.fd_trial_ref" u 1:3 title "Im{epsr_ls22}: Reference model order 2" w p,\
     "ls22.2.eps.fd_trial" u 1:2 title "Re{epsr_ls22}: Model order 2" w l,\
     "ls22.2.eps.fd_trial" u 1:3 title "Im{epsr_ls22}: Model order 2" w l
     
#PAUSE pause -1
