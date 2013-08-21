#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_6_RCS_dielectric_IBC_sphere_rcs_freq_ref.jpg"
set xlabel "Frequency (Hz)"
set ylabel "RCS (dBsm)"
plot "dielectric_ibc.rcs" u 1:4 title "RCS: GGI_RCS_sphere" w l,\
     "PROBLEM_SPECIFICATION_FILES/dielectric_ibc.rcs_ref" u 1:4 title "RCS: GGI_RCS_sphere reference" w p
     
#PAUSE pause -1
