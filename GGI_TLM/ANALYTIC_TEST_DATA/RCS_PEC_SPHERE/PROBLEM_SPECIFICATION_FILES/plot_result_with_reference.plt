#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_1_RCS_sphere_rcs_ref.jpg"
set xlabel "Angle (degrees)"
set ylabel "RCS (dBsm)"
plot "pec.rcs" u 2:4 title "RCS: GGI_RCS_sphere" w l,\
     "PROBLEM_SPECIFICATION_FILES/pec.rcs_ref" u 2:4 title "RCS: GGI_RCS_sphere reference" w p
     
#PAUSE pause -1
