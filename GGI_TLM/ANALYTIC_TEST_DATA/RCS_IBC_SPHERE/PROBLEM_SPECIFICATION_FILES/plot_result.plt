#JPG set term jpeg
#OUTPUT_TO_FILE set output "RCS_IBC_SPHERE.jpg"

set xlabel "Frequency(Hz)"
set ylabel "Surface Impedance (Ohms/square)"
plot "z11.fout" u 1:2 title "Re{Zsibc}" w l,\
     "z11.fout" u 1:3 title "Im{Zsibc}" w l
     
pause -1

#OUTPUT_TO_FILE set output "Test_case_1_RCS_sphere_rcs.jpg"
set xlabel "Angle (degrees)"
set ylabel "RCS (dBsm)"
plot "sphere.rcs" u 2:4 title "RCS: GGI_RCS_sphere" w l
     
pause -1
