#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_1_RCS_sphere_epsr.jpg"
set xlabel "Frequency(Hz)"
set ylabel "Relative permittivity"
plot "epsr.fout" u 1:2 title "Re{epsr}" w l,\
     "epsr.fout" u 1:3 title "Im{epsr}" w l
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_1_RCS_sphere_mur.jpg"
set xlabel "Frequency(Hz)"
set ylabel "Relative permeability"
plot "mur.fout" u 1:2 title "Re{epsr}" w l,\
     "mur.fout" u 1:3 title "Im{epsr}" w l
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_1_RCS_sphere_rcs.jpg"
set xlabel "Angle (degrees)"
set ylabel "RCS (dBsm)"
plot "pec.rcs" u 2:4 title "RCS: GGI_RCS_sphere" w l
     
#PAUSE pause -1
