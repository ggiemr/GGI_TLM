#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_5_RCS_coated_IBC_sphere_freq_epsr.jpg"
set xlabel "Frequency(Hz)"
set ylabel "Relative permittivity"
plot "epsr.fout" u 1:2 title "Re{epsr}" w l,\
     "epsr.fout" u 1:3 title "Im{epsr}" w l
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_5_RCS_coated_IBC_sphere_freq_mur.jpg"
set xlabel "Frequency(Hz)"
set ylabel "Relative permeability"
plot "mur.fout" u 1:2 title "Re{epsr}" w l,\
     "mur.fout" u 1:3 title "Im{epsr}" w l
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_5_RCS_coated_IBC_sphere_freq_rcs.jpg"
set xlabel "Frequency(Hz)"
set ylabel "RCS (dBsm)"
plot "coated_ibc.rcs" u 1:4 title "RCS: GGI_RCS_sphere" w l
     
#PAUSE pause -1
