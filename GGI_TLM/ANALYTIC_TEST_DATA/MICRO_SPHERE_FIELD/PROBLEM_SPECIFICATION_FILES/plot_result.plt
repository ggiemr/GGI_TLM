#JPG set term jpeg

##OUTPUT_TO_FILE set output "Test_case_3_MICRO_SPHERE_FIELD_epsr.jpg"
#set xlabel "Frequency(Hz)"
#set ylabel "Relative permittivity"
#plot "epsr.fout" u 1:2 title "Re{epsr}" w l,\
#     "epsr.fout" u 1:3 title "Im{epsr}" w l
#     
##PAUSE pause -1
#
##OUTPUT_TO_FILE set output "Test_case_3_MICRO_SPHERE_FIELD_mur.jpg"
#set xlabel "Frequency(Hz)"
#set ylabel "Relative permeability"
#plot "mur.fout" u 1:2 title "Re{epsr}" w l,\
#     "mur.fout" u 1:3 title "Im{epsr}" w l
#     
##PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_3_MICRO_SPHERE_FIELD_rcs.jpg"
set xlabel "Time (s)"
set ylabel "Electric Field (V/m)"
plot "micro_sphere.rcs" u 1:4 title "E field: GGI_RCS_sphere" w l
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "Test_case_3_MICRO_SPHERE_FIELD_field.jpg"
set xlabel "Time (s)"
set ylabel "Electric Field (V/m)"
plot "micro_sphere.field" u 1:2 title "E field: GGI_RCS_sphere" w l
     
#PAUSE pause -1
