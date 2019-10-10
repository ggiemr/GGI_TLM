#JPG set term jpeg

#OUTPUT_TO_FILE set output "parallel_conductors_current.jpg"
set xlabel "Time (s)"
set ylabel "A"
plot "parallel_conductors.cable_current.tout" u 1:3 title "Port 1 current" w l ,\
     "current_from_H.tout" u 1:3 title "Conductor current from H field" w l ,\
     "current_from_sigma_E.tout" u 1:3 title "Conductor current from integral of sigma* E field" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "parallel_conductors_field.jpg"
set xlabel "Time (s)"
set ylabel "A/m"
set autoscale x
set autoscale y
plot "Hy.tout" u 1:3 title "Hy average" w l,\
     "Hy_wire_centre.tout" u 1:3 title "Hy centre" w l
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "parallel_conductors_Force.jpg"
set xlabel "Time (s)"
set ylabel "Force per unit length (N/m)"
set autoscale x
set autoscale y
plot "Fx.tout" u 1:3 title "Fx" w l ,\
     "Fx_wire_centre.tout" u 1:3 title "Fx" w l 
     
#PAUSE pause -1


