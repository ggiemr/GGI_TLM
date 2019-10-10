#JPG set term jpeg

#OUTPUT_TO_FILE set output "parallel_conductors_current.jpg"
set xlabel "Time (s)"
set ylabel "A"
plot "PROBLEM_SPECIFICATION_FILES/parallel_conductors.cable_current.tout_ref" u 1:3 title "Port 1 current reference" w p ,\
     "PROBLEM_SPECIFICATION_FILES/current_from_H.tout_ref" u 1:3 title "Conductor current from H fieldreference" w p,\
     "parallel_conductors.cable_current.tout" u 1:3 title "Port 1 current" w l ,\
     "current_from_H.tout" u 1:3 title "Conductor current from H field" w l,\
     
#PAUSE pause -1

#OUTPUT_TO_FILE set output "parallel_conductors_field.jpg"
set xlabel "Time (s)"
set ylabel "A/m"
set autoscale x
set autoscale y
plot "PROBLEM_SPECIFICATION_FILES/Hy.tout_ref" u 1:3 title "Hy reference" w p ,\
     "Hy.tout" u 1:3 title "Hy" w l
#PAUSE pause -1

#OUTPUT_TO_FILE set output "parallel_conductors_Force.jpg"
set xlabel "Time (s)"
set ylabel "Force per unit length (N/m)"
set autoscale x
set autoscale y
plot "PROBLEM_SPECIFICATION_FILES/Fx.tout_ref" u 1:3 title "Fx reference" w p ,\
     "Fx.tout" u 1:3 title "Fx" w l 
#PAUSE pause -1


