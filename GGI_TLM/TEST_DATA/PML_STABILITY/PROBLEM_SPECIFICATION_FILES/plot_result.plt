#JPG set term jpeg

set autoscale x
set autoscale y

set xlabel "Time (s)"
set ylabel "V/m"

#OUTPUT_TO_FILE set output "PML_stability.jpg"

plot "cavity.field.tout" u 1:3 title "E field: GGI_TLM" w l 
     
#PAUSE pause -1
