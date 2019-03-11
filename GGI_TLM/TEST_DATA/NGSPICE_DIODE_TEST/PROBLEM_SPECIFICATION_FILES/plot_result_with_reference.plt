#JPG set term jpeg

#OUTPUT_TO_FILE set output "Source_voltage.jpg"

plot "PROBLEM_SPECIFICATION_FILES/source_voltage.tout_ref" u 1:3 title "Source end voltage: GGI_TLM Reference " w p,\
     "source_voltage.tout" u 1:3 title "source voltage: GGI_TLM" w l 
     
#PAUSE pause -1

