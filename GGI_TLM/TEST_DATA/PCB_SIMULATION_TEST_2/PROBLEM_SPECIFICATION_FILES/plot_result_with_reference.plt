#JPG set term jpeg

#OUTPUT_TO_FILE set output "Source_voltage.jpg"

plot "PROBLEM_SPECIFICATION_FILES/source_voltage1.tout_ref" u 1:3 title "Source end voltage1: GGI_TLM Reference " w p,\
     "source_voltage1.tout" u 1:3 title "source voltage1: GGI_TLM" w l ,\
     "PROBLEM_SPECIFICATION_FILES/load_voltage1.tout_ref" u 1:3 title "Source end voltage1: GGI_TLM Reference " w p,\
     "load_voltage1.tout" u 1:3 title "load voltage1: GGI_TLM" w l 
     
#PAUSE pause -1

plot "PROBLEM_SPECIFICATION_FILES/source_voltage2.tout_ref" u 1:3 title "Source end voltage2: GGI_TLM Reference " w p,\
     "source_voltage2.tout" u 1:3 title "source voltage2: GGI_TLM" w l ,\
     "PROBLEM_SPECIFICATION_FILES/load_voltage2.tout_ref" u 1:3 title "Source end voltage2: GGI_TLM Reference " w p,\
     "load_voltage2.tout" u 1:3 title "load voltage2: GGI_TLM" w l 
     
#PAUSE pause -1

