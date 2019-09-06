#JPG set term jpeg

#OUTPUT_TO_FILE set output "Source_voltage.jpg"

plot "PROBLEM_SPECIFICATION_FILES/v1.tout_ref" u 1:3 title "v1: GGI_TLM Reference " w p,\
                                 "v1.tout" u 1:3     title "v1: GGI_TLM" w l ,\
     "PROBLEM_SPECIFICATION_FILES/v2.tout_ref" u 1:3 title "v2: GGI_TLM Reference " w p,\
                                 "v2.tout" u 1:3     title "v2: GGI_TLM" w l ,\
     "PROBLEM_SPECIFICATION_FILES/v3.tout_ref" u 1:3 title "v3: GGI_TLM Reference " w p,\
                                 "v3.tout" u 1:3     title "v3: GGI_TLM" w l ,\
     "PROBLEM_SPECIFICATION_FILES/v4.tout_ref" u 1:3 title "v4: GGI_TLM Reference " w p,\
                                 "v4.tout" u 1:3     title "v4: GGI_TLM" w l 
     
#PAUSE pause -1
