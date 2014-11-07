plot "monopole_original_TLM.cable_current.tout" u 1:3 w p,\
     "../monopole.cable_current.tout" u 1:3 w p
pause -1

plot "monopole_original_TLM.field.tout" u 1:3 w l,\
     "../monopole.field.tout" u 1:3 w l
pause -1
