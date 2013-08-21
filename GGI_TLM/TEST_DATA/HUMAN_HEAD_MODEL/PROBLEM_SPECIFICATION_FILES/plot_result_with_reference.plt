#JPG set term jpeg

#OUTPUT_TO_FILE set output "Test_case_59_HEAD_field_ref.jpg"
set autoscale x
set autoscale y
set xlabel "Time (s)"
set ylabel "E field"
plot "head.field.tout" u 1:3 title "GGI_TLM" w l,\
     "PROBLEM_SPECIFICATION_FILES/head.field.tout_ref" u 1:3 title "GGI_TLM reference" w p 
#PAUSE pause -1

