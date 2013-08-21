#set term jpeg
#set output "TCS_circular_aperture.jpg"
set xrange [1e9:1e10]
set yrange [-80:-10]
set xlabel "Frequency (Hz)"
set ylabel "Transmission Cross Section (dB square metres)"
plot "TCS_freq_stir_2.5mm.fout" u 1:7 title "TLM" w l,\
"PROBLEM_SPECIFICATION_FILES/TCS_analytic.fout" u 1:2 title "Analytic formula" w l
pause -1
