set term X11 enhanced
set termoption enhanced
set xlabel "time(s)"
set ylabel "V"

#set autoscale x
#set autoscale y
#plot "ce336_600w.dat" u 1:2 title "600W" w l
#pause -1

set xlabel "Frequency(Hz)"
set ylabel "V"

set autoscale x
set autoscale y
set logscale x
plot "ce336_600w.favg" u 1:6 title "600W" w l
pause -1

quit

