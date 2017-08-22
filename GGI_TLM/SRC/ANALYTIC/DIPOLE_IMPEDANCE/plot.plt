set yrange[-1000:1000]
set title "Dipole self impedance"
set xlabel "Frequency(Hz)"
set ylabel "Impedance (ohms)"
plot "Zs.dat" u 1:2 title "Re{Z}" w l, "Zs.dat" u 1:3 title "Im{Z}" w l
pause -1

set title "Dipole mutual impedance"
set autoscale x
set autoscale y
set xlabel "Separation (m)"
set ylabel "Impedance (ohms)"
plot "Zm.dat" u 1:2 title "Re{Z}" w l, "Zm.dat" u 1:3 title "Im{Z}" w l
pause -1


