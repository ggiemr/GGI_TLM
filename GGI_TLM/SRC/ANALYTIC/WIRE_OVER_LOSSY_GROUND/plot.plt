
set grid

set key top left

set title "Attenuation constant"
set xlabel 'Frequency (Hz)'
set ylabel 'Np/m'
set logscale x
#set xrange[1e4:1e7]
#set yrange[0:4]
plot 'Full_Carson_wire_over_lossy_ground.out' u 1:2 title 'Re{gamma} (Full Carson)' w l,\
     'Approx_Carson_wire_over_lossy_ground.out' u 1:2 title 'Re{gamma} (Approx Carson)' w l,\
     'Wait_wire_over_lossy_ground_all_modes.out' u 1:2 title 'Re{gamman} (Wait)' w p,\
     'Wait_wire_over_lossy_ground.out' u 1:2 title 'Re{gamma1} (Wait)' w p,\
     'Wait_wire_over_lossy_ground.out' u 1:5 title 'Re{gamma2} (Wait)' w p
     
pause -1

set title "Propagation constant"
set nologscale x
set ylabel 'rad/m'
set autoscale x
set autoscale y
set yrange[0:]
plot 'k0.out' u 1:2 title 'k0' w l,\
     'Full_Carson_wire_over_lossy_ground.out' u 1:3 title 'Im{gamma} (Full Carson)' w l,\
     'Approx_Carson_wire_over_lossy_ground.out' u 1:3 title 'Im{gamma} (Approx Carson)' w l,\
     'Wait_wire_over_lossy_ground_all_modes.out' u 1:4 title 'Im{gamman} (Wait)' w p ,\
     'Wait_wire_over_lossy_ground.out' u 1:3 title 'Im{gamma1} (Wait)' w p ,\
     'Wait_wire_over_lossy_ground.out' u 1:6 title 'Im{gamma2} (Wait)' w p
pause -1


set title "attenuation over 1km"
set logscale x
set ylabel 'dB'
set autoscale x
set autoscale y
#set yrange[-150:0.0]
plot 'Full_Carson_wire_over_lossy_ground.out' u 1:4 title 'attenuation over 1km (Full Carson)' w l,\
     'Approx_Carson_wire_over_lossy_ground.out' u 1:4 title 'attenuation over 1km (Approx Carson)' w l,\
     'Wait_wire_over_lossy_ground_all_modes.out' u 1:4 title 'attenuation over 1km (Wait mode n)' w p,\
     'Wait_wire_over_lossy_ground.out' u 1:4 title 'attenuation over 1km (Wait mode 1)' w p,\
     'Wait_wire_over_lossy_ground.out' u 1:7 title 'attenuation over 1km (Wait mode 2)' w p

pause -1


