
set term jpeg
set nologscale x
set autoscale x
set autoscale y

set output "input_signal.jpg"
set title "Input signal"
set xlabel "time (s)"
set ylabel "Voltage (V)"
plot "input_signal.dat" u 1:2 w l
#pause -1

set autoscale x
set autoscale y
set xlabel "time (s)"
set ylabel "Voltage (V)"
set output "sub_sampled_signals.jpg"
set title "sub-sampled signal 1"
plot "sub_segment_time_1.dat" u 1:2 w l,\
     "sub_segment_time_2.dat" u 1:2 w l,\
     "sub_segment_time_3.dat" u 1:2 w l,\
     "sub_segment_time_4.dat" u 1:2 w l
#pause -1

set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "dBuV"
set xrange [1000:2e6]
#set yrange [50:120]
set output "dBuV_lin.jpg"
plot "V.fft" u 1:6 w l,\
     "V.favg" u 1:4 w l,\
     "square_wave_fourier_series.fout" u 1:4 w p
#pause -1


set autoscale x
set autoscale y
set xlabel "Frequency (Hz)"
set ylabel "dBuV"
set xrange [1000:2e6]
#set yrange [50:120]
set output "dBuV_log.jpg"
set logscale x
plot "V.fft" u 1:6 w l,\
     "V.favg" u 1:4 w l,\
     "square_wave_fourier_series.fout" u 1:4 w p
#pause -1


