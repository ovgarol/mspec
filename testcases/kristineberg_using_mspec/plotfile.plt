reset
set terminal pdfcairo font ',7' size 6,4
set output 'output.pdf'
set multiplot layout 3,3
SCALE = 0.5

set style line 2  lc rgb '#0025ad' lt 1 lw 1.5 # --- blue
set style line 3  lc rgb '#0042ad' lt 1 lw 1.5 #      .
set style line 4  lc rgb '#0060ad' lt 1 lw 1.5 #      .
set style line 5  lc rgb '#007cad' lt 1 lw 1.5 #      .
set style line 6  lc rgb '#0099ad' lt 1 lw 1.5 #      .
set style line 7  lc rgb '#00ada4' lt 1 lw 1.5 #      .
set style line 8  lc rgb '#00ad88' lt 1 lw 1.5 #      .
set style line 9  lc rgb '#00ad6b' lt 1 lw 1.5 #      .
set style line 10 lc rgb '#00ad4e' lt 1 lw 1.5 #      .
set style line 11 lc rgb '#00ad31' lt 1 lw 1.5 #      .
set style line 12 lc rgb '#00ad14' lt 1 lw 1.5 #      .
set style line 13 lc rgb '#09ad00' lt 1 lw 1.5 # --- green
set style line 14 lc rgb '#08ad00' lt 1 lw 1.5 # ---?

set xlabel 'day of experiment'
set key top right
plot 'outputlow.dat'  u ($0*SCALE):54  w l t 'low CO2' lc rgb 'blue' lw 1.5,'observations/lowCO2.dat' u ($0*2):3 w p pt 6 ps 0.5 lc rgb 'blue' notitle, 'outputhigh.dat'  u ($0*SCALE):54  w l t 'high CO2' lc rgb 'red' lw 1.5,'observations/highCO2.dat' u ($0*2):3 w p pt 6 ps 0.5 lc rgb 'red' notitle
set key top right
plot 'outputlow.dat' u ($0*SCALE):55 w l t 'low CO2' lc rgb 'blue' lw 2, 'observations/lowCO2.dat' u ($0*2):4 w p pt 6 ps 0.5 lc rgb 'blue' notitle,'outputhigh.dat' u ($0*SCALE):55 w l t 'high CO2' lc rgb 'red' lw 2,'observations/highCO2.dat' u ($0*2):4 w p pt 6 ps 0.5 lc rgb 'red' notitle
unset y2tics
unset y2label
plot 'outputlow.dat' u ($0*SCALE):58 w l lw 2 lc rgb 'blue' t 'low CO2', 'observations/lowCO2.dat' u ($0*2):7 w p pt 6 ps 0.5 lc rgb 'blue' notitle, 'outputhigh.dat' u ($0*SCALE):58 w l lw 2 lc rgb 'red' t 'high CO2', 'observations/highCO2.dat' u ($0*2):7 w p pt 6 ps 0.5 lc rgb 'red' notitle
plot 'outputlow.dat' u ($0*SCALE):62 w l lc rgb 'blue' lw 2 t 'low CO2','observations/lowCO2.dat' u ($0*2):9 w p pt 6 ps 0.5 lc rgb 'blue' notitle, 'outputhigh.dat' u ($0*SCALE):62 w l lc rgb 'red' lw 2 t 'high CO2','observations/highCO2.dat' u ($0*2):9 w p pt 6 ps 0.5 lc rgb 'red' notitle
set title 'low CO2'
set ylabel ''
set autoscale y

plot 'outputlow.dat' u ($0*SCALE):(0*$0):(-$171) w filledcu lc rgb 'web-green' title 'growth',\
 'outputlow.dat' u ($0*SCALE):(0*$0):172 w filledcu lc rgb 'khaki' title 'respiration',\
 'outputlow.dat' u ($0*SCALE):172:($172+$173) w filledcu lc rgb 'light-blue' title 'sinking',\
 'outputlow.dat' u ($0*SCALE):($172+$173):($172+$173+$174) w filledcu lc rgb 'dark-turquoise' title 'aggregation',\
 'outputlow.dat' u ($0*SCALE):($172+$173+$174):($172+$173+$174+$175) w filledcu lc rgb 'sienna1' title 'grazing',\


set title 'high CO2'
set autoscale y
plot 'outputhigh.dat' u ($0*SCALE):(0*$0):(-$171) w filledcu lc rgb 'web-green' title 'growth',\
 'outputhigh.dat' u ($0*SCALE):(0*$0):172 w filledcu lc rgb 'khaki' title 'respiration',\
 'outputhigh.dat' u ($0*SCALE):172:($172+$173) w filledcu lc rgb 'light-blue' title 'sinking',\
 'outputhigh.dat' u ($0*SCALE):($172+$173):($172+$173+$174) w filledcu lc rgb 'dark-turquoise' title 'aggregation',\
 'outputhigh.dat' u ($0*SCALE):($172+$173+$174):($172+$173+$174+$175) w filledcu lc rgb 'sienna1' title 'grazing',\


set title 'biomass (low CO2)'
unset key
set ylabel 'size (log_{ESD})'
set ylabel offset -1.5,0
set pm3d map
set logscale cb
set xtics out
set yrange [0:17]
set ytics ('-0.5' 0, '0.0' 1, '0.5' 2, '1.0' 3, '1.3' 4, '1.6' 5, '1.9' 6, '2.2' 7, '2.5' 8, '3.0' 9, '3.5' 10, '4.0' 11, '4.5' 12, '5.0' 13, '5.5' 14, '6.0' 15, '6.2' 16)
set ytics offset 0,0.5
splot 'outputlow.dat' matrix u ($2*SCALE):($1-175):3 every ::175 w pm3d
set title 'biomass (high CO2)'
unset key
set ylabel 'size (log_{ESD})'
set ylabel offset -1.5,0
set pm3d map
set yrange [0:17]
splot 'outputhigh.dat' matrix u ($2*SCALE):($1-175):3 every ::175 w pm3d
unset multiplot
exit
