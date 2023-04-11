#
set title " "
#set xrange [0:2000]
#set yrange [0:300000]

set style data histogram
set style histogram rowstacked
#set style fill solid border -1
#set boxwidth 10
#set xtic rotate by -45 scale 0
#set bmargin 10 
plot 'Vp.dat' using 5:xtic(1)

set terminal pdf color
set output 'historial_Vp.pdf'
replot 
exit
