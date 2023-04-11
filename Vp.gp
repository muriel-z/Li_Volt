#Grafico del problema 3b

set xlabel 'Particula'
set ylabel ' Valor del voltaje'

set xzeroaxis


set title ' '

set key box top left
#set logscale y

#set xrange [-0.3:0.3]
#set yrange [0.0001:1000]


set grid xtic
set grid ytic

plot 'Vp.dat' u 1:5 t 'Vp' w p pt 7 ps 0.6

set terminal pdf color
set output 'grafico_Vp.pdf'
replot 
exit
