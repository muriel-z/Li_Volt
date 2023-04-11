#Grafico del problema 1c

set xlabel 'Vp(bins)'
set ylabel ''

set xzeroaxis

set title 'Histograma Vtop = 10'

set key top right

set style data histograms

set grid xtic
set grid ytic

set logscale y


set xtics 

plot 'hist_Vpp10.00.d' u 1:2 t 'Voltaje de las particulas' w boxes 

set terminal pdf color
set output 'Histograma_vp.pdf'
replot 
exit
