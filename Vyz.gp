set xlabel 'malla en z'
set ylabel 'malla en y '

set xzeroaxis


set title 'Valores del voltaje'

set key bottom left

set pm3d map


splot 'Vyz.dat' u 2:1:3

set terminal png
set output 'grafico_Vyz.png'
replot 
exit
