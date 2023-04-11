set xlabel 'malla en z'
set ylabel 'malla en x '

set xzeroaxis


set title 'Valores del voltaje'

set key bottom left

set pm3d map


splot 'Vxz.dat' u 2:1:3

set terminal png
set output 'grafico_Vxz_surf.png'
replot 
exit
